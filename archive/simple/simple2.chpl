use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
use GPU;

config const debug = false,
            L = 100.0, // size of the simulation box in x and y
            nsteps = 50000,
            dt = 0.001,
            save_interval = 5000,
            numParticles = 625, // number of particles
            thermo = false,
            kbt = 0.5; //reduced temperature

var hxo2 = L/2,
    sigma = 1.0,
    rcut = 2.5*sigma,
    r2cut = rcut**2,
    rcutsmall = sigma*2.0**(1.0/6.0),
    r2cutsmall = rcutsmall*rcutsmall,
    pi = 4.0*atan(1.0),
    gamma = 3.0, //1.5
    print_interval = 100;

var ptc_init_counter = 1;

record Particle {
    var id: int;
    var x,y,z: real;
    var vx,vy,vz: real;
    var vxave,vyave,vzave: real;
    var fx,fy,fz: real;
    var fxold,fyold,fzold: real;
    var m: real; //mass
    var ptype: int; //  ptype=1 (active), ptype=2(solvent), ptype=3 (boundary), ptype=-1 (unassigned)
    proc init() {
        this.id = ptc_init_counter;
        ptc_init_counter += 1;
        this.x = 0.0;
        this.y = 0.0;
        this.vx = 0.0;
        this.vy = 0.0;
        this.vxave = 0.0;
        this.vyave = 0.0;
        this.vzave = 0.0;
        this.fx = 0.0;
        this.fy = 0.0;
        this.fxold = 0.0;
        this.fyold = 0.0;
        this.ptype = -1;
    }
    proc info() {
        var typestring: string;
        if (this.ptype == 1) {
            typestring = "active";
        } else if (this.ptype == 2) {
            typestring = "solvent";
        } else if (this.ptype == 3) {
            typestring = "boundary";
        } else {
            typestring = "unassigned";
        }
        var s:string = "id# %i \t type: %s \t mass: %s \n pos: %r \t %r \t %r \n vel: %r \t %r \t %r \n force: %r \t %r \t %r".format(this.id,typestring,this.m,this.x,this.y,this.vx,this.vy,this.fx,this.fy);
        return s;
    }
    proc p(px: real, py: real) {
        this.x = px;
        this.y = py;
    }
    proc p() {
        return (this.x,this.y,this.z);
    }
    proc v(velx: real, vely:real) {
        this.vx = velx;
        this.vy = vely;
    }
    proc v() {
        return (this.vx,this.vy);
    }
    proc f(forcex:real,forcey:real) {
        this.fx = forcex;
        this.fy = forcey;
    }
    proc f() {
        return (this.fx,this.fy);
    }
    proc set(p: Particle) {
        this.x = p.x;
        this.y = p.y;
        this.vx = p.vx;
        this.vy = p.vy;
        this.vxave = p.vxave;
        this.vyave = p.vyave;
        this.fx = p.fx;
        this.fy = p.fy;
        this.fxold = p.fxold;
        this.fyold = p.fyold;
    }
}

var bin_init_counter = 1;
record Bin {
    //var id: (int,int); //id of each bin
    var id: int; //id of each bin
    var atoms: [1..64] int; // list of particle id's in each bin (max 64)
    var neighbors: [1..4] int; // indices of each bin's neighboring bin
    var ncount: int; // count of number of particles in each bin
    var x: (real,real); // precalculate the max_x and min_x values for the space the box occupies
    var y: (real,real); // this is for easier neighbor list creation

    proc init() { // record initializer
        this.id = bin_init_counter;
        bin_init_counter += 1;
        this.ncount = 0;
        this.x = (0.0,0.0);
        this.y = (0.0,0.0);
        for i in 1..4 {
            this.neighbors[i] = -1;
        }
        for i in 1..64 {
          this.atoms[i] = -1;
        }
    }
}
if here.gpus.isEmpty() {
    writeln("no gpus");
} else {
    writeln("gpus",here.gpus);
}
const device = here.gpus[1];
writeln(device);

proc main() {
  on device {
    var randStream = new RandomStream(real); // creating random number generator

    // variables
    var solvent: [1..numParticles] Particle;
    var KE: [1..numParticles] real;
    var numBins = ceil(L/rcut):int;
    const binSpace = {1..numBins*numBins};
    var bins : [1..numBins*numBins] Bin;
    var binSpaceiodd : [1..(numBins*numBins)/2] int;
    var binSpacejodd : [1..(numBins*numBins)/2] int;
    var binSpaceieven : [1..(numBins*numBins)/2] int;
    var binSpacejeven : [1..(numBins*numBins)/2] int;

    // init bins
    // making lists of binids for different
    var icount = 0;
    for ibin in 1..numBins by 2 {
        for jbin in 1..numBins {
            icount += 1;
            binSpaceiodd[icount] = (jbin-1)*numBins+ibin;
        }
    }
    icount = 0;
    for ibin in 2..numBins by 2 {
        for jbin in 1..numBins {
            icount +=1;
            binSpaceieven[icount] = (jbin-1)*numBins+ibin;
        }
    }
    icount = 0;
    for ibin in 1..numBins {
        for jbin in 1..numBins by 2 {
            icount += 1;
            binSpacejodd[icount] = (jbin-1)*numBins+ibin;
        }
    }
    icount = 0;
    for ibin in 1..numBins {
        for jbin in 2..numBins by 2 {
            icount += 1;
            binSpacejeven[icount] = (jbin-1)*numBins+ibin;
        }
    }
    init_bins(bins,numBins);

    // initializing particles
    // init in lattice
    if (numParticles != sqrt(numParticles)*sqrt(numParticles)) {writeln("non-square numParticles");halt();}
    var row_length = sqrt(numParticles):int;
    var row,col,center,spacing :real;

    // generating gaussian random numbers;
    var grns : [1..numParticles] real;
    grns = gaussRand(numParticles,0.0,1.0);
    spacing = 1.0;
    center = hxo2 - spacing*(row_length/2);
    foreach i in 1..numParticles {
        row = i % row_length;
        col = ((i - row)/row_length)-1;
        solvent[i].x = center + spacing*rcutsmall*row;
        solvent[i].y = center + spacing*rcutsmall*col;
        solvent[i].vx = 0.1*grns[i];
        solvent[i].vy = 0.1*grns[i];
        solvent[i].ptype = 2;
        //KE[i] = 0.0;
        KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
    var vxave = (+ reduce solvent.vx)/numParticles;
    var vyave = (+ reduce solvent.vy)/numParticles;
    foreach i in 1..numParticles {
        solvent[i].vx = solvent[i].vx - vxave;
        solvent[i].vy = solvent[i].vy - vyave;
    }
    //
    update_cells(bins,solvent,numBins);
    calc_forces(solvent,bins,numBins,binSpaceiodd,binSpacejodd,binSpaceieven,binSpacejeven);
    var xt:stopwatch;
    var total_time :real;
    xt.start();
    for istep in 1..nsteps {
      if (istep % print_interval == 0) {
        xt.stop();
        total_time = xt.elapsed();
        writeln("Step: ",istep,"\t",
                istep/total_time," iter/s\tElapsed ",
                total_time," s");
        xt.start();
      }
      update_position(solvent);
      update_cells(bins,solvent,numBins);
      calc_forces(solvent,bins,numBins,binSpaceiodd,binSpacejodd,binSpaceieven,binSpacejeven);
      update_velocities(solvent,KE);
    }
  }
  writeln("Done!");
}

proc init_bins(bins,numBins:int) {
  var binid,ibinnab,jbinnab,binidnbor:int;
   foreach ibin in 1..numBins {
      for jbin in 1..numBins {
          binid = (jbin-1)*numBins+ibin;
          //  4   3   2
          //     *   1
          //
          /*  1 = i+1,j
              2 = i+1,j+1
              3 = i,j+1
              4 = i-1,j+1
          */
          ibinnab=ibin+1;
          jbinnab=jbin;
          if (ibinnab <= numBins) {
              binidnbor=(jbinnab-1)*numBins+ibinnab;
              bins[binid].neighbors[1] = binidnbor;
          } else {
              bins[binid].neighbors[1] = -1;
          }

          ibinnab=ibin+1;
          jbinnab=jbin+1;
          if (ibinnab <= numBins) && (jbinnab <= numBins) {
              binidnbor=(jbinnab-1)*numBins+ibinnab;
              bins[binid].neighbors[2] = binidnbor;
          } else {
              bins[binid].neighbors[2] = -1;
          }

          ibinnab=ibin;
          jbinnab=jbin+1;
          if (jbinnab <= numBins) {
              binidnbor=(jbinnab-1)*numBins+ibinnab;
              bins[binid].neighbors[3] = binidnbor;
          } else {
              bins[binid].neighbors[3] = -1;
          }

          ibinnab=ibin-1;
          jbinnab=jbin+1;
          if (ibinnab > 0) && (jbinnab <= numBins){
              binidnbor=(jbinnab-1)*numBins+ibinnab;
              bins[binid].neighbors[4] = binidnbor;
          } else {
              bins[binid].neighbors[4] = -1;
          }

      }
  }
  //return bins;
}

proc update_cells(bins,solvent,numBins:int) {
   foreach binid in 1..numBins*numBins {
    if bins[binid].ncount > 0 {
      // resetting the ids in bins
      for i in 1..bins[binid].ncount {
        bins[binid].atoms[i] = -1;
      }
      bins[binid].ncount = 0;
    }
  }
  // populdating particles into bins
  for i in 1..numParticles { // this must be run in serial
    var ibin,jbin,binid : int;
    ibin=ceil(solvent[i].x/rcut):int;
    jbin=ceil(solvent[i].y/rcut):int;
    binid=(jbin-1)*numBins+ibin;
    if (binid > numBins*numBins) {
      halt("binid > numBins*numBins");
    } else if (binid < 1) {
      halt("binid < 1");
    }
    bins[binid].ncount += 1;
    bins[binid].atoms[bins[binid].ncount] = solvent[i].id;
  }
  //return bins;
}

proc calc_forces(solvent,bins,numBins,binSpaceiodd,binSpacejodd,binSpaceieven,binSpacejeven) {
   foreach binid in 1..numBins*numBins {
    var i,j:int;
    var dx,dy,r2,ffor :real;
    if (bins[binid].ncount > 1) {
      for icount in 1..bins[binid].ncount -1 {
        i = bins[binid].atoms[icount];
        for jcount in (icount+1)..bins[binid].ncount {
          j = bins[binid].atoms[jcount];
          //if (i == -1) || (j == -1) {halt("1 i or j = -1 in atom list");}
          dx = solvent[j].x - solvent[i].x;
          dy = solvent[j].y - solvent[i].y;
          r2 = (dx*dx + dy*dy);
          if (r2 <= r2cut) {
              // LJ force
              ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
              solvent[i].fx += ffor*dx;
              solvent[i].fy += ffor*dy;
              solvent[j].fx -= ffor*dx;
              solvent[j].fy -= ffor*dy;
          }
        }
      }
    }
  }
  // odd neighbors to the east (1)
   foreach binid in binSpaceiodd {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[1];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("2 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }

  // even neighbors to the east (1)
   foreach binid in binSpaceieven {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[1];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("3 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }

  // odd neighbors to the NE (2) (i+1,j+1)
   foreach binid in binSpaceiodd {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[2];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("4 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }

  // even neighbors to the NE (2)
   foreach binid in binSpaceieven {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[2];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("5 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }

  // odd neighbors to the N (3)
   foreach binid in binSpacejodd {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[3];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("6 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }
  // even neighbors to the N (3)
   foreach binid in binSpacejeven {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[3];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("7 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }
  // odd neighbors to the NW (4)
   foreach binid in binSpaceiodd {
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[4];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("8 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }
  // even neighbors to the NW (4)
   foreach binid in binSpaceieven{
    var binidnbor,i,j :int; // create vars inside parfor
    var dx,dy,r2,ffor :real;
    binidnbor = bins[binid].neighbors[4];
    if (binidnbor != -1) { // test if neighbor is valid
      if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
        for icount  in 1..bins[binid].ncount {
          i = bins[binid].atoms[icount];
          for jcount in 1..bins[binid].ncount {
            j = bins[binidnbor].atoms[jcount];
            //if (i == -1) || (j == -1) {halt("9 i or j = -1 in atom list");}
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                solvent[i].fx += ffor*dx;
                solvent[i].fy += ffor*dy;
                solvent[j].fx -= ffor*dx;
                solvent[j].fy -= ffor*dy;
            }
          }
        }
      }
    }
  }
}

proc update_position(solvent) {
   foreach i in 1..numParticles {
    solvent[i].x += solvent[i].vx*dt + (solvent[i].fx/2)*(dt**2);
    solvent[i].y += solvent[i].vy*dt + (solvent[i].fy/2)*(dt**2);
    //periodic boundary conditions
    if solvent[i].x > L {
        solvent[i].x = solvent[i].x - L;
    } else if solvent[i].x < 0.0 {
        solvent[i].x = solvent[i].x + L;
    }
    if solvent[i].y > L {
        solvent[i].y = solvent[i].y - L;
    } else if solvent[i].y < 0.0 {
        solvent[i].y = solvent[i].y + L;
    }

    solvent[i].vx += 0.5*dt*solvent[i].fx;
    solvent[i].vy += 0.5*dt*solvent[i].fy;
    solvent[i].fx = 0.0;
    solvent[i].fy = 0.0;
  }
}

proc update_velocities(solvent,KE) {
   foreach i in 1..numParticles {
    solvent[i].vx += 0.5*dt*solvent[i].fx;
    solvent[i].vy += 0.5*dt*solvent[i].fy;
    // calculating kinetic energy here too
    KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
  }
}
// HELPER FUNCTIONS
proc gaussRand(n:int, mean: real, stddev: real) {
  var gauss : [1..n] real;
  var u1 : [1..n] real;
  var u2 : [1..n] real;
  fillRandom(u1);
  fillRandom(u2);
  // generating gaussian random numbers;
   foreach i in 1..numParticles {
    /* var z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);  // Box-Muller transform
    gauss[i] = mean + stddev * z0; */
    gauss[i] = mean + stddev * (sqrt(-2.0 * log(u1[i])) * cos(2.0 * pi * u2[i]));
    if (gauss[i] > 2.5) {
        gauss[i] = 2.5;
    } else if (gauss[i] < -2.5) {
        gauss[i] = -2.5;
    }
  }
  return gauss;
}
