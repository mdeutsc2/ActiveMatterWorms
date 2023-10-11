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
    var atoms: list(int); // list of particle id's in each bin
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
    }
}
var solvent: [1..numParticles] Particle;
var KE: [1..numParticles] real;
var KE_total: [1..nsteps/print_interval] real;
var amom: [1..nsteps/print_interval] real;
var xmom: [1..nsteps/print_interval] real;
var ymom: [1..nsteps/print_interval] real;

var numBins = ceil(L/rcut):int;
//const binSpace = {1..numBins, 1..numBins};
const binSpace = {1..numBins*numBins};
var bins : [1..numBins*numBins] Bin;
var binSpaceiodd : [1..(numBins*numBins)/2] int;
var binSpacejodd : [1..(numBins*numBins)/2] int;
var binSpaceieven : [1..(numBins*numBins)/2] int;
var binSpacejeven : [1..(numBins*numBins)/2] int;
var randStream = new RandomStream(real); // creating random number generator
const numTasks = here.numPUs();
if here.gpus.isEmpty() {
    writeln("no gpus");
} else {
    writeln("gpus",here.gpus);
}

proc main () {
    writeln(numTasks);
    writeln("numBins: ",numBins);
    init_bins();

    // init particles
    // init in lattice
    if (numParticles != sqrt(numParticles)*sqrt(numParticles)) {writeln("non-square numParticles");halt();}
    var row_length = sqrt(numParticles):int;
    var row,col,center,spacing :real;
    spacing = 1.0;
    center = hxo2 - spacing*(row_length/2);
    for i in 1..numParticles {
        row = i % row_length;
        col = ((i - row)/row_length)-1;
        solvent[i].x = center + spacing*rcutsmall*row;
        solvent[i].y = center + spacing*rcutsmall*col;
        solvent[i].vx = 0.1*gaussRand(0.0,1.0);
        solvent[i].vy = 0.1*gaussRand(0.0,1.0);
        solvent[i].ptype = 2;
        //KE[i] = 0.0;
        KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
    var vxave = (+ reduce solvent.vx)/numParticles;
    var vyave = (+ reduce solvent.vy)/numParticles;
    for i in 1..numParticles {
        solvent[i].vx = solvent[i].vx - vxave;
        solvent[i].vy = solvent[i].vy - vyave;
    }

    // init randomly


    update_cells(0);

    var t = 0.0;
    var total_time = 0.0;
    var ct: stopwatch, wt:stopwatch, xt:stopwatch; //calc time, io time, totaltime
    write_xyz(0);
    //calc_forces_old();
    calc_forces(here.gpus[0]);
    var vel_mag = 0.0;
    //setting up stopwatch
    xt.start();
    for istep in 1..nsteps {
        if (istep % print_interval == 0) {
	    xt.stop();
	    total_time = xt.elapsed();
            /*
	    amom[(istep/print_interval):int] = 0.0;
            xmom[(istep/print_interval):int] = 0.0;
            ymom[(istep/print_interval):int] = 0.0;
            var xcom,ycom,rx,ry: real;
            // calculating the center of mass
            xcom = 0.0;
            ycom = 0.0;
            for i in 1..numParticles {
                xcom = xcom + solvent[i].x;
                ycom = ycom + solvent[i].y;
            }
            xcom = xcom / numParticles;
            ycom = ycom / numParticles;
            // calculating the angular velocity
            for i in 1..numParticles {
                rx = solvent[i].x - xcom;
                ry = solvent[i].y - ycom;
                amom[(istep/print_interval):int] += (rx*solvent[i].vy - ry*solvent[i].vx);
                xmom[(istep/print_interval):int] += solvent[i].vx;
                ymom[(istep/print_interval):int] += solvent[i].vy;
            }
            amom[(istep/print_interval):int] = amom[(istep/print_interval):int]/numParticles;
            xmom[(istep/print_interval):int] = xmom[(istep/print_interval):int]/numParticles;
            ymom[(istep/print_interval):int] = ymom[(istep/print_interval):int]/numParticles;
            KE_total[(istep/print_interval):int] = (+ reduce KE);
	    */
            vel_mag = sqrt((max reduce solvent.vx)**2 + (max reduce solvent.vy)**2);
            writeln("Step: ",istep,"\t",
                    istep/total_time,"iter/s\tCalc:",
                    (ct.elapsed()/total_time)*100," %\tElapsed:",
                    total_time," s\t",
                    KE_total[(istep/print_interval):int],"\t",
                    (max reduce solvent.x),"\t",
                    (max reduce solvent.y),"\t",
                    vel_mag);
        xt.start();
        }
	    ct.start();
        // update positions
        update_position(here);
        update_cells(istep);

        //calc_forces_old();
        calc_forces(here.gpus[0]);

        update_velocities(here);


        ct.stop();
        if (istep % save_interval == 0){
            wt.start();
            write_xyz(istep);
            wt.stop();
        }
    }
    //write_macro(nsteps);
    xt.stop();
    writeln("Total Time:",xt.elapsed()," s");
}

proc update_position(loc) {
    on loc {
    // update positions
    forall i in 1..numParticles {
        solvent[i].x += solvent[i].vx*dt + (solvent[i].fx/2)*(dt**2);
        solvent[i].y += solvent[i].vy*dt + (solvent[i].fy/2)*(dt**2);
        //periodic boundary conditions
        // if solvent[i].x > L {
        //     solvent[i].x = solvent[i].x - L;
        // } else if solvent[i].x < 0.0 {
        //     solvent[i].x = solvent[i].x + L;
        // }
        // if solvent[i].y > L {
        //     solvent[i].y = solvent[i].y - L;
        // } else if solvent[i].y < 0.0 {
        //     solvent[i].y = solvent[i].y + L;
        // }

        solvent[i].vx += 0.5*dt*solvent[i].fx;
        solvent[i].vy += 0.5*dt*solvent[i].fy;
        // collide with boundaries
        if isclose(solvent[i].x, L, rtol = 1e-3, atol = 0.0) {
            solvent[i].vx = -1.0*solvent[i].vx;
        }
        if solvent[i].x <= 0.01 {
            solvent[i].vx = -1.0*solvent[i].vx;
        }
        if isclose(solvent[i].y, L, rtol = 1e-3, atol = 0.0) {
            solvent[i].vy = -1.0*solvent[i].vy;
        }
        if (solvent[i].y <= 0.01) {
            solvent[i].vy = -1.0*solvent[i].vy;
        }
        solvent[i].fx = 0.0;
        solvent[i].fy = 0.0;
    }
    }
}

proc lj(i,j) {
    var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,r,rhatx,rhaty,omega,fdissx,fdissy,gauss,frand :real;
    dx = solvent[j].x - solvent[i].x;
    dy = solvent[j].y - solvent[i].y;
    r2 = (dx*dx + dy*dy);
    //if (debug) {writeln("icount: ",icount,"\t",i,"\tjcount: ",jcount,"\t",j,"\t",r2,"\t",r2cut);}
    if (r2 <= r2cut) {
        // LJ force
        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        ffx = ffor*dx;
        ffy = ffor*dy;
        if (debug) {
            if (i == 1) || (j == 1) {
                writeln("LJ ",i," ",j," ",sqrt(r2)," ",sqrt(r2cut)," ",ffx," ",ffy," ",dx," ",dy);
            }
        }
        solvent[i].fx += ffx;
        solvent[i].fy += ffy;
        solvent[j].fx -= ffx;
        solvent[j].fy -= ffy;
        if (thermo) {
            // DPD thermostat
            //adding dissipative force
            dvx = solvent[j].vx - solvent[i].vx;
            dvy = solvent[j].vy - solvent[i].vy;
            r = sqrt(r2);
            rhatx = dx/r;
            rhaty = dy/r;
            omega = (1.0-r2/r2cut);
            fdissx = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhatx; //gamma = 1/damp (proportional to friction force)
            fdissy = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhaty;
            solvent[i].fx -= fdissx;
            solvent[i].fy -= fdissy;
            solvent[j].fx += fdissx;
            solvent[j].fy += fdissy;
            // adding random forces
            gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
            //gauss = randStream.getNext();
            frand = (1.0/sqrt(dt))*sqrt(omega)*gauss*sqrt(2.0*kbt*gamma);
            solvent[i].fx += frand*rhatx;
            solvent[i].fy += frand*rhaty;
            solvent[j].fx -= frand*rhatx;
            solvent[j].fy -= frand*rhaty;
        }
    }
}

proc calc_forces_old () {
    // if (thermo) {
    //     // needed for idiot thermostat
    //     var temp_KE_total = (+ reduce KE);
    //     alpha = sqrt(kbt*numParticles)/temp_KE_total;
    //     if (temp_KE_total <= 0.00001) {
    //         alpha = 1.0;
    //     }
    //     writeln(alpha);
    // }
    var Lo2 = L/2.0;
    for i in 1..numParticles {
        var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,r,frand,gauss,fdissx,fdissy,omega,alpha,lb,dp,col_prob:real;
        // calculating the force
        for j in (i+1)..numParticles {
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            // if (dx > Lo2) {
            //     dx = dx - L;
            // }
            // if (dx < -Lo2) {
            //     dx = dx + L;
            // }
            // if (dy > Lo2) {
            //     dy = dy - L;
            // }
            // if (dy < -Lo2) {
            //     dy = dy + L;
            // }
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                ffx = ffor*dx;
                ffy = ffor*dy;
                solvent[i].fx += ffx;
                solvent[i].fy += ffy;
                solvent[j].fx -= ffx;
                solvent[j].fy -= ffy;
                if (thermo) {
                    // DPD thermostat
                    //adding dissipative force
                    dvx = solvent[j].vx - solvent[i].vx;
                    dvy = solvent[j].vy - solvent[i].vy;
                    r = sqrt(r2);
                    rhatx = dx/r;
                    rhaty = dy/r;
                    omega = (1.0-r2/r2cut);
                    fdissx = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhatx; //gamma = 1/damp (proportional to friction force)
                    fdissy = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhaty;
                    solvent[i].fx -= fdissx;
                    solvent[i].fy -= fdissy;
                    solvent[j].fx += fdissx;
                    solvent[j].fy += fdissy;
                    // adding random forces
                    gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
                    //gauss = randStream.getNext();
                    frand = (1.0/sqrt(dt))*sqrt(omega)*gauss*sqrt(2.0*kbt*gamma);
                    solvent[i].fx += frand*rhatx;
                    solvent[i].fy += frand*rhaty;
                    solvent[j].fx -= frand*rhatx;
                    solvent[j].fy -= frand*rhaty;
                }
            }
        }
    }
}

proc calc_forces(loc) {
    // loop over all bins
    on loc {
    forall binid in binSpace {
        if (bins[binid].ncount > 1) {
            // calculate the forces between atoms inside each bin
            for icount in 0..bins[binid].ncount-2 {
                var i = bins[binid].atoms[icount];
                for jcount in (icount+1)..bins[binid].ncount-1 {
                    var j = bins[binid].atoms[jcount];
                    lj(i,j); // lennard-jones interaction between particles i and j;
                }
            }
        }
    }
    }
    // odd neighbors to the east (1)
    on loc {
    forall binid in binSpaceiodd {
        //var binid=(jbin-1)*numBins+ibin;
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj(i,j);
                    }

                }
            }
        }
    }
    }
    // even neighbors to the east (1)
    on loc {
    forall binid in binSpaceieven {
    // for ibin in 2..numBins by 2 {
    //     for jbin in 1..numBins {
            // var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[1];
            if (binidnbor != -1) {
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }

                    }
                }
            }
        //}
    }
    }
    // odd neighbors to the NE (2) (i+1,j+1)
    on loc {
    forall binid in binSpaceiodd {
    // for ibin in 1..numBins by 2 {
    //     for jbin in 1..numBins {
            // var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[2];
            if (binidnbor != -1) {
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }
                    }
                }
            }
        //}
    }
    }
    // even neighbors to the NE (2)
    on loc {
    forall binid in binSpaceieven {
    // for ibin in 2..numBins by 2 { // can we do this with every i, but every other j?
    //     for jbin in 1..numBins {
    //        var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[2];
            if (binidnbor != -1) { // check if neighbor is valid
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) { // check if neighbor has any atoms at all
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }

                    }
                }
            }
        //}
    }
    }
    // odd neighbors to the N (3)
    on loc {
    forall binid in binSpacejodd {
    // for ibin in 1..numBins {
    //     for jbin in 1..numBins by 2 {
    //        var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[3];
            if (binidnbor != -1) {
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }

                    }
                }
            }
        //}
    }
    }
    // even neighbors to the N (3)
    on loc {
    forall binid in binSpacejeven {
    // for ibin in 1..numBins {
    //     for jbin in 2..numBins by 2 {
    //         var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[3];
            if (binidnbor != -1) {
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }

                    }
                }
            }
        //}
    }
    }
    // odd neighbors to the NW (4)
    on loc {
    forall binid in binSpaceiodd {
    // for ibin in 1..numBins by 2 {
    //     for jbin in 1..numBins {
    //         var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[4];
            if (binidnbor != -1) {
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }

                    }
                }
            }
        //}
    }
    }
    // even neighbors to the NW (4)
    on loc {
    forall binid in binSpaceieven {
    // for ibin in 2..numBins by 2 {
    //     for jbin in 1..numBins {
    //         var binid=(jbin-1)*numBins+ibin;
            var binidnbor = bins[binid].neighbors[4];
            if (binidnbor != -1) {
                if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                    for icount in 0..bins[binid].ncount-1 {
                        var i = bins[binid].atoms[icount];
                        for jcount in 0..bins[binidnbor].ncount-1 {
                            var j = bins[binidnbor].atoms[jcount];
                            lj(i,j);
                        }

                    }
                }
            }
        //}
    }
    }
}

proc update_velocities(loc) {
    //update velocities
    on loc {
    forall i in 1..numParticles {
        //solvent[i].vx += 0.5*dt*(solvent[i].fxold + solvent[i].fx);
        //solvent[i].vy += 0.5*dt*(solvent[i].fyold + solvent[i].fy);
        //solvent[i].fxold = solvent[i].fx;
        //solvent[i].fyold = solvent[i].fy;
        solvent[i].vx += 0.5*dt*solvent[i].fx;
        solvent[i].vy += 0.5*dt*solvent[i].fy;
        // calculating kinetic energy here too
        KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
    }
}

// CELL LIST FUNCTIONS
proc update_cells(istep:int) {
    // ncount array set to zero
    // particle list set to empty
    for binid in binSpace{
        if bins[binid].ncount > 0{
            bins[binid].ncount = 0;
            bins[binid].atoms.clear();
        }
    }
    var ibin,jbin,binid :int;
    for i in 1..numParticles { // populating particles into bins
        ibin=ceil(solvent[i].x/rcut):int;
        jbin=ceil(solvent[i].y/rcut):int;
        binid=(jbin-1)*numBins+ibin;

        if (binid > numBins*numBins) {
            writeln(solvent[i].info());
            write_xyz(istep);
            dump_particles();
        } else if (binid < 1) {
            writeln(solvent[i].info());
            write_xyz(istep);
            dump_particles();
        }
        //append i to particle_list for binid
        bins[binid].ncount += 1;
        bins[binid].atoms.pushBack(solvent[i].id);
        //if (debug) {writeln(i,"\t",solvent[i].x,"\t",solvent[i].y,"\t",binposx,"\t",binposy);}


    }
    if (debug) {
        writeln("recalculating bins");
        for ibin in binSpace {
            if bins[ibin].ncount > 0 {
            writeln(bins[ibin].id,"\t",bins[ibin].atoms);
            }
        }
    }
}

proc init_bins() {
    // writeln((4/numBins):int +1); //row
    // writeln(4%numBins); // col
    var binid,ibinnab,jbinnab,binidnbor:int;
    for ibin in 1..numBins {
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
        // bins[ibin].x[0] = (ibin[0]-1)*rcut;
        // bins[ibin].x[1] = ibin[0]*rcut;
        // bins[ibin].y[0] = (ibin[1]-1)*rcut;
        // bins[ibin].y[1] = ibin[1]*rcut;

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
}

// I/O FUNCTIONS
proc write_xyz(istep:int) {
    var ic:int;
    var filename:string = "simple%{08u}.xyz".format(istep);
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
        var xyzfile = open(filename, ioMode.cw);
        var myFileWriter = xyzfile.writer();
        myFileWriter.writeln(numParticles + 4);
        myFileWriter.writeln("# 0");
        ic = 2;
        for i in 1..numParticles {
            myFileWriter.writeln("S ",solvent[i].x," ",solvent[i].y," ",0.0);
            ic += 1;
        }
        myFileWriter.writeln("E ",0.0," ",0.0," ",0.0);
        myFileWriter.writeln("E ",L," ",0.0," ",0.0);
        myFileWriter.writeln("E ",0.0," ",L," ",0.0);
        myFileWriter.writeln("E ",L," ",L," ",0.0);
        ic += 4;
        writeln(filename,"\t",ic," lines written");
    } catch e: Error {
        writeln(e);
    }
}

proc write_macro(nsteps: int) {
    var filename:string = "energies.dat";
    //var vel_mag = sqrt((max reduce solvent.vx)**2 + (max reduce solvent.vy)**2);
    try{
        var datfile = open(filename, ioMode.cw);
        var myFileWriter = datfile.writer();
        myFileWriter.writeln("step \t KE \t maxX \t maxY \t vel_mag");
        for istep in 1..nsteps/print_interval {
            myFileWriter.writeln(istep,"\t",
                                KE_total[istep],"\t",
                                amom[istep],"\t",
                                xmom[istep],"\t",
                                ymom[istep],"\t");
                                //(max reduce solvent.x),"\t",
                                //(max reduce solvent.y),"\t",
                                //vel_mag);
        }
        datfile.fsync();
        writeln("energies.dat written");
    } catch e: Error {
        writeln(e);
    }
}

proc dump_particles() {
    var filename:string = "particles.dump";
    try {
        var dfile = open(filename, ioMode.cw);
        var myFileWriter = dfile.writer();
        for i in 1..numParticles {
            myFileWriter.writeln(solvent[i].info());
        }
    } catch e: Error {
        writeln(e);
    }
}

// HELPER FUNCTIONS
inline proc gaussRand(mean: real, stddev: real): real {
    var u1 = randStream.getNext();
    var u2 = randStream.getNext();
    var z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);  // Box-Muller transform
    var gauss = mean + stddev * z0;
    if (gauss > 2.5) {
        gauss = 2.5;
    } else if (gauss < -2.5) {
        gauss = -2.5;
    }
    return gauss;
}
