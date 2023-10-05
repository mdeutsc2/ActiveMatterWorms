use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
// user-defined modules
//use Structs; //imports Particle and Bin structures
//use Worms;

const numTasks = here.numPUs();
// configuration
config const np = 40,
            nworms = 250,
            nsteps = 2000000    ,//00,
            fdogic = 0.06,
            walldrive = true;
            fdogicwall = 0.0,
            fdep = 1.0, // TODO: change to 4.0?
            fdepwall = 0.0,
            diss = 0.08,
            dt = 0.001, //0.02
            kspring = 57.146436,
            kbend = 40.0,
            length0 = 0.8, //particle spacing on worms
            rcut = 2.5,
            save_interval = 250,
            boundary = 1, // 1 = circle, 2 = cardioid, 3 = channel
            fluid_cpl = true,
            debug = false,
            thermo = true, // turn thermostat on?
            kbt = 0.5,
            sigma = 2.0;

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
        this.z = 0.0;
        this.vx = 0.0;
        this.vy = 0.0;
        this.vz = 0.0;
        this.vxave = 0.0;
        this.vyave = 0.0;
        this.vzave = 0.0;
        this.fx = 0.0;
        this.fy = 0.0;
        this.fz = 0.0;
        this.fxold = 0.0;
        this.fyold = 0.0;
        this.fzold = 0.0;
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
        var s:string = "id# %i \t type: %s \t mass: %s \n pos: %r \t %r \t %r \n vel: %r \t %r \t %r \n force: %r \t %r \t %r".format(this.id,typestring,this.m,this.x,this.y,this.z,this.vx,this.vy,this.vz,this.fx,this.fy,this.fz);
        return s;
    }
    proc p(px: real, py: real, pz: real) {
        this.x = px;
        this.y = py;
        this.z = pz;
    }
    proc p() {
        return (this.x,this.y,this.z);
    }
    proc v(velx: real, vely:real, velz:real) {
        this.vx = velx;
        this.vy = vely;
        this.vz = velz;
    }
    proc v() {
        return (this.vx,this.vy,this.vz);
    }
    proc f(forcex:real,forcey:real,forcez:real) {
        this.fx = forcex;
        this.fy = forcey;
        this.fz = forcez;
    }
    proc f() {
        return (this.fx,this.fy,this.fz);
    }
    proc set(p: Particle) {
        this.x = p.x;
        this.y = p.y;
        this.z = p.z;
        this.vx = p.vx;
        this.vy = p.vy;
        this.vz = p.vz;
        this.vxave = p.vxave;
        this.vyave = p.vyave;
        this.vzave = p.vzave;
        this.fx = p.fx;
        this.fy = p.fy;
        this.fz = p.fz;
        this.fxold = p.fxold;
        this.fyold = p.fyold;
        this.fzold = p.fzold;
    }
}

var bin_init_counter = 1;
record Bin {
    //var id: (int,int); //id of each bin
    var id: int; //id of each bin
    var atoms: list(int); // list of particle id's in each bin (wormid = (iw-1)*np+ip;)
    var types: list(int); // corresponding particle type in each bin
    var neighbors: [1..4] int; // indices of each bin's neighboring bin
    var ncount: int; // count of number of particles in each bin
    var wcount: int; // count of number of worm particles in each bin
    var scount: int; // count of number of solvent particles in each bin
    var bcount: int; // count of number of boundary particles in each bin
    var x: (real,real); // precalculate the max_x and min_x values for the space the box occupies
    var y: (real,real); // this is for easier neighbor list creation

    proc init() { // record initializer
        this.id = bin_init_counter;
        bin_init_counter += 1;
        this.ncount = 0;
        this.wcount = 0;
        this.scount = 0;
        this.bcount = 0;
        this.x = (0.0,0.0);
        this.y = (0.0,0.0);
        for i in 1..4 {
            this.neighbors[i] = -1;
        }
    }
}


// variables
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1.0/6.0),
      r2cutsmall = rcutsmall*rcutsmall,
      rwall = 74,//125.0*rcutsmall*sqrt(2.0),
      pi = 4.0*atan(1.0),
      twopi = 2*pi,
      pio4 = pi*0.25,
      density = nworms*np/(pi*rwall**2), //density for a circle
      hx = 2.0*rwall + 1.0,
      hy = hx,
      hyo2 = hy/2,
      hxo2 = hx/2,
      nxcell = ((hx/rcutsmall):int - 1),
      nycell = ((hy/rcutsmall):int - 1),
      dcell = hx/(nxcell:real),
      ncells = nxcell*nycell,
      nstepso500 = nsteps/500,
      gnoise = 0.80/sqrt(10.0)*0.8,
      dt2o2 = dt*dt*0.50,
      dto2 = dt*0.50,
      length2 = 2.0*length0,
      lengthmax = (length0)*((np - 1):real),
      r2inside = (rwall - rcutsmall) * (rwall-rcutsmall),
      a = 0.24, // layer spacing of worms in init_worms?
      gamma = 3.0, // frictional constant for dissipative force (~1/damp)
      numPoints = 420, //number of boundary points (for circle w/ r-75)
      numSol = 3200, // number of solution particles (3200 for circular, 800 for cardioid)
      fluid_offset = r2cutsmall-0.1;//3.0; // z-offset of fluid




var wormsDomain: domain(2) = {1..nworms,1..np};
var worms: [wormsDomain] Particle;
ptc_init_counter = 1;
var solvent: [1..numSol] Particle;
ptc_init_counter = 1;
var bound: [1..numPoints] Particle;
var savex: [1..np] real(64);
var savey: [1..np] real(64);
var ireverse: [1..nworms] int;
var ddx: [1..9] int;
var ddy: [1..9] int;
var hhead: [1..ncells] int; //
var ipointto: [1..nworms*np+numPoints] int; // linked list, every particle points to another particle
var nnab: [wormsDomain] int;
var KEworm: [1..nworms*np] real;
var KEsol: [1..numSol] real;
var KEworm_total: [1..nsteps] real;
var KEsol_total: [1..nsteps] real; //workaround https://stackoverflow.com/questions/59753273/how-to-append-data-to-an-existing-file

var numBins = ceil(hx/rcut):int;
writeln("numBins:\t",numBins);
//const binSpace = {1..numBins, 1..numBins};
const binSpace = {1..numBins*numBins};
var bins : [1..numBins*numBins] Bin;
var binSpaceiodd : [1..(numBins*numBins)/2] int;
var binSpacejodd : [1..(numBins*numBins)/2] int;
var binSpaceieven : [1..(numBins*numBins)/2] int;
var binSpacejeven : [1..(numBins*numBins)/2] int;

var randStream = new RandomStream(real); // creating random number generator

var t = 0.0;
var total_time = 0.0;
var ct: stopwatch, wt:stopwatch, xt:stopwatch; //calc time, io time, totaltime
//main
proc main() {
    writeln("starting...",numTasks);
    // save params to file
    write_params(); 
    // initialize the alternating loops over bins
    init_binspace();
    // initialize bins and neighboring bin lists
    init_bins();
    // initialize worms
    init_worms();

    if (fluid_cpl) {
        init_fluid();
    }

    // populate the bins with lists of atoms
    update_cells(0);
    // equilibrate the fluid
    if (fluid_cpl) {
    for istep in 1..5000 {
        var ioper = 5000/10;
        if (istep%ioper == 0) {
            writeln("fluid equilibration...",istep);
        }
        fluid_step(0,dt);
    }
    }
    writeln("fluid equilibrated...5000dt");
    write_xyzv(0);
    //setting up stopwatch
    xt.start();
    for itime in 1..nsteps {
        //writeln(itime);
        t = (itime:real) *dt;

        if (itime % 100 == 0) {
            xt.stop();
            total_time = xt.elapsed();
            writeln("Step: ",itime,"\t",itime/total_time,"iter/s\tCalc:",(ct.elapsed()/total_time)*100,"%\tIO:",(wt.elapsed()/total_time)*100," %\tElapsed:",total_time," s\tEst:",(nsteps/itime)*total_time," s");
            xt.start();
            }
	    ct.start();
        // first update positions and store old forces
        update_pos(itime);

        intraworm_forces();

        calc_forces(itime,dt);

        if (fluid_cpl) {
            KEsol_total[itime] = (+ reduce KEsol);
        } else {
            KEsol_total[itime] = 0.0;
        }

        update_vel();
        KEworm_total[itime] = (+ reduce KEworm);

	    ct.stop();
        if (itime % save_interval == 0){
            wt.start();
            write_xyzv(itime);
            wt.stop();
        }
    }
    //finalize();
    xt.stop();
    write_macro(nsteps);
    writeln("Total Time:",xt.elapsed()," s");
}


// //functions
proc init_worms() {
    //array for locating neighbor cells
    var rand1:real; // temp fix for randSteam
    var r:real, dth:real, xangle:real;
    ddx[1] = 1;
    ddy[1] = 0;
    ddx[2] = 1;
    ddy[2] = 1;
    ddx[3] = 0;
    ddy[3] = 1;
    ddx[4] = -1;
    ddy[4] = 1;
    ddx[5] = -1;
    ddy[5] = 0;
    ddx[6] = -1;
    ddy[6] = -1;
    ddx[7] = 0;
    ddy[7] = -1;
    ddx[8] = 1;
    ddy[8] = -1;
    ddx[9] = 0;
    ddy[9] = 0;
    var thetanow = 5.0*pi :real;
    var rmin = a*thetanow;
    writeln("nworms\t", nworms);
    writeln("np\t", np);
    writeln("rwall\t", rwall);
    writeln("nxcells\t", nxcell, "\t nycells\t", nycell);
    writeln("density\t", density);
    writeln("dt\t",dt);
    //setting up the worms
    if (boundary == 1){
        // circular bound
        var equidistantThetaValues: [1..numPoints] real = 0.0;
        var deltaTheta = 2 * pi / numPoints;

        for i in 1..numPoints {
            equidistantThetaValues[i] = i * deltaTheta;
        }

        for i in 1..numPoints {
            bound[i].x = rwall * cos(equidistantThetaValues[i])+hxo2;
            bound[i].y = rwall * sin(equidistantThetaValues[i])+hyo2;
            bound[i].z = 0.0;
            bound[i].ptype = 3;
        }
        for iw in 1..nworms {
            ireverse[iw] = 0;
            rand1 = randStream.getNext();
            if (rand1 <= 0.5) {
                ireverse[iw] = 1;
            }
            for i in 1..np {
                r = a*thetanow;
                dth = length0/r;
                thetanow += dth;
                worms[iw,i].x = hxo2 + r*cos(thetanow);
                worms[iw,i].y = hyo2 + r*sin(thetanow);
                xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                //TODO give them an initial velocity going around the circle
                worms[iw,i].ptype = 1;
                // vx[iw,i] = 0.0; // NOTE: this is all set by initializer
                // vy[iw,i] = 0.0;
                // vxave[iw,i] = 0.0;
                // vyave[iw,i] = 0.0;
                // fx[iw,i] = 0.0;
                // fy[iw,i] = 0.0;
                // fxold[iw,i] = 0.0;
                // fyold[iw,i] = 0.0;
            }
            thetanow += 4.0*dth;
        }
    } else if (boundary == 2) {
        // cardioid boundary
        var equidistantArcLengths: [1..numPoints] real;
        var thetaValues: [1..numPoints] real;
        var ca = 1.5*(rwall/2);
        var totalArcLength = 8 * ca;
        for i in 1..numPoints {
            equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
            thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
        }

        for i in 1..numPoints {
            bound[i].x = ca * (1 - cos(thetaValues[i])) * cos(thetaValues[i]) + hxo2 + ca;
            bound[i].y = ca * (1 - cos(thetaValues[i])) * sin(thetaValues[i]) + hyo2;
        }
        //now place worm particles
        for iw in 1..nworms {
            ireverse[iw] = 0;
            rand1 = randStream.getNext();
            if (rand1 <= 0.5) {
                ireverse[iw] = 1;
            }
            for i in 1..np {
                r = (0.65*a)*thetanow;
                dth = length0/r;
                thetanow += dth;
                worms[iw,i].x = hxo2 + r*cos(thetanow);
                worms[iw,i].y = hyo2 + r*sin(thetanow);
                xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                //TODO give them an initial velocity going around the circle
                worms[iw,i].ptype = 1;
                // vx[iw,i] = 0.0;
                // vy[iw,i] = 0.0;
                // vxave[iw,i] = 0.0;
                // vyave[iw,i] = 0.0;
                // fx[iw,i] = 0.0;
                // fy[iw,i] = 0.0;
                // fxold[iw,i] = 0.0;
                // fyold[iw,i] = 0.0;
            }
            thetanow += 4.0*dth;
        }

    } else if (boundary == 3) {
        // channel boundary
        writeln("channel boundary not implemented");
        halt();
    }
    //reverse some of the worms and give them crazy colors
    for iw in 1..nworms {
        if (ireverse[iw] == 1) {
            for i in 1..np {
                savex[i] = worms[iw,i].x;
                savey[i] = worms[iw,i].y;
            }
            for ip in 1..np {
                worms[iw,ip].x = savex[np+1-ip];
                worms[iw,ip].y = savey[np+1-ip];
            }
        }
    }

    //init macro variable arrays
    for i in 1..nworms*np{
        KEworm[i] = 0.0;
    }
    writeln("init done");
}

proc update_pos(itime:int) {
    //forall iw in 1..nworms {
    forall iw in 1..nworms {
       foreach i in 1..np {
            worms[iw,i].x = worms[iw,i].x + worms[iw,i].vx*dt + worms[iw, i].fx*dt2o2;
            worms[iw,i].y = worms[iw,i].y + worms[iw,i].vy*dt + worms[iw, i].fy*dt2o2;
            worms[iw,i].fxold = worms[iw,i].fx;
            worms[iw,i].fyold = worms[iw,i].fy;
            worms[iw,i].fx = 0.0;
            worms[iw,i].fy = 0.0;
            worms[iw,i].fz = 0.0;
        }
    }
    fluid_pos(dt);
}

proc intraworm_forces() {
    var rsq:real,rand1:real,rand2:real,v1:real,v2:real,fac:real,g1:real,th:real;
    //zero out the force arrays and add Gaussian noise
    rsq = 0.0;
    /* for iw in 1..nworms {
        for i in 1..np {
            while ((rsq >= 0.999) || (rsq <= 0.001)) {
                rand1 = randStream.getNext();
                rand2 = randStream.getNext();
                v1 = 2.0*rand1 - 1.0;
                v2 = 2.0*rand2 - 1.0;
                rsq = v1*v1 + v2*v2;
            }
            rand1 = randStream.getNext();
            fac = sqrt(-2.0*log(rsq)/rsq);
            g1 = v1*fac*gnoise;
            th = rand1*twopi;
            worms[iw,i].fx = g1*cos(th) - diss*worms[iw,i].vx;
            worms[iw,i].fy = g1*sin(th) - diss*worms[iw,i].vy;
        }
    } */
    //first set of springs nearest neighbor springs
    forall iw in 1..nworms {
        var ip1:int,r:real,ff:real,ffx:real,ffy:real,dx:real,dy:real;
        for i in 1..np-1 {
            ip1 = i + 1;
            dx = worms[iw,ip1].x - worms[iw,i].x;
            dy = worms[iw,ip1].y - worms[iw,i].y;
            r = sqrt(dx*dx + dy*dy);

            ff = -kspring*(r - length0)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            worms[iw,ip1].fx += ffx;
            worms[iw,i].fx -= ffx;
            worms[iw,ip1].fy += ffy;
            worms[iw,i].fy -= ffy;
        }
    }
    //bond bending terms
    forall iw in 1..nworms {
        var i3:int,i4:int,x2:real,x3:real,x4:real,y2:real,y3:real,y4:real,y23:real,y34:real,x23:real,x34:real,r23:real,r34:real,cosvalue:real;
	    var sinvalue:real,ff:real,dot:real,fac:real,f2x:real,f2y:real,f3x:real,f3y:real,f4x:real,f4y:real;
        for i2 in 1..(np-2) {
            i3 = i2 + 1;
            i4 = i2 + 2;
            //print*, i2,i3,i4
            x2 = worms[iw,i2].x;
            y2 = worms[iw,i2].y;
            x3 = worms[iw, i3].x;
            y3 = worms[iw, i3].y;
            x4 = worms[iw, i4].x;
            y4 = worms[iw, i4].y;
            y23 = y3 - y2;
            y34 = y4 - y3;
            x23 = x3 - x2;
            x34 = x4 - x3;
            r23 = sqrt(x23*x23 + y23*y23);
            r34 = sqrt(x34*x34 + y34*y34);

            cosvalue = abs(x23*x34 + y23*y34)/(r23*r34);
            //cosvalue = (x23*x34 + y23*y34)/(r23*r34)
            //print*, x23,x34,y23,y34,r23,r34
            if (cosvalue < 0.0) {
               cosvalue = 0.0;
               //print *, x23, x34, y23, y34, r23, r34
               //print *, "cosvalue .lt. 0.0d0"
               //stop
            }
            if (cosvalue > 1.0) {
               cosvalue = 1.0;
            }
            sinvalue = sqrt(1.0 - cosvalue*cosvalue);

            ff = -kbend*sinvalue/(r23*r34);

            dot = x23*x34 + y23*y34;
            fac = dot/(r23*r23);

            f2x = ff*(x34 - fac*x23);
            f2y = ff*(y34 - fac*y23);

            fac = dot/(r34*r34);
            f4x = ff*(fac*x34 - x23);
            f4y = ff*(fac*y34 - y23);
            f3x = -f2x - f4x;
            f3y = -f2y - f4y;

            worms[iw, i2].fx += f2x;
            worms[iw, i2].fy += f2y;

            worms[iw, i3].fx += f3x;
            worms[iw, i3].fy += f3y;

            worms[iw, i4].fx += f4x;
            worms[iw, i4].fy += f4y;
        }
    }
}

proc calc_forces() {
    forall binid in binSpace {
        if (bins[binid].ncount > 1) {
            // calculate the forces between atoms inside each bin
            for icount in 0..bins[binid].ncount-2 {
                var i = bins[binid].atoms[icount];
                var itype = bins[binid].types[icount];
                for jcount in (icount+1)..bins[binid].ncount-1 {
                    var j = bins[binid].atoms[jcount];
                    var jtype = bins[binidnbor].types[jcount];
                    cell_forces(i,j,itype,jtype);
                }
            }
        }
    }
    // odd neighbors to the east (1)
    forall binid in binSpaceiodd {
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }
    // even neighbors to the east (1)
    forall binid in binSpaceieven {
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != 1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // odd neighbors to the NE (2) (i+1,j+1)
    forall binid in binSpaceiodd {
        var binidnbor = bins[binid].neighbors[2];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // even neighbors to the NE (2)
    forall binid in binSpaceieven {
        var binidnbor = bins[binid].neighbors[2];
        if (binidnbor != -1) { // check if neighbor is valid
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) { // check if neighbor has any atoms at all
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // odd neighbors to the N (3)
    forall binid in binSpacejodd {
        var binidnbor = bins[binid].neighbors[3];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // even neighbors to the N (3)
    forall binid in binSpacejeven {
        var binidnbor = bins[binid].neighbors[3];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }

                }
            }
        }
    }

    // odd neighbors to the NW (4)
    forall binid in binSpaceiodd {
        var binidnbor = bins[binid].neighbors[4];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // even neighbors to the NW (4)
    forall binid in binSpaceieven {
        var binidnbor = bins[binid].neighbors[4];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }
}

proc cell_forces(i:int,j:int,itype:int,jtype:int) {
    //worm-worm interaction
    if (itype == 1) && (jtype == 1) {
        var iworm,jworm,ip,jp,inogo : int;
        iworm = 1 + ((i - 1)/np):int;
        ip = i - np*(iworm - 1);
        jworm = 1 + ((j - 1)/np):int;
        jp = j - np*(jworm - 1);

        inogo = 0;
        if ((iworm == jworm) && (abs(ip-jp) <= 2)) {
            // on the same worm and close means no interaction calculated here
            inogo = 1;
        }
        if (inogo == 0) {
            dddx = worms[jworm, jp].x - worms[iworm, ip].x;
            dddy = worms[jworm, jp].y - worms[iworm, ip].y;
            r2 = dddx**2 + dddy**2;
            riijj = sqrt(r2);
            //add attractive force fdep between all pairs
            if (r2 <= r2cutsmall) {
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
                ffx = ffor*dddx;
                ffy = ffor*dddy;
                worms[iworm,ip].fx += ffx;
                worms[jworm,jp].fx -= ffx;
                worms[iworm,ip].fy += ffy;
                worms[jworm,jp].fy -= ffy;

                //add 'dogic drive' to interacting pairs
                //first calculate unit vectors along each worm
                ip1 = ip + 1;
                if (ip1 <= np) {
                    dxi = worms[iworm,ip1].x - worms[iworm, ip].x;
                    dyi = worms[iworm,ip1].y - worms[iworm, ip].y;
                } else {
                    dxi = worms[iworm, ip].x - worms[iworm, ip - 1].x;
                    dyi = worms[iworm, ip].y - worms[iworm, ip - 1].y;
                }

                jp1 = jp + 1;
                if (jp1 <= np) {
                    dxj = worms[jworm, jp1].x - worms[jworm, jp].x;
                    dyj = worms[jworm, jp1].y - worms[jworm, jp].x;
                } else {
                    dxj = worms[jworm, jp].x - worms[jworm, jp - 1].x;
                    dyj = worms[jworm, jp].y - worms[jworm, jp - 1].y;
                }


                //if the two vectors have any component pointing in opposite directions
                if (dxi*dxj + dyi*dyj <= 0.0) {
                    //normalize those vectors to make them unit vectors
                    ri = sqrt(dxi*dxi + dyi*dyi);
                    dxi = dxi/ri;
                    dyi = dyi/ri;

                    rj = sqrt(dxj*dxj + dyj*dyj);
                    dxj = dxj/rj;
                    dyj = dyj/rj;
                    //now they are both unit vectors. Find the direction for the force...

                    dx = (dxi - dxj)/2.0;
                    dy = (dyi - dyj)/2.0;

                    //normalize

                    r = sqrt(dx*dx + dy*dy);
                    dx = dx/r;
                    dy = dy/r;

                    //add an extra attractive component where kinesin drive is present

                    ffx = fdogic*(dx) + 0.7*dddx/riijj;
                    ffy = fdogic*(dy) + 0.7*dddy/riijj;

                    //ffx=fdogic*(dx)
                    //ffy=fdogic*(dy)

                    worms[iworm,ip].fx += ffx;
                    worms[jworm,jp].fx -= ffx;
                    worms[iworm,ip].fy += ffy;
                    worms[jworm,jp].fy -= ffy;
                }
            }
        }
    }
    //worm-boundary interaction
    if ((itype == 1) && (jtype == 3)) {
        var iw,ip:int;
        iw = 1 + ((i - 1)/np):int; // find which worm j is in
        ip = i - np*(iw - 1); // which particle in the worm is j?
        dogic_wall(iw,ip,j);
    } else if ((itype == 3) && (jtype == 1)) {
        var iw,ip:int;
        iw = 1 + ((j - 1)/np):int; // find which worm j is in
        ip = j - np*(iw - 1); // which particle in the worm is j?
        dogic_wall(iw,ip,i);
    }
    //solvent-worm interaction
    if ((itype == 2) && (jtype == 1)) {
        var dx,dy,dz,r2,ffor,ffx,ffy:real;
        var iw,ip:int;
        iw = 1 + ((j - 1)/np):int; // find which worm j is in
        ip = j - np*(iworm - 1); // which particle in the worm is j?
        dz = fluid_offset;
        dx = solvent[i].x - worms[iw,ip].x;
        dy = solvent[i].y - worms[iw,ip].y;
        r2 = (dx*dx + dy*dy + dz*dz);
        if (r2 <= r2cutsmall) {
            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
            ffx = ffor*dx;
            ffy = ffor*dy;
            solvent[i].fx += ffx;
            solvent[i].fy += ffy;
            worms[iw,ip].fx -= ffx;
            worms[iw,ip].fy -= ffy;
        }
    } else if ((itype == 1) && (jtype == 2)) {
        var dx,dy,dz,r2,ffor,ffx,ffy:real;
        var iw,ip:int;
        iw = 1 + ((i - 1)/np):int; // find which worm j is in
        ip = i - np*(iworm - 1); // which particle in the worm is j?
        dz = fluid_offset;
        dx = solvent[j].x - worms[iw,ip].x;
        dy = solvent[j].y - worms[iw,ip].y;
        r2 = (dx*dx + dy*dy + dz*dz);
        if (r2 <= r2cutsmall) {
            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
            ffx = ffor*dx;
            ffy = ffor*dy;
            solvent[j].fx += ffx;
            solvent[j].fy += ffy;
            worms[iw,ip].fx -= ffx;
            worms[iw,ip].fy -= ffy;
        }
    }
    //solvent-solvent interaction
    if (itype == 2) && (jtype == 2) {
        lj_thermo(i,j,r2cut);
    }
    //solvent-boundary interaction
    if ((itype == 2) && (jtype == 3)) {
        lj(i,j,sigma*r2cutsmall);
    } else if ((itype == 3) && (jtype == 2)) {
        lj(j,i,sigma*r2cutsmall);
    }
}

proc dogic_wall(iw:int,ip:int,ib:int){
    var dx:real, dy:real, r:real, r2:real, th:real,
        xwall:real, ywall:real, rr2:real, ffor:real,
        dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
    var ip1,ib1:int;
    //calculate distance to the wall
    dx = worms[iw,i].x-bound[ib].x;
    dy = worms[iw,i].y-bound[ib].y;
    r2 = (dx*dx + dy*dy);
    //if close enough to the wall, calculate wall forces
    //use the short cut-off
    if (r2 <= r2cut) {
        ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        worms[iw,i].fx += ffor*dx;
        worms[iw,i].fy += ffor*dy;
        if (walldrive) {
            //first calculate unit vector along the worm
            ip1 = i + 1;
            if (ip1 <= np) {
                dxi = worms[iw, ip1].x - worms[iw, i].x;
                dyi = worms[iw, ip1].y - worms[iw, i].y;
            } else {
                dxi = worms[iw, i].x - worms[iw, i - 1].x;
                dyi = worms[iw, i].y - worms[iw, i - 1].y;
            }
            //make it a unit vector
            ri = sqrt(dxi*dxi + dyi*dyi);
            dxi = dxi/ri;
            dyi = dyi/ri;
            //calculate the unit vector along the wall
            ib1 = ib + 1;
            if (ib1 <= numPoints) {
                dxj = bound[ib1].x - bound[ib].x;
                dyj = bound[ib1].y - bound[ib].y;
            } else {
                dxj = bound[ib].x - bound[ib-1].x;
                dyj = bound[ib].y - bound[ib-1].y;
            }
            ri = sqrt(dxj*dxj + dyj*dyj);
            dxj = dxj/ri;
            dyj = dyj/ri;

            //if the vectors are not antiparallel, reverse the vector along the wall
            if (dxi*dxj + dyi*dyj > 0.0) {
                dxj = -dxj;
                dyj = -dyj;
            }
            //if the two vectors have any component pointing in opposite directions
            if (dxi*dxj + dyi*dyj < 0.0) {
                //Find the direction for the force...
                dx = (dxi - dxj)/2.0;
                dy = (dyi - dyj)/2.0;
                //normalize the direction vector
                ri = sqrt(dx*dx + dy*dy);
                dx = dx/ri;
                dy = dy/ri;

                //turn on extra-strong driving force
                ffx = fdogicwall*dx;
                ffy = fdogicwall*dy;
                worms[iw,i].fx += ffx;
                worms[iw,i].fy += ffy;
            }
        }
    }
}

proc update_vel() {
    forall iw in 1..nworms {
        foreach i in 1..np {
            worms[iw, i].vx += dto2*(worms[iw,i].fx + worms[iw,i].fxold);
            worms[iw, i].vy += dto2*(worms[iw,i].fy + worms[iw,i].fyold);
            //vxave[iw, i] = vxave[iw, i]/nnab[iw, i];
            //vyave[iw, i] = vyave[iw, i]/nnab[iw, i];
            KEworm[iw*i] = 0.5*(worms[iw,i].vx * worms[iw,i].vx + worms[iw,i].vy * worms[iw,i].vy);
        }
    }
    fluid_vel(dt);
}

//FLUID FUNCTIONS
proc init_fluid() {
    var random_placement = false;
    if (random_placement) {
        // put the solvent particles in a random x and random y
        if (boundary == 1) {
            // circular boundary
            var numTries = 1000;
            var rand1,rand2,alpha,x,y,r,r2,dx,dy : real;
            solvent[1].x = hxo2;
            solvent[1].y = hyo2; // TODO add solvent[1].z
            for i in 2..numSol {
                var valid = false;
                var try_count = 0;
                while (valid == false) {
                    rand1 = randStream.getNext();
                    rand2 = randStream.getNext();
                    // random angle
                    alpha = 2*pi*rand1;
                    // random radius
                    r = rwall*0.95 * sqrt(rand2);
                    // calculating coordinates
                    x = r * cos(alpha) + hxo2;
                    y = r * sin(alpha) + hyo2;
                    var invalid_count = 0;
                    for j in (i-1)..1 by -1 {
                        dx = abs(x - solvent[i].x);
                        dy = abs(y - solvent[i].y);
                        r2 = (dx*dx + dy*dy);
                        if (r2 <= r2cut) {
                            invalid_count += 1;
                        }
                    }
                    if (invalid_count == 0) {
                        // no conflicts found with any of the existing particles
                        valid = true;
                        //writeln(i,"\t",x,"\t",y);
                        solvent[i].x = x;
                        solvent[i].y = y;
                        //add small gaussian velocity here
                        break;
                    }
                    try_count += 1;
                    if (try_count > numTries+1) {
                        writeln(i,"\t",try_count,"\t",x,"\t",y);
                        writeln("particle placement failed, increase numTries");
                        halt();
                    }
                }
            }
        } else if (boundary == 2) {
            // cardioid boundary
            writeln("haven't put in cardioid fluid yet");
            halt();
        } else if (boundary == 3) {
            // channel boundary
            writeln("haven't put in channel fluid");
            halt();
        } else {
            writeln("something went wrong in solvent init. invalid boundary value");
            halt();
        }
    } else {
        // place particles in a square lattice centered in the boundary
        if (boundary == 1) {
            // circular boundary
            var fluid_a = sqrt(2)*(0.95*rwall); // box size of fluid
            var fluid_px = fluid_a + hxo2;
            var fluid_mx = fluid_a - hxo2;
            var fluid_py = fluid_a + hxo2;
            var fluid_my = fluid_a - hyo2;
            var spacing = 1.5;
            var row_length = numSol/(floor(fluid_a/(rcutsmall*spacing)):int);
            writeln("Row:",row_length,"\t",row_length**2,"\t",fluid_a,"\t",fluid_a/(rcutsmall*spacing));
            if (row_length**2 > numSol) {
                writeln("fluid density too high, fixme!");
                halt();
            }
            var row,col:real;
            for i in 1..numSol {
                row = i % row_length;
                col = ((i - row)/row_length)-1;
                solvent[i].x = fluid_mx + spacing*rcutsmall*row;
                solvent[i].y = fluid_my + spacing*rcutsmall*col;
                solvent[i].z = 0.0;
                solvent[i].vx = 0.0; //taken care of by type init
                solvent[i].vy = 0.0;
                solvent[i].vz = 0.0;
                //solfx[i] = 0.0;
                //solfy[i] = 0.0;
                //solfxold[i] = 0.0;
                //solfyold[i] = 0.0;
            }
        } else if (boundary == 2) {
            // cardioid boundary
            // get the fluid init bounds from the max x and y coords of the worms
            var max_x = 0.0;
            var max_y = 0.0;
            for iw in 1..nworms {
                for i in 1..np {
                    if (worms[iw,i].x > max_x) {
                        max_x = worms[iw,i].x;
                    }
                    if (worms[iw,i].y > max_y) {
                        max_y = worms[iw,i].y;
                    }
                }
            }
            // circular boundary
            var fluid_a = sqrt(2)*0.75*(max_x-hxo2); // box size of fluid
            var fluid_px = fluid_a + hxo2;
            var fluid_mx = fluid_a - hxo2;
            var fluid_py = fluid_a + hyo2;
            var fluid_my = fluid_a - hyo2;
            var spacing = 1.5;
            var row_length = numSol/(floor(fluid_a/(rcutsmall*spacing)):int);
            writeln("Row:",row_length,"\t",row_length**2,"\t",fluid_a,"\t",fluid_a/(rcutsmall*spacing));
            if (row_length**2 > numSol) {
                writeln("fluid density too high, fixme!");
                halt();
            }
            var row,col:real;
            for i in 1..numSol {
                row = i % row_length;
                col = ((i - row)/row_length)-1;
                solvent[i].x = hxo2-0.75*(max_x-hxo2) + spacing*rcutsmall*row;
                solvent[i].y = hyo2-0.75*(max_y-hyo2) + spacing*rcutsmall*col;
                solvent[i].z = 0.0;
                solvent[i].vx = 0.0; //taken care of by type init
                solvent[i].vy = 0.0;
                solvent[i].vz = 0.0;
                //solfx[i] = 0.0;
                //solfy[i] = 0.0;
                //solfxold[i] = 0.0;
                //solfyold[i] = 0.0;
            }
        } else if (boundary == 3) {
            writeln("haven't put in channel fluid");
            halt();
        } else {
            writeln("invalid boundary");
            halt();
        }
    }
    // init macro variable arrays
    for i in 1..numSol {
        KEsol[i] = 0.0;
        solvent[i].ptype = 2;
    }
    writeln("fluid init done");
}

proc lj_thermo(i:int,j:int,r2cut_local:real) {
    var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,r,rhatx,rhaty,omega,fdissx,fdissy,gauss,frand :real;
    dx = solvent[j].x - solvent[i].x;
    dy = solvent[j].y - solvent[i].y;
    r2 = (dx*dx + dy*dy);
    //if (debug) {writeln("icount: ",icount,"\t",i,"\tjcount: ",jcount,"\t",j,"\t",r2,"\t",r2cut);}
    if (r2 <= r2cut_local) {
        // LJ force
        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        ffx = ffor*dx;
        ffy = ffor*dy;
        if (debug) {
            if (i == 1) || (j == 1) {
                writeln("LJ ",i," ",j," ",sqrt(r2)," ",sqrt(r2cut_local)," ",ffx," ",ffy," ",dx," ",dy);
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
            omega = (1.0-r2/r2cut_local);
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

proc lj(i:int,j:int,r2cut_local:real) {
    var dx,dy,r2,ffor,ffx,ffy,sigma12,sigma6:real;
    //calculate distance to the wall
    dx = solvent[j].x - bound[i].x;
    dy = solvent[j].y - bound[i].y;
    r2 = (dx*dx + dy*dy);
    //if close enough to the wall, calculate wall forces
    //use the short cut-off
    if (r2 <= r2cut_local) {
        //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        sigma12 = sigma**12;
        sigma6 = sigma**6;
        ffor = 48.0*sigma12*r2**(-7.0) -24*sigma6*r2**(-4.0);
        solvent[j].fx += ffor*dx;
        solvent[j].fy += ffor*dy;
    }
}

proc fluid_pos(dt_fluid:real) {
    forall i in 1..numSol {
        solvent[i].x += solvent[i].vx*dt_fluid + (solvent[i].fx/2)*(dt_fluid**2);
        solvent[i].y += solvent[i].vy*dt_fluid + (solvent[i].fy/2)*(dt_fluid**2);
        // updating velocity and position
        // solvx[i] += solfx[i] * dt;
        // solvy[i] += solfy[i] * dt;
        // solx[i] += solvx[i] * dt;
        // soly[i] += solvy[i] * dt;

        // doing a velocity half-step here
        solvent[i].vx += 0.5*dt*solvent[i].fx;
        solvent[i].vy += 0.5*dt*solvent[i].fy;
        solvent[i].vz = 0.0; // just in case
        solvent[i].fx = 0.0;
        solvent[i].fy = 0.0;
        solvent[i].fz = 0.0; // just in case
    }
}

proc fluid_force_old() {
    //     // fluid-fluid
    //     var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,r,frand,gauss,fdissx,fdissy,omega:real;
    //     for i in 1..numSol {
    //         // calculating the force
    //         for j in (i+1)..numSol {
    //             lj_thermo(i,j,r2cut);
    //         }
    //     }
    //     //fluid-boundary
    //     // var r2cutsol = r2cut
    //     var r2cutsol = sigma*2.0**(1.0/6.0);
    //     forall i in 1..numSol {
    //         // calculate the force on the boundaries.
    //         for ib in 1..numPoints  {
    //             lj(ib,i,r2cutsol);

    //         }
    //     }

    //     // fluid-worms
    //     var dz = fluid_offset;
    //     for i in 1..numSol{
    //         var dx,dy,r2,ffor,ffx,ffy:real;
    //         for iw in 1..nworms {
    //             for ip in 1..np {
    //                 dx = solvent[i].x - worms[iw,i].x;
    //                 dy = solvent[i].y - worms[iw,i].y;
    //                 r2 = (dx*dx + dy*dy + dz*dz);
    //                 if (r2 <= r2cutsmall) {
    //                     ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
    //                     ffx = ffor*dx;
    //                     ffy = ffor*dy;
    //                     solvent[i].fx += ffx;
    //                     solvent[i].fy += ffy;
    //                     worms[iw,ip].fx -= ffx;
    //                     worms[iw,ip].fy -= ffy;
    //                 }
    //             }
    //         }
    //     }
}

proc fluid_vel(dt_fluid:real) {
    //update velocities
    forall i in 1..numSol {
        //solvent[i].vx += 0.5*dt_fluid*(solvent[i].fxold + solvent[i].fx);
        //solvent[i].vy += 0.5*dt_fluid*(solvent[i].fyold + solvent[i].fy);
        //solvent[i].fxold = solvent[i].fx;
        //solvent[i].fyold = solvent[i].fy;
        solvent[i].vx += 0.5*dt_fluid*solvent[i].fx;
        solvent[i].vy += 0.5*dt_fluid*solvent[i].fy;
        // calculating kinetic energy here too
        KEsol[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
}

proc fluid_step(istep:int,dt_fluid:real) {
    // update positions
    fluid_pos(dt_fluid);
    // update_fluid_cells(istep);
    fluid_force_old();

    fluid_vel(dt_fluid);
}

// CELL LIST FUNCTIONS
proc update_cells(istep:int) {
    // ncount array set to zero
    // particle list set to empty

    // resetting all other bins lists for worms and solvent
    for binid in binSpace{
        if bins[binid].ncount > 0 {
            bins[binid].ncount = 0; // number of boundary atoms should not change
            bins[binid].scount = 0;
            bins[binid].wcount = 0;
            bins[binid].bcount = 0;
            bins[binid].atoms.clear();
            bins[binid].types.clear();
        }
    }
    var ibin,jbin,binid :int;
    for i in 1..numSol { // populating solvent particles into bins
        ibin=ceil(solvent[i].x/rcut):int;
        jbin=ceil(solvent[i].y/rcut):int;
        binid=(jbin-1)*numBins+ibin;
        if (binid > numBins*numBins) {
            //writeln(solvent[i].info());
            sudden_halt(istep);
        } else if (binid < 1) {
            //writeln(solvent[i].info());
            sudden_halt(istep);
        }
        //append i to particle_list for binid
        bins[binid].scount += 1;
        bins[binid].ncount += 1;
        bins[binid].atoms.pushBack(solvent[i].id);
        bins[binid].types.pushBack(solvent[i].ptype);s
        //if (debug) {writeln(i,"\t",solvent[i].x,"\t",solvent[i].y,"\t",binposx,"\t",binposy);}
    }
    var wormid : int;
    for iw in 1..nworms { // populating worm particles into bins
        for ip in 1..np {
            ibin=ceil(worms[iw,ip].x/rcut):int;
            jbin=ceil(worms[iw,ip].y/rcut):int;
            binid=(jbin-1)*numBins+ibin;
            wormid = (iw-1)*np+ip;
            if (binid > numBins*numBins) {
                //writeln(solvent[i].info());
                sudden_halt(istep);
            } else if (binid < 1) {
                //writeln(solvent[i].info());
                sudden_halt(istep);
            }
            //append i to particle_list for binid
            bins[binid].wcount += 1;
            bins[binid].ncount += 1;
            bins[binid].atoms.pushBack(wormid);
            bins[binid].types.pushBack(worms[iw,ip].ptype);
            //if (debug) {writeln(i,"\t",solvent[i].x,"\t",solvent[i].y,"\t",binposx,"\t",binposy);}
        }
    }
    for ib in 1..numPoints {
        ibin=ceil(bound[ib].x/rcut):int;
        jbin=ceil(bound[ib].y/rcut):int;
        binid=(jbin-1)*numBins+ibin;
        //append i to particle_list for binid
        bins[binid].bcount += 1;
        bins[binid].ncount += 1;
        bins[binid].atoms.pushBack(bound[ib].id);
        bins[binid].types.pushBack(bound[ib].ptype);
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

}

proc init_binspace() {
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

// IO FUNCTIONS
proc write_xyz(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename:string = "amatter%{07u}.xyz".format(istep);
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer();
    // number of active particles + 4 edge-defining particles + boundary + solvent (optional)
    if (fluid_cpl) {
        myFileWriter.writeln(nworms * np + 4 + numPoints + numSol);
    } else {
        myFileWriter.writeln(nworms*np + 4 + numPoints);
    }
    myFileWriter.writeln("# 0");
    ic = 2;
    for iw in 1..nworms {
        dx = worms[iw,1].x - hxo2;
        dy = worms[iw,1].y - hyo2;
        xang = atan2(dy,dx);
        rx = -sin(xang);
        ry = cos(xang);
        dot = (worms[iw,1].x - worms[iw,np].x)*rx + (worms[iw,1].y - worms[iw,np].y)*ry;
        if (dot >= 0.0) {
            for i in 1..np {
                myFileWriter.writeln("A ",worms[iw,i].x," ",worms[iw,i].y," ", 0.0);
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",worms[iw,i].x," ",worms[iw,i].y," ", 0.0);
                //myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        }
    }
    for i in 1..numPoints {
        myFileWriter.writeln("I ",bound[i].x," ",bound[i].y," ",0.0);
        ic += 1;
    }
    if (fluid_cpl) {
        for i in 1..numSol {
            myFileWriter.writeln("S ",solvent[i].x," ",solvent[i].y," ",0.0);
            ic += 1;
        }
    }
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0);
    //ic += 4;
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    writeln(filename,"\t",ic," lines written");
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}

proc write_xyzv(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename:string = "amatter%{07u}.xyz".format(istep);
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer();
    // number of active particles + 4 edge-defining particles + boundary + solvent (optional)
    if (fluid_cpl) {
        myFileWriter.writeln(nworms * np + 4 + numPoints + numSol);
    } else {
        myFileWriter.writeln(nworms*np + 4 + numPoints);
    }
    myFileWriter.writeln("# 0");
    ic = 2;
    for iw in 1..nworms {
        dx = worms[iw,1].x - hxo2;
        dy = worms[iw,1].y - hyo2;
        xang = atan2(dy,dx);
        rx = -sin(xang);
        ry = cos(xang);
        dot = (worms[iw,1].x - worms[iw,np].x)*rx + (worms[iw,1].y - worms[iw,np].y)*ry;
        if (dot >= 0.0) {
            for i in 1..np {
                myFileWriter.writeln("A ",worms[iw,i].x," ",worms[iw,i].y," ", 0.0," ",worms[iw,i].vx," ",worms[iw,i].vy," ",0.0);
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",worms[iw,i].x," ",worms[iw,i].y," ", 0.0," ",worms[iw,i].vx," ",worms[iw,i].vy," ",0.0);
                //myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        }
    }
    for i in 1..numPoints {
        myFileWriter.writeln("I ",bound[i].x," ",bound[i].y," ",0.0," ",0.0," ",0.0," ",0.0);
        ic += 1;
    }
    if (fluid_cpl) {
        for i in 1..numSol {
            myFileWriter.writeln("S ",solvent[i].x," ",solvent[i].y," ",0.0," ",solvent[i].vx," ",solvent[i].vy," ",0.0);
            ic += 1;
        }
    }
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    //ic += 4;
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    writeln(filename,"\t",ic," lines written");
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}

proc write_macro(nsteps: int) {
    var filename:string = "energies.dat";
    try{
        var datfile = open(filename, ioMode.cw);
        var myFileWriter = datfile.writer();
        myFileWriter.writeln("step \t KEworm \t KEsol \n");
        for istep in 1..nsteps {
            myFileWriter.writeln(istep,"\t",KEworm_total[istep],"\t",KEsol_total[istep]);
        }
        datfile.fsync();
        writeln("energies.dat written");
    } catch e: Error {
        writeln(e);
    }
}

proc write_params() {
    var filename:string = "params.dat";
    try {
        var paramsfile = open(filename, ioMode.cw);
        var myFileWriter = paramsfile.writer();
        myFileWriter.writeln("np\t",np.type,"\t",np,"\n",
                             "nworms\t",nworms.type,"\t",nworms,"\n",
                             "nsteps\t",nsteps.type,"\t",nsteps,"\n",
                             "fdogic\t",fdogic.type,"\t",fdogic,"\n",
                             "walldrive\t",walldrive.type,"\t",walldrive,"\n",
                             "fdogicwall\t",fdogicwall.type,"\t",fdogicwall,"\n",
                             "fdep\t",fdep.type,"\t",fdep,"\n",
                             "fdepwall\t",fdepwall.type,"\t",fdepwall,"\n",
                             "diss\t",diss.type,"\t",diss,"\n",
                             "dt\t",dt.type,"\t",dt,"\n",
                             "kspring\t",kspring.type,"\t",kspring,"\n",
                             "kbend\t",kbend.type,"\t",kbend,"\n",
                             "length0\t",length0.type,"\t",length0,"\n",
                             "rcut\t",rcut.type,"\t",rcut,"\n",
                             "save_interval\t",save_interval.type,"\t",save_interval,"\n",
                             "boundary\t",boundary.type,"\t",boundary,"\n",
                             "fluid_cpl\t",fluid_cpl.type,"\t",fluid_cpl,"\n",
                             "debug\t",debug.type,"\t",debug,"\n",
                             "thermo\t",thermo.type,"\t",thermo,"\n",
                             "kbt\t",kbt.type,"\t",kbt,"\n",
                             "sigma\t",sigma.type,"\t",sigma,"\n",
                             "gamma\t",gamma,type,"\t",gamma,"\n",
                             "numPoints\t",numPoints.type,"\t",numPoints,"\n",
                             "numSol\t",numSol.type,"\t",numSol,"\n",
                             "fluid_offset\t",fluid_offset.type,"\t",fluid_offset,"\n"
                             )
        paramsfile.fsync();
        writeln("params.dat written");
    } catch e: Error {
        writeln(e);
    }
}

// HELPER FUNCTIONS
proc gaussRand(mean: real, stddev: real): real {
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

proc sudden_halt(istep:int) {
    write_xyzv(istep);
    write_macro(nsteps);
}


// OLD FUNCTIONS

proc worm_wall() {
    var dx:real, dy:real, r:real, r2:real, th:real,
           xwall:real, ywall:real, rr2:real, ffor:real,
           dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
    var ip1:int;
    //put worm-wall interactions here, and dissipation force proportional to velocity
    for iw in 1..nworms{
        for i in 1..np{
            //dissipation proportional to v relative to local average
            // TODO swap vxave with intvx
            worms[iw,i].fx = worms[iw,i].fx - diss*(worms[iw,i].vx - worms[iw,i].vxave);
            worms[iw,i].fy = worms[iw,i].fy - diss*(worms[iw,i].vy - worms[iw,i].vyave);
            //now that we have used them, zero out vxave and vyave, recalculate below
            worms[iw,i].vxave = worms[iw,i].vx;
            worms[iw,i].vyave = worms[iw,i].vy;
            nnab[iw, i] = 1;
            //calculate distance to the center
            dx = worms[iw, i].x - hxo2;
            dy = worms[iw, i].y - hyo2;
            r2 = (dx*dx + dy*dy);
            //if close enough to the wall, calculate wall forces
            //use the short cut-off
            if (r2 >= r2inside) {
                //find the nearest spot on the wall
                th = atan2(dy, dx);
                xwall = hxo2 + rwall*cos(th);
                ywall = hyo2 + rwall*sin(th);
                dx = xwall - worms[iw,i].x;
                dy = ywall - worms[iw,i].y;
                rr2 = dx*dx + dy*dy;
                r = sqrt(rr2);
                //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
                ffor = -48.0*rr2**(-7.0) + 24.0*rr2**(-4.0);
                worms[iw, i].fx += ffor*dx;
                worms[iw, i].fy += ffor*dy;

                //iwalldrive = 1; // made this a const parameter
                //Turning on dogic drive with the wall!!!
                if (iwalldrive == 1) {
                    //first calculate unit vector along the worm
                    ip1 = i + 1;
                    if (ip1 <= np) {
                        dxi = worms[iw, ip1].x - worms[iw, i].x;
                        dyi = worms[iw, ip1].y - worms[iw, i].y;
                    } else {
                        dxi = worms[iw, i].x - worms[iw, i - 1].x;
                        dyi = worms[iw, i].y - worms[iw, i - 1].y;
                    }
                    //make it a unit vector
                    ri = sqrt(dxi*dxi + dyi*dyi);
                    dxi = dxi/ri;
                    dyi = dyi/ri;
                    //calculate the unit vector along the wall
                    dxj = -sin(th); //for cardioid, make this a vector pointing to the next wall particle
                    dyj = cos(th);

                    //if the vectors are not antiparallel, reverse the vector along the wall
                    if (dxi*dxj + dyi*dyj > 0.0) {
                        dxj = -dxj;
                        dyj = -dyj;
                    }

                    //if the two vectors have any component pointing in opposite directions
                    if (dxi*dxj + dyi*dyj < 0.0) {
                        //Find the direction for the force...
                        dx = (dxi - dxj)/2.0;
                        dy = (dyi - dyj)/2.0;
                        //normalize the direction vector
                        ri = sqrt(dx*dx + dy*dy);
                        dx = dx/ri;
                        dy = dy/ri;

                        //turn on extra-strong driving force
                        ffx = fdogicwall*dx;
                        ffy = fdogicwall*dy;
                        worms[iw,i].fx += ffx;
                        worms[iw,i].fy += ffy;
                    }
                }
            }
        }
    }
}

proc worm_wall_new() {
    //put worm-wall interactions here, and dissipation force proportional to velocity
    forall iw in 1..nworms{
        var dx:real, dy:real, r:real, r2:real, th:real,
           xwall:real, ywall:real, rr2:real, ffor:real,
           dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
        var ip1,ib1:int;
        for i in 1..np{
            //dissipation proportional to v relative to local average
            // TODO swap vxave with intvx
            //worms[iw,i].fx = worms[iw,i].fx - diss*(worms[iw,i].vx - worms[iw,i].vxave);
            //worms[iw,i].fy = worms[iw,i].fy - diss*(worms[iw,i].vy - worms[iw,i].vyave);
            //now that we have used them, zero out vxave and vyave, recalculate below
            //worms[iw,i].vxave = worms[iw,i].vx;
            //worms[iw,i].vyave = worms[iw,i].vy;
            nnab[iw, i] = 1;

            // calculate the force on the boundaries.
            for ib in 1..numPoints  {
                //calculate distance to the wall
                dx = worms[iw,i].x-bound[ib].x;
                dy = worms[iw,i].y-bound[ib].y;
                r2 = (dx*dx + dy*dy);
                //if close enough to the wall, calculate wall forces
                //use the short cut-off
                if (r2 <= r2cut) {
                    ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
                    //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                    worms[iw,i].fx += ffor*dx;
                    worms[iw,i].fy += ffor*dy;
                    if (walldrive) {
                        //first calculate unit vector along the worm
                        ip1 = i + 1;
                        if (ip1 <= np) {
                            dxi = worms[iw, ip1].x - worms[iw, i].x;
                            dyi = worms[iw, ip1].y - worms[iw, i].y;
                        } else {
                            dxi = worms[iw, i].x - worms[iw, i - 1].x;
                            dyi = worms[iw, i].y - worms[iw, i - 1].y;
                        }
                        //make it a unit vector
                        ri = sqrt(dxi*dxi + dyi*dyi);
                        dxi = dxi/ri;
                        dyi = dyi/ri;
                        //calculate the unit vector along the wall
                        ib1 = ib + 1;
                        if (ib1 <= numPoints) {
                            dxj = bound[ib1].x - bound[ib].x;
                            dyj = bound[ib1].y - bound[ib].y;
                        } else {
                            dxj = bound[ib].x - bound[ib-1].x;
                            dyj = bound[ib].y - bound[ib-1].y;
                        }
                        ri = sqrt(dxj*dxj + dyj*dyj);
                        dxj = dxj/ri;
                        dyj = dyj/ri;

                        //if the vectors are not antiparallel, reverse the vector along the wall
                        if (dxi*dxj + dyi*dyj > 0.0) {
                            dxj = -dxj;
                            dyj = -dyj;
                        }
                        //if the two vectors have any component pointing in opposite directions
                        if (dxi*dxj + dyi*dyj < 0.0) {
                            //Find the direction for the force...
                            dx = (dxi - dxj)/2.0;
                            dy = (dyi - dyj)/2.0;
                            //normalize the direction vector
                            ri = sqrt(dx*dx + dy*dy);
                            dx = dx/ri;
                            dy = dy/ri;

                            //turn on extra-strong driving force
                            ffx = fdogicwall*dx;
                            ffy = fdogicwall*dy;
                            worms[iw,i].fx += ffx;
                            worms[iw,i].fy += ffy;
                        }
                }
                }
            }
        }
    }
}

proc cell_sort_old(itime:int) {
    var dddx:real,dddy:real,r2:real,riijj:real,ffor:real,ffx:real,ffy:real,dxi:real,dxj:real,
        ri:real,rj:real,r:real,dx:real,dy:real,dyi:real,dyj:real;
    var iworm:int,jworm:int,ip:int,jp:int,ii:int,jj:int,kk:int,i:int,ip1:int,
        jp1:int,scell:int,scnab:int,inogo:int,icnab:int,jcnab:int,icell:int,jcell:int;

    for iworm in 1..nworms{
        for ip in 1..np {
            ii = (iworm-1)*np+ip; //unique particle id 1<=ii<=nworms*np
            icell = 1+floor(worms[iworm,ip].x/dcell):int; //finds where particle is in grid
            jcell = 1+floor(worms[iworm,ip].y/dcell):int;
            if ((icell > nxcell) || (icell < 1)) {
                writeln("nxcell=",nxcell," icell=",icell);
                writeln(worms[iworm,ip].x/dcell);
                writeln(floor(worms[iworm,ip].x/dcell):int);
                writeln(1+floor(worms[iworm,ip].x/dcell):int);
                writeln("icell out of bounds\t",iworm," ",ip," ",worms[iworm,ip].x," ",worms[iworm,ip].vx," ",worms[iworm,ip].fx);
                write_xyz(itime);
                halt();
            }
            if ((jcell > nycell) || (jcell < 1)) {
                writeln("nxcell=",nycell," icell=",jcell);
                writeln(worms[iworm,ip].y/dcell);
                writeln(floor(worms[iworm,ip].y/dcell));
                writeln("jcell out of bounds\t",iworm," ",ip," ",worms[iworm,ip].y," ",worms[iworm,ip].vx," ",worms[iworm,ip].fx);
                write_xyz(itime);
                halt();
            }
            scell = icell + (jcell-1)*nxcell; // 1d-ndexing for 2d cells
            if ((scell > ncells) || (scell < 1)) {
                writeln("scell out of bounds",scell,"\t",icell,"\t",jcell);
            }
            ipointto[ii] = hhead[scell]; //
            hhead[scell] = ii;
        }
    }
    for icell in 1..nxcell {
        for jcell in 1..nycell {
            scell = icell + (jcell-1)*nxcell;
            if (hhead[scell] != -1){
                //there are particles in the cell called scell so
                //lets check all the neighbor cells
                for idir in 1..9 {
                    icnab = icell + ddx[idir];
                    if (icnab > nxcell) {break;}
                    if (icnab == 0) {break;}

                    jcnab = jcell + ddy[idir];
                    if (jcnab > nycell) {break;}
                    if (jcnab == 0) {break;}
                    scnab = icnab + (jcnab-1)*nxcell; //1d neighbor
                    if (hhead[scnab] != -1) {
                        //there are particles in the cell called scnab
                        ii = hhead[scell]; // ii is the # of the head particle
                        while (ii > 0) {

                            iworm = 1 + ((ii - 1)/np):int; // find which worm ii is in
                            ip = ii - np*(iworm - 1); // which particle in the worm is ii?
                            jj = hhead[scnab];//  head particle of neighboring cell

                            while (jj > 0) {
                                jworm = 1 + ((jj - 1)/np):int;
                                jp = jj - np*(jworm - 1);
                                inogo = 0;
                                if ((iworm == jworm) && (abs(ip-jp) <= 2)) {
                                    // on the same worm and close means no interaction calculated here
                                    inogo = 1;
                                }

                                if ((ii < jj) && (inogo == 0)) {
                                    dddx = worms[jworm, jp].x - worms[iworm, ip].x;
                                    dddy = worms[jworm, jp].y - worms[iworm, ip].y;
                                    r2 = dddx**2 + dddy**2;
                                    riijj = sqrt(r2);
                                    //add attractive force fdep between all pairs
                                    if (r2 <= r2cutsmall) {
                                        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
                                        ffx = ffor*dddx;
                                        ffy = ffor*dddy;
                                        worms[iworm,ip].fx += ffx;
                                        worms[jworm,jp].fx -= ffx;
                                        worms[iworm,ip].fy += ffy;
                                        worms[jworm,jp].fy -= ffy;

                                        //take these neighbors into account in calculating vxave and vyave
                                        //worms[iworm,ip].vxave += worms[jworm,jp].vx;
                                        //worms[iworm,ip].vyave += worms[jworm,jp].vy;
                                        nnab[iworm, ip] = nnab[iworm, ip] + 1;
                                        //worms[jworm,jp].vxave += worms[iworm,ip].vx;
                                        //worms[jworm,jp].vyave += worms[iworm,ip].vy;
                                        nnab[jworm, jp] = nnab[jworm, jp] + 1;

                                        //add 'dogic drive' to interacting pairs
                                        //first calculate unit vectors along each worm
                                        ip1 = ip + 1;
                                        if (ip1 <= np) {
                                            dxi = worms[iworm,ip1].x - worms[iworm, ip].x;
                                            dyi = worms[iworm,ip1].y - worms[iworm, ip].y;
                                        } else {
                                            dxi = worms[iworm, ip].x - worms[iworm, ip - 1].x;
                                            dyi = worms[iworm, ip].y - worms[iworm, ip - 1].y;
                                        }

                                        jp1 = jp + 1;
                                        if (jp1 <= np) {
                                            dxj = worms[jworm, jp1].x - worms[jworm, jp].x;
                                            dyj = worms[jworm, jp1].y - worms[jworm, jp].x;
                                        } else {
                                            dxj = worms[jworm, jp].x - worms[jworm, jp - 1].x;
                                            dyj = worms[jworm, jp].y - worms[jworm, jp - 1].y;
                                        }


                                        //if the two vectors have any component pointing in opposite directions
                                        if (dxi*dxj + dyi*dyj <= 0.0) {
                                            //normalize those vectors to make them unit vectors
                                            ri = sqrt(dxi*dxi + dyi*dyi);
                                            dxi = dxi/ri;
                                            dyi = dyi/ri;

                                            rj = sqrt(dxj*dxj + dyj*dyj);
                                            dxj = dxj/rj;
                                            dyj = dyj/rj;
                                            //now they are both unit vectors. Find the direction for the force...

                                            dx = (dxi - dxj)/2.0;
                                            dy = (dyi - dyj)/2.0;

                                            //normalize

                                            r = sqrt(dx*dx + dy*dy);
                                            dx = dx/r;
                                            dy = dy/r;

                                            //add an extra attractive component where kinesin drive is present

                                            ffx = fdogic*(dx) + 0.7*dddx/riijj;
                                            ffy = fdogic*(dy) + 0.7*dddy/riijj;

                                            //ffx=fdogic*(dx)
                                            //ffy=fdogic*(dy)

                                            worms[iworm,ip].fx += ffx;
                                            worms[jworm,jp].fx -= ffx;
                                            worms[iworm,ip].fy += ffy;
                                            worms[jworm,jp].fy -= ffy;
                                        }
                                    }
                                }
                                jj = ipointto[jj];
                                }
                            ii = ipointto[ii];
                        }
                    }
                }
            }
        }
    }
}

proc fluid_force() {
    forall binid in binSpace {
        if (bins[binid].ncount > 1) {
            // calculate the forces between atoms inside each bin
            for icount in 0..bins[binid].ncount-2 {
                var i = bins[binid].atoms[icount];
                for jcount in (icount+1)..bins[binid].ncount-1 {
                    var j = bins[binid].atoms[jcount];
                    lj_thermo(i,j,r2cut);
                }
            }
        }
    }
    // odd neighbors to the east (1)
    forall binid in binSpaceiodd {
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }
                }
            }
        }
    }
    // even neighbors to the east (1)
    forall binid in binSpaceieven {
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != 1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }
                }
            }
        }
    }

    // odd neighbors to the NE (2) (i+1,j+1)
    forall binid in binSpaceiodd {
        var binidnbor = bins[binid].neighbors[2];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }
                }
            }
        }
    }

    // even neighbors to the NE (2)
    forall binid in binSpaceieven {
        var binidnbor = bins[binid].neighbors[2];
        if (binidnbor != -1) { // check if neighbor is valid
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) { // check if neighbor has any atoms at all
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }
                }
            }
        }
    }

    // odd neighbors to the N (3)
    forall binid in binSpacejodd {
        var binidnbor = bins[binid].neighbors[3];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }

                }
            }
        }
    }

    // even neighbors to the N (3)
    forall binid in binSpacejeven {
        var binidnbor = bins[binid].neighbors[3];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }

                }
            }
        }
    }

    // odd neighbors to the NW (4)
    forall binid in binSpaceiodd {
        var binidnbor = bins[binid].neighbors[4];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }

                }
            }
        }
    }

    // even neighbors to the NW (4)
    forall binid in binSpaceieven {
        var binidnbor = bins[binid].neighbors[4];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        lj_thermo(i,j,r2cut);
                    }

                }
            }
        }
    }

        //fluid-boundary
    // var r2cutsol = r2cut
    var r2cutsol = sigma*2.0**(1.0/6.0);
    forall i in 1..numSol {
        // calculate the force on the boundaries.
        for ib in 1..numPoints  {
            lj(ib,i,r2cutsol);

        }
    }

    // fluid-worms
    var dz = fluid_offset;
    for i in 1..numSol{
        var dx,dy,r2,ffor,ffx,ffy:real;
        for iw in 1..nworms {
            for ip in 1..np {
                dx = solvent[i].x - worms[iw,ip].x;
                dy = solvent[i].y - worms[iw,ip].y;
                r2 = (dx*dx + dy*dy + dz*dz);
                if (r2 <= r2cutsmall) {
                    ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                    ffx = ffor*dx;
                    ffy = ffor*dy;
                    solvent[i].fx += ffx;
                    solvent[i].fy += ffy;
                    worms[iw,ip].fx -= ffx;
                    worms[iw,ip].fy -= ffy;
                }
            }
        }
    }
}

proc fluid_multistep() {
    var iterations = 1;
    var timestep = dt/iterations;
    for i in 1..iterations {
        fluid_step(timestep);
    }
}