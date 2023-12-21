use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
use CTypes; // for chpl string -> c string pointer conversion
use DynamicIters;
// user-defined modules
import Structs;
import C; // import C-extension module for logging/appending


const numTasks = here.numPUs();
// configuration
config const np = 100,//16,
            nworms = 300,//625,
            nsteps = 12000000    ,//00,
            fdogic = 0.0,
            walldrive = true,
            fdogicwall = 0.001,
            fdep = 1.0,// TODO: change to 4.0?
            fdepwall = 0.0,
            diss = 0.02,
            dt = 0.005, //0.02
            kspring = 57.146436,
            kbend = 40.0,
            length0 = 0.8, //particle spacing on worms
            rcut = 2.5,
            save_interval = 4000,
            boundary = 1, // 1 = circle, 2 = cardioid, 3 = channel, 4 = torus
            fluid_cpl = true,
            debug = false,
            thermo = true, // turn thermostat on?
            kbt = 0.00001, //0.25
            //numSol = 7000, // cardiod number of solution particles
            fluid_rho = 0.05,//8000, // disk number of solution particles
            sigma = 2.0,
            worm_particle_mass = 4.0,
            L = 2.0; // thickness of cell

const io_interval = 500;
// variables
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1.0/6.0),
      r2cutsmall = rcutsmall*rcutsmall,
      rwall = 199,//125.0*rcutsmall*sqrt(2.0),
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
      nstepso1e5 = nsteps/10000,
      gnoise = 0.80/sqrt(10.0)*0.8,
      dt2o2 = dt*dt*0.50,
      dto2 = dt*0.50,
      length2 = 2.0*length0,
      lengthmax = (length0)*((np - 1):real),
      r2inside = (rwall - rcutsmall) * (rwall-rcutsmall),
      a = 0.24, // layer spacing of worms in init_worms?
      gamma = 3.0, // frictional constant for dissipative force (~1/damp)
      dpd_ratio = 1.0,
      numPoints = 2000,//1200//589, //number of boundary points (for circle w/ r-75)
      fluid_offset = rcutsmall*sigma;//3.0; // z-offset of fluid


var wormsDomain: domain(2) = {1..nworms,1..np};
var worms: [wormsDomain] Structs.Particle;
Structs.ptc_init_counter = 1;

var numSol:int;
numSol = init_fluid_count();
var solvent: [1..numSol] Structs.Particle;
Structs.ptc_init_counter = 1;

var bound: [1..numPoints] Structs.Particle;
var savex: [1..np] real(64);
var savey: [1..np] real(64);
var savez: [1..np] real(64);
var ireverse: [1..nworms] int;
var ddx: [1..9] int;
var ddy: [1..9] int;
var hhead: [1..ncells] int; //
var ipointto: [1..nworms*np+numPoints] int; // linked list, every particle points to another particle
var nnab: [wormsDomain] int = 1;
var KEworm: [1..nworms*np] real; // Kinetic energy of individual worm particles
var KEworm_local: [1..nworms*np] real; // this is the for velocity minus the average velocity (tries to measure thermal fluctuations)
var AMworm: [1..nworms*np] real; // angular momentum of individual worm particles
var KEsol: [1..numSol] real; // Kinetic energy of individual fluid particles
var AMsol: [1..numSol] real; // angular momentum of individual fluid particles
var KEworm_total: [1..nsteps] real;
var KEworm_local_total: [1..nsteps] real;
var AMworm_total: [1..nsteps] real;
var KEsol_total: [1..nsteps] real;
var AMsol_total: [1..nsteps] real;

var numBins = ceil(hx/rcut):int;
writeln("numBins:\t",numBins);
//const binSpace = {1..numBins, 1..numBins};
const binSpace = {1..numBins*numBins};
var bins : [1..numBins*numBins] Structs.Bin;
var binSpaceiodd : [1..(numBins*numBins)/2] int;
var binSpacejodd : [1..(numBins*numBins)/2] int;
var binSpaceieven : [1..(numBins*numBins)/2] int;
var binSpacejeven : [1..(numBins*numBins)/2] int;

var randStream = new RandomStream(real); // creating random number generator

var logfile:string = "amatter.log";
var t = 0.0;
var total_time = 0.0;
var ct: stopwatch, wt:stopwatch, xt:stopwatch; //calc time, io time, totaltime
//main
proc main() {
    init_log(logfile);
    write_log(logfile, "=============== AMATTER3D ===============");
    write_log(logfile,date.today():string+":Matt Deutsch, Kent State University");
    write_log(logfile,"starting... "+numTasks:string+" procs");
    
    // save params to file
    write_params();
    // initialize the alternating loops over bins
    init_binspace();
    // initialize bins and neighboring bin lists
    init_bins();
    // initialize worms
    init_worms();

    //setting mass of worms to be differents
    for iw in 1..nworms {
        for i in 1..np {
            worms[iw,i].m = worm_particle_mass;
        }
    }

    if (fluid_cpl) {
        solvent = init_fluid(solvent, numSol);
    }
    // populate the bins with lists of atoms
    update_cells(0);
    //equilibrate the fluid
   //  if (fluid_cpl) {
   //    for istep in 1..5000 {
   //        var ioper = 5000/10;
   //        if (istep%ioper == 0) {
   //            write_log(logfile,"fluid equilibration..."+istep:string);
   //        }
   //        fluid_step(0,dt);
   //    }
   //  }
    update_cells(0); //again after fluid
    write_log(logfile,"numSol\t"+numSol:string);
    write_log(logfile,"fluid equilibrated...5000dt");

    var macro_filename:string = "energies.dat";
    init_macro(macro_filename);
    write_xyzv(0);
    
    //setting up stopwatch
    var log_str:string;
    xt.start();
    for itime in 1..nsteps {
        //writeln(itime);
        t = (itime:real) *dt;

        if (itime % io_interval == 0) {
            xt.stop();
            total_time += xt.elapsed();
            var out_str:string = "Step: "+itime:string+"\t"+
                      (io_interval/xt.elapsed()):string+"iter/s\tCalc:"+
                      ((ct.elapsed()/total_time)*100):string+"%\tIO:"+
                      ((wt.elapsed()/total_time)*100):string+" %\tElapsed:"+
                      total_time:string+" s\t Est:"+
                      ((nsteps-itime)*(total_time/itime)):string+" s";
            write_log(logfile,out_str);
            //writeln((+ reduce solvent.vx),"\t",(+ reduce solvent.vy));
            //writeln((+ reduce solvent.x),"\t",(+ reduce solvent.y));
            xt.restart();
            }
	    ct.start();
        // first update positions and store old forces
        update_pos(itime);

        update_cells(itime);

        intraworm_forces();

        calc_forces(itime,dt);
        if (fluid_cpl) {
            KEsol_total[itime] = (+ reduce KEsol);
            AMsol_total[itime] = (+ reduce AMsol);
        } else {
            KEsol_total[itime] = 0.0;
        }

        update_vel();
        KEworm_total[itime] = (+ reduce KEworm);
        KEworm_local_total[itime] = (+ reduce KEworm_local);
        AMworm_total[itime] = (+ reduce AMworm);

	    ct.stop();
        if (itime % nstepso1e5 == 0){
            write_macro(macro_filename,itime);
        }
        if (itime % save_interval == 0){
            wt.start();
            write_xyzv(itime);
            wt.stop();
        }
    }
    //finalize();
    xt.stop();
    //write_macro(nsteps);
    write_log(logfile,"Total Time:"+total_time:string+" s");
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
    var thetanow = 85.0*pi :real; // changes the initial radius of annulus
    var rmin = a*thetanow;
    write_log(logfile,"nworms\t"+nworms:string);
    write_log(logfile,"np\t"+np:string);
    write_log(logfile,"rwall\t"+rwall:string);
    write_log(logfile,"nxcells\t"+nxcell:string+"\t nycells\t"+nycell:string);
    write_log(logfile,"density\t"+density:string);
    write_log(logfile,"dt\t"+dt:string);
    write_log(logfile,"opt numPoints\t"+(ceil(2*pi*rwall)/rcutsmall):string);
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
                worms[iw,i].z = randStream.getNext()*(L);
                xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                //TODO give them an initial velocity going around the circle
                worms[iw,i].ptype = 1;
                worms[iw,i].m = 1.0; // setting mass
                // vx[iw,i] = 0.0; // NOTE: this is all set by initializer
                // vy[iw,i] = 0.0;
                worms[iw,i].vxave = 0.0;
                worms[iw,i].vyave = 0.0;
                worms[iw,i].vzave = 0.0;
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
        write_log(logfile,"cardioid cirum: "+totalArcLength:string+"\t"+(totalArcLength/r2cutsmall):string);
        for i in 1..numPoints {
            equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
            thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
        }

        for i in 1..numPoints {
            bound[i].x = ca * (1 - cos(thetaValues[i])) * cos(thetaValues[i]) + hxo2 + ca;
            bound[i].y = ca * (1 - cos(thetaValues[i])) * sin(thetaValues[i]) + hyo2;
            bound[i].ptype = 3;
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
                worms[iw,i].z = 0.0;
                xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                //TODO give them an initial velocity going around the circle
                worms[iw,i].ptype = 1;
                worms[iw,i].m = 1.0; // setting mass
                // vx[iw,i] = 0.0;
                // vy[iw,i] = 0.0;
                worms[iw,i].vxave = 0.0;
                worms[iw,i].vyave = 0.0;
                worms[iw,i].vzave = 0.0;
                // fx[iw,i] = 0.0;
                // fy[iw,i] = 0.0;
                // fxold[iw,i] = 0.0;
                // fyold[iw,i] = 0.0;
            }
            thetanow += 2.0*dth;
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
                savez[i] = worms[iw,i].z;
            }
            for ip in 1..np {
                worms[iw,ip].x = savex[np+1-ip];
                worms[iw,ip].y = savey[np+1-ip];
                worms[iw,ip].z = savez[np+1-ip];
            }
        }
    }

    //init macro variable arrays
    for i in 1..nworms*np{
        KEworm[i] = 0.0;
    }
    write_log(logfile,"init done");
}

proc update_pos(itime:int) {
    //forall iw in 1..nworms {
    forall iw in 1..nworms {
       forall i in 1..np {
            worms[iw,i].fx = worms[iw,i].fx/worms[iw,i].m;
            worms[iw,i].fy = worms[iw,i].fy/worms[iw,i].m;
            worms[iw,i].fz = worms[iw,i].fz/worms[iw,i].m;

            //dissipation proportional to v relative to local average
            //worms[iw,i].fx = worms[iw,i].fx - diss*(worms[iw,i].vx - worms[iw,i].vxave);
            //worms[iw,i].fy = worms[iw,i].fy - diss*(worms[iw,i].vy - worms[iw,i].vyave);

            worms[iw,i].x = worms[iw,i].x + worms[iw,i].vx*dt + worms[iw, i].fx*dt2o2;
            worms[iw,i].y = worms[iw,i].y + worms[iw,i].vy*dt + worms[iw, i].fy*dt2o2;
            worms[iw,i].z = worms[iw,i].z + worms[iw,i].vz*dt + worms[iw, i].fz*dt2o2;

            worms[iw,i].fxold = worms[iw,i].fx;
            worms[iw,i].fyold = worms[iw,i].fy;
            worms[iw,i].fzold = worms[iw,i].fz;

            worms[iw,i].fx = 0.0;
            worms[iw,i].fy = 0.0;
            worms[iw,i].fz = 0.0;

            // periodic boundary conditions top and bottom for worms
            //periodic boundary conditions
            if worms[iw,i].z > L {
                //worms[iw,i].z = worms[iw,i].z - L;
                worms[iw,i].z = L;
            } else if worms[iw,i].z < 0.0 {
                //worms[iw,i].z = worms[iw,i].z + L;
                worms[iw,i].z = 0.0;
            }
        }
    }
    fluid_pos(dt);
}

proc intraworm_forces() {
    //first set of springs nearest neighbor springs
    forall iw in 1..nworms {
        var ip1:int,r:real,ff:real,ffx:real,ffy:real,ffz:real,dx:real,dy:real,dz:real;
        for i in 1..np-1 {
            ip1 = i + 1;
            dx = worms[iw,ip1].x - worms[iw,i].x;
            dy = worms[iw,ip1].y - worms[iw,i].y;
            dz = worms[iw,ip1].z - worms[iw,i].z;
            r = sqrt(dx*dx + dy*dy + dz*dz);

            ff = -kspring*(r - length0)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            ffz = ff*dz;
            worms[iw,ip1].fx += ffx;
            worms[iw,i].fx -= ffx;
            worms[iw,ip1].fy += ffy;
            worms[iw,i].fy -= ffy;
            worms[iw,ip1].fz += ffz;
            worms[iw,i].fz -= ffz;
        }
    }
    //bond bending terms
    forall iw in 1..nworms {
        var i3:int,i4:int,x2:real,x3:real,x4:real,y2:real,y3:real,y4:real,z2:real,z3:real,z4:real,z23:real,z34:real,y23:real,y34:real,x23:real,x34:real,r23:real,r34:real,cosvalue:real;
	    var sinvalue:real,ff:real,dot:real,fac:real,f2x:real,f2y:real,f2z:real,f3x:real,f3y:real,f3z:real,f4x:real,f4y:real,f4z:real;
        for i2 in 1..(np-2) {
            i3 = i2 + 1;
            i4 = i2 + 2;
            //print*, i2,i3,i4
            x2 = worms[iw, i2].x;
            y2 = worms[iw, i2].y;
            z2 = worms[iw, i2].z;

            x3 = worms[iw, i3].x;
            y3 = worms[iw, i3].y;
            z3 = worms[iw, i3].z;

            x4 = worms[iw, i4].x;
            y4 = worms[iw, i4].y;
            z4 = worms[iw, i4].z;

            z23 = z3-z2;
            z34 = z4-z2;

            y23 = y3 - y2;
            y34 = y4 - y3;

            x23 = x3 - x2;
            x34 = x4 - x3;

            r23 = sqrt(x23*x23 + y23*y23 + z23*z23);
            r34 = sqrt(x34*x34 + y34*y34 + z34*z34);

            cosvalue = abs(x23*x34 + y23*y34 + z23*z34)/(r23*r34);
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

            dot = x23*x34 + y23*y34 + z23*z34;
            fac = dot/(r23*r23);

            f2x = ff*(x34 - fac*x23);
            f2y = ff*(y34 - fac*y23);
            f2z = ff*(z34 - fac*z23);

            fac = dot/(r34*r34);
            f4x = ff*(fac*x34 - x23);
            f4y = ff*(fac*y34 - y23);
            f4z = ff*(fac*z34 - z23);

            f3x = -f2x - f4x;
            f3y = -f2y - f4y;
            f3z = -f2z - f4z;

            worms[iw, i2].fx += f2x;
            worms[iw, i2].fy += f2y;
            worms[iw, i2].fz += f2z;

            worms[iw, i3].fx += f3x;
            worms[iw, i3].fy += f3y;
            worms[iw, i3].fz += f3z;

            worms[iw, i4].fx += f4x;
            worms[iw, i4].fy += f4y;
            worms[iw, i4].fz += f3z;
        }
    }
}

proc calc_forces(istep:int,dt:real) {
    var binDomain: domain(1) = {1..(numBins*numBins)/2};
    forall binid in adaptive(binSpace) {
        if (bins[binid].ncount > 1) {
            // calculate the forces between atoms inside each bin
            for icount in 0..bins[binid].ncount-2 {
                var i = bins[binid].atoms[icount];
                var itype = bins[binid].types[icount];
                for jcount in (icount+1)..bins[binid].ncount-1 {
                    var j = bins[binid].atoms[jcount];
                    var jtype = bins[binid].types[jcount];
                    cell_forces(i,j,itype,jtype);
                }
            }
        }
    }
    // odd neighbors to the east (1)
    forall i in adaptive(binDomain) {
        var binid = binSpaceiodd[i];
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
    forall i in adaptive(binDomain) {
        var binid = binSpaceieven[i];
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

    // odd neighbors to the NE (2) (i+1,j+1)
    forall i in adaptive(binDomain) {
        var binid = binSpaceiodd[i];
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
    forall i in adaptive(binDomain) {
        var binid = binSpaceieven[i];
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
    forall i in adaptive(binDomain) {
        var binid = binSpacejodd[i];
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
    forall i in adaptive(binDomain) {
        var binid = binSpacejeven[i];
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
    forall i in adaptive(binDomain) {
        var binid = binSpaceiodd[i];
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
    forall i in adaptive(binDomain) {
        var binid = binSpaceieven[i];
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

inline proc cell_forces(i:int,j:int,itype:int,jtype:int) {
    //worm-worm interaction 3D
    if (itype == 1) && (jtype == 1) {
        worm_cell_forces(i,j);
    }
    //worm-boundary interaction 2D
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
    //solvent-worm interaction 2D
    if ((itype == 2) && (jtype == 1)) {
        var dx,dy,dz,r2,ffor,ffx,ffy:real;
        var iw,ip:int;
        iw = 1 + ((j - 1)/np):int; // find which worm j is in
        ip = j - np*(iw - 1); // which particle in the worm is j?
        dz = fluid_offset;
        dx = solvent[i].x - worms[iw,ip].x;
        dy = solvent[i].y - worms[iw,ip].y;
        r2 = (dx*dx + dy*dy + dz*dz);
        if (r2 <= r2cut) {
            //writeln("solvent worm ",sqrt(r2),"\t",sqrt(r2cut));
            //writeln("solvent worm ",r2,"\t",r2cut);
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
        ip = i - np*(iw - 1); // which particle in the worm is j?
        dz = fluid_offset;
        dx = solvent[j].x - worms[iw,ip].x;
        dy = solvent[j].y - worms[iw,ip].y;
        r2 = (dx*dx + dy*dy + dz*dz);
        if (r2 <= r2cut) {
            //writeln("worm solvent ",sqrt(r2),"\t",sqrt(r2cut));
            //writeln("worm solvent",r2,"\t",r2cut);
            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
            ffx = ffor*dx;
            ffy = ffor*dy;
            solvent[j].fx += ffx;
            solvent[j].fy += ffy;
            worms[iw,ip].fx -= ffx;
            worms[iw,ip].fy -= ffy;
        }
    }
    //solvent-solvent interaction 2D
    if (itype == 2) && (jtype == 2) {
        lj_thermo(i,j,r2cutsmall);
    }
    //solvent-boundary interaction 2D
    if ((itype == 2) && (jtype == 3)) {
        lj(j,i,sigma*r2cutsmall);
    } else if ((itype == 3) && (jtype == 2)) {
        lj(i,j,sigma*r2cutsmall);
    }
}

inline proc worm_cell_forces(i:int,j:int) {
   var dddx:real,dddy:real,dddz:real,r2:real,riijj:real,ffor:real,ffx:real,ffy:real,ffz:real,dxi:real,dxj:real,
       ri:real,rj:real,r:real,dx:real,dy:real,dz:real,dyi:real,dyj:real,dzi:real,dzj:real,r2shift:real;
   var dvx:real,dvy:real,dvz:real,rhatx:real,rhaty:real,rhatz:real,omega:real,fdissx:real,fdissy:real,fdissz:real,gauss:real,frand:real;
   var iworm:int,jworm:int,ip:int,jp:int,ip1:int,jp1:int,inogo:int;
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
      dddz = worms[jworm, jp].z - worms[iworm, ip].z;
      r2 = dddx**2 + dddy**2 + dddz**2;
      riijj = sqrt(r2);
      //add attractive force fdep between all pairs
      if (r2 <= r2cutsmall) {
            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
            //ffor = -48.0*(r2-r2shift)**(-7.0) + 24.0*(r2-r2shift)**(-4.0) + fdep/riijj;
            ffx = ffor*dddx;
            ffy = ffor*dddy;
            ffz = ffor*dddz;
            worms[iworm,ip].fx += ffx;
            worms[jworm,jp].fx -= ffx;
            worms[iworm,ip].fy += ffy;
            worms[jworm,jp].fy -= ffy;
            worms[iworm,ip].fz += ffz;
            worms[jworm,jp].fz -= ffz;
            

            // DPD thermostat
            //adding dissipative force
            dvx = worms[jworm,jp].vx - worms[iworm,ip].vx;
            dvy = worms[jworm,jp].vy - worms[iworm,ip].vy;
            dvz = worms[jworm,jp].vz - worms[iworm,ip].vz;

            r = sqrt(r2);
            rhatx = dddx/r;
            rhaty = dddy/r;
            rhatz = dddz/r;

            omega = (1.0-r2/r2cutsmall);

            fdissx = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty + dvz*rhatz)*rhatx; //gamma = 1/damp (proportional to friction force)
            fdissy = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty + dvz*rhatz)*rhaty;
            //fdissz = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty + dvz*rhatz)*rhatz;

            worms[iworm,ip].fx -= fdissx;
            worms[iworm,ip].fy -= fdissy;
            //worms[iworm,ip].fz -= fdissz;

            worms[jworm,jp].fx += fdissx;
            worms[jworm,jp].fy += fdissy;
            //worms[jworm,jp].fz += fdissz;

            // adding random forces
            gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
            gauss = randStream.getNext();
            frand = dpd_ratio*(1.0/sqrt(dt))*sqrt(omega)*gauss*sqrt(2.0*kbt*gamma);

            worms[iworm,ip].fx += frand*rhatx;
            worms[iworm,ip].fy += frand*rhaty;
            //worms[iworm,ip].fz += frand*rhatz;

            worms[jworm,jp].fx -= frand*rhatx;
            worms[jworm,jp].fy -= frand*rhaty;
            //worms[jworm,jp].fz -= frand*rhatz;

            //vxave,vyave
            worms[iworm,ip].vxave += worms[jworm,jp].vx;
            worms[iworm,ip].vyave += worms[jworm,jp].vy;
            worms[iworm,ip].vzave += worms[jworm,jp].vz;

            worms[jworm,jp].vxave += worms[iworm,ip].vx;
            worms[jworm,jp].vyave += worms[iworm,ip].vy;
            worms[jworm,jp].vzave += worms[iworm,ip].vz;

            nnab[iworm,ip] += 1;
            nnab[jworm,jp] += 1;

            //add 'dogic drive' to interacting pairs
            //first calculate unit vectors along each worm
            ip1 = ip + 1;
            if (ip1 <= np) {
               dxi = worms[iworm,ip1].x - worms[iworm, ip].x;
               dyi = worms[iworm,ip1].y - worms[iworm, ip].y;
               dzi = worms[iworm,ip1].z - worms[iworm, ip].z;
            } else {
               dxi = worms[iworm, ip].x - worms[iworm, ip - 1].x;
               dyi = worms[iworm, ip].y - worms[iworm, ip - 1].y;
               dzi = worms[iworm, ip].z - worms[iworm, ip - 1].z;
            }

            jp1 = jp + 1;
            if (jp1 <= np) {
               dxj = worms[jworm, jp1].x - worms[jworm, jp].x;
               dyj = worms[jworm, jp1].y - worms[jworm, jp].y;
               dzj = worms[jworm, jp1].z - worms[jworm, jp].z;
            } else {
               dxj = worms[jworm, jp].x - worms[jworm, jp - 1].x;
               dyj = worms[jworm, jp].y - worms[jworm, jp - 1].y;
               dzj = worms[jworm, jp].z - worms[jworm, jp - 1].z;
            }


            //if the two vectors have any component pointing in opposite directions
            if (dxi*dxj + dyi*dyj + dzi*dzj<= 0.0) {
               //normalize those vectors to make them unit vectors
               ri = sqrt(dxi*dxi + dyi*dyi *dzi*dzi);
               dxi = dxi/ri;
               dyi = dyi/ri;
               dzi = dzi/ri;

               rj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);
               dxj = dxj/rj;
               dyj = dyj/rj;
               dzj = dzj/rj;
               //now they are both unit vectors. Find the direction for the force...

               dx = (dxi - dxj)/2.0;
               dy = (dyi - dyj)/2.0;
               dz = (dzi - dzj)/2.0;

               //normalize

               r = sqrt(dx*dx + dy*dy + dz*dz);
               dx = dx/r;
               dy = dy/r;
               dz = dz/r;

               //add an extra attractive component where kinesin drive is present

               ffx = fdogic*(dx) + 0.7*dddx/riijj;
               ffy = fdogic*(dy) + 0.7*dddy/riijj;
               ffz = fdogic*(dz) + 0.7*dddz/riijj;

               worms[iworm,ip].fx += ffx;
               worms[jworm,jp].fx -= ffx;
               worms[iworm,ip].fy += ffy;
               worms[jworm,jp].fy -= ffy;
               worms[iworm,ip].fz += ffz;
               worms[jworm,jp].fz -= ffz;
            }
      }
   }
}

inline proc dogic_wall(iw:int,ip:int,ib:int){
    var dx:real, dy:real, r:real, r2:real, th:real,
        xwall:real, ywall:real, rr2:real, ffor:real,
        dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
    var ip1,ib1:int;
    //calculate distance to the wall
    dx = worms[iw,ip].x-bound[ib].x;
    dy = worms[iw,ip].y-bound[ib].y;
    r2 = (dx*dx + dy*dy);
    //if close enough to the wall, calculate wall forces
    //use the short cut-off
    if (r2 <= r2cut) {
        r = sqrt(r2);
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdepwall/r;
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        ffor = (1/r)**4.0; //TODO raise this to a higher power to get the worms closer to the wall? try ^6 or ^8
        worms[iw,ip].fx += ffor*dx;
        worms[iw,ip].fy += ffor*dy;
        if (walldrive) {
            //first calculate unit vector along the worm
            ip1 = ip + 1;
            if (ip1 <= np) {
                dxi = worms[iw, ip1].x - worms[iw, ip].x;
                dyi = worms[iw, ip1].y - worms[iw, ip].y;
            } else {
                dxi = worms[iw, ip].x - worms[iw, ip - 1].x;
                dyi = worms[iw, ip].y - worms[iw, ip - 1].y;
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
                worms[iw,ip].fx += ffx;
                worms[iw,ip].fy += ffy;
            }
        }
    }
}

proc update_vel() {
    forall iw in 1..nworms {
        forall i in 1..np {
            worms[iw, i].vx += dto2*(worms[iw,i].fx + worms[iw,i].fxold);
            worms[iw, i].vy += dto2*(worms[iw,i].fy + worms[iw,i].fyold);
            //worms[iw, i].vz += dto2*(worms[iw,i].fz + worms[iw,i].fzold);

            worms[iw, i].vxave = worms[iw, i].vxave/nnab[iw, i]:real;
            worms[iw, i].vyave = worms[iw, i].vyave/nnab[iw, i]:real;
            worms[iw, i].vzave = worms[iw, i].vzave/nnab[iw, i]:real;

            KEworm[iw*i] = 0.5*worms[iw,i].m*(worms[iw,i].vx * worms[iw,i].vx + worms[iw,i].vy * worms[iw,i].vy);
            KEworm_local[iw*i] = 0.5*worms[iw,i].m*((worms[iw,i].vx-worms[iw,i].vxave) * (worms[iw,i].vx-worms[iw,i].vxave) + (worms[iw,i].vy-worms[iw,i].vyave) * (worms[iw,i].vy-worms[iw,i].vyave));
            AMworm[iw*i] = (worms[iw,i].x-hxo2)*worms[iw,i].vy - (worms[iw,i].y-hyo2)*worms[iw,i].vx; // angular momentum calculation

            nnab[iw,i] = 1;
            worms[iw,i].vxave = 0.0;
            worms[iw,i].vyave = 0.0;
            worms[iw,i].vzave = 0.0;
        }
    }
    fluid_vel(dt);
}

//FLUID FUNCTIONS

proc init_fluid_count():int {
   var numSol:int;
   if (boundary == 1) {
      var fluid_a = 2*rwall;
      var nSol_row = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
      numSol = nSol_row*nSol_row; // calcuating the number of particles in the box (exceeds rwall)
      var pos : [1..numSol,1..3] real;
      var fluid_px = fluid_a;
      var fluid_mx = 0;
      var fluid_py = fluid_a;
      var fluid_my = 0;
      var spacing = fluid_a/nSol_row;
      var row_length = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
      var row,col:real;
      for i in 1..numSol {
            row = i % row_length;
            col = ((i - row)/row_length)+1;

            pos[i,1] = fluid_mx + spacing*row + spacing; // x
            pos[i,2] = fluid_my + spacing*col; // y
            pos[i,3] = 0.0;
      }
      // mark all solvent particles that are out of bounds
      var rm_count = 0;
      for i in 1..numSol {
         var dx = pos[i,1] - hxo2;
         var dy = pos[i,2] - hyo2;
         var r = sqrt(dx*dx + dy*dy);
         if r > (0.95*rwall) {
            pos[i,3] = -1.0;
            rm_count +=1;
         }
      }
      // recounting
      writeln(rm_count);
      numSol -= rm_count;
   } else {
      halt("no other boundaries supported yet");
   }
   return numSol;
}

proc init_fluid(ref solvent: [] Structs.Particle,ref numSol: int) {
   var random_placement = false;
   // place particles in a square lattice centered in the boundary
   if (boundary == 1) {
      // circular boundary
      var fluid_a = 2*rwall;
      var nSol_row = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
      var numSol_tmp = nSol_row*nSol_row; // calcuating the number of particles in the box (exceeds rwall)
      var pos : [1..numSol_tmp,1..3] real;
      var fluid_px = fluid_a;
      var fluid_mx = 0;
      var fluid_py = fluid_a;
      var fluid_my = 0;
      var spacing = fluid_a/nSol_row;
      var row_length = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
      var row,col:real;
      for i in 1..numSol_tmp {
            row = i % row_length;
            col = ((i - row)/row_length)+1;

            pos[i,1] = fluid_mx + spacing*row + spacing; // x
            pos[i,2] = fluid_my + spacing*col; // y
            pos[i,3] = 0.0;
      }
      // mark all solvent particles that are out of bounds
      var count = 1;
      for i in 1..numSol_tmp {
         var dx = pos[i,1] - hxo2;
         var dy = pos[i,2] - hyo2;
         var r = sqrt(dx*dx + dy*dy);
         if r < (0.95*rwall) {
            solvent[count].x = pos[i,1];
            solvent[count].y = pos[i,2];
            solvent[count].z = 0.0;
            solvent[count].vx = 0.0;
            solvent[count].vy = 0.0;
            solvent[count].vz = 0.0;
            solvent[count].fx = 0.0;
            solvent[count].fy = 0.0;
            solvent[count].fz = 0.0;
            solvent[count].ptype = 2;
            solvent[count].m = 1.0;
            count +=1;
         }
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
      if (np*nworms != 10000) {
      max_x = 127.427;
      max_y = 126.69;
      }
      writeln("Worm extent for fluid:\t",max_x,"\t",max_y);
      var fluid_a = sqrt(2)*(max_x-hxo2); // box size of fluid
      var fluid_px = fluid_a + hxo2;
      var fluid_mx = fluid_a - hxo2;
      var fluid_py = fluid_a + hyo2;
      var fluid_my = fluid_a - hyo2;
      var spacing = 1.0;
      var row_length = numSol/(floor(fluid_a/(rcutsmall*spacing)):int);
      writeln("numSol ",numSol);
      writeln("Row:",row_length,"\t",row_length**2,"\t",fluid_a,"\t",fluid_a/(rcutsmall*spacing));
      if (row_length**2 > numSol) {
            writeln("fluid density too high, fixme!");
            //halt();
      }
      var row,col:real;
      for i in 1..numSol {
            row = i % row_length;
            col = ((i - row)/row_length)-1;
            solvent[i].x = hxo2-0.45*(max_x-hxo2) + spacing*rcutsmall*col + 1;
            solvent[i].y = hyo2-1.1*(max_y-hyo2) + spacing*rcutsmall*row - 2;
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
   // init macro variable arrays
   for i in 1..numSol {
      KEsol[i] = 0.0;
   }
   if (boundary == 1){
      writeln("fluid density=",numSol/(pi*rwall**2));
   } else if (boundary == 2){
      writeln("fluid density=",numSol/(6*pi*(1.5*(rwall/2))**2));
   }
   writeln("fluid init done");
   return solvent;
}

inline proc lj_thermo(i:int,j:int,r2cut_local:real) {
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

inline proc lj(i:int,j:int,r2cut_local:real) {
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
        //ffor = 48.0*sigma12*r2**(-7.0) -24*sigma6*r2**(-4.0);
        ffor = (1/sqrt(r2))**4.0;
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
        solvent[i].fx = solvent[i].fx/solvent[i].m;
        solvent[i].fy = solvent[i].fy/solvent[i].m;
        solvent[i].vx += 0.5*dt_fluid*solvent[i].fx;
        solvent[i].vy += 0.5*dt_fluid*solvent[i].fy;
        solvent[i].vz = 0.0; // just in case
        solvent[i].fx = 0.0;
        solvent[i].fy = 0.0;
        solvent[i].fz = 0.0; // just in case
    }
}

proc fluid_force_old() {
         // fluid-fluid
       var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,r,frand,gauss,fdissx,fdissy,omega:real;
         for i in 1..numSol {
             // calculating the force
             for j in (i+1)..numSol {
                 lj_thermo(i,j,r2cut);
             }
         }
         //fluid-boundary
         // var r2cutsol = r2cut
         var r2cutsol = sigma*r2cutsmall;
         forall i in 1..numSol {
             // calculate the force on the boundaries.
             for ib in 1..numPoints  {
                 lj(ib,i,r2cutsol);

             }
            //lj(1,i,r2cutsol);
         }

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
        solvent[i].fx = solvent[i].fx/solvent[i].m;
        solvent[i].fy = solvent[i].fy/solvent[i].m;
        solvent[i].vx += 0.5*dt_fluid*solvent[i].fx;
        solvent[i].vy += 0.5*dt_fluid*solvent[i].fy;
        // calculating kinetic energy here too
        KEsol[i] = 0.5*solvent[i].m*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
        AMsol[i] = (solvent[i].x-hxo2)*solvent[i].vy - (solvent[i].y-hyo2)*solvent[i].vx; // angular momentum calculation
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
            write_log(logfile,solvent[i].info());
            sudden_halt(istep);
        } else if (binid < 1) {
            write_log(logfile,solvent[i].info());
            sudden_halt(istep);
        }
        //append i to particle_list for binid
        bins[binid].scount += 1;
        bins[binid].ncount += 1;
        bins[binid].atoms.pushBack(solvent[i].id);
        bins[binid].types.pushBack(solvent[i].ptype);
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
                write_log(logfile,worms[iw,ip].info());
                sudden_halt(istep);
            } else if (binid < 1) {
                write_log(logfile,worms[iw,ip].info());
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
}

proc init_bins() {
    // writeln((4/numBins):int +1); //row
    // writeln(4%numBins); // col
    var binid,ibinnab,jbinnab,binidnbor:int;
    for ibin in 1..numBins {
        for jbin in 1..numBins {
            binid = (jbin-1)*numBins+ibin;
            //  4   3   2
            //     (*)   1
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
    var filename:string = "amatter%{010u}.xyz".format(istep);
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
                myFileWriter.writeln("A ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z);
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z);
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
    //writeln(filename,"\t",ic," lines written");
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}

proc write_xyzv(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename:string = "amatter%{010u}.xyz".format(istep);
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
                myFileWriter.writeln("A ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z," ",worms[iw,i].vx," ",worms[iw,i].vy," ",worms[iw,i].vz);
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",worms[iw,i].x," ",worms[iw,i].y," ",worms[iw,i].z," ",worms[iw,i].vx," ",worms[iw,i].vy," ",worms[iw,i].vz);
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
    write_log(logfile,filename+"\t"+ic:string+" lines written",true);
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}

proc init_macro(filename:string) {
    if C.appendFileExists(filename.c_str()) != 0 {
        var error_str:string = filename+" already exists";
        halt(error_str);
    } else {
        try {
            var datfile = open(filename, ioMode.cw);
            var myFileWriter = datfile.writer();
            myFileWriter.writeln("step \t KEworm \t KEworm-vave \t KEsol");
            datfile.fsync();
            write_log(logfile,filename+" header written");
        } catch e: Error {
            writeln(e);
        }
    }
}

proc write_macro(filename: string,istep: int) {
    var ret_int:int;
    var out_str:string = istep:string+"\t"+KEworm_total[istep]:string+"\t"+KEworm_local_total[istep]:string+"\t"+KEsol_total[istep]:string+"\t"+AMworm_total[istep]:string+"\t"+AMsol_total[istep]:string;
    ret_int = C.appendToFile(filename.c_str(),out_str.c_str());
    if ret_int != 0 {
        var err_str:string = "error in appending to "+filename;
        halt(err_str);
    }
}

proc init_log(filename:string) {
    if C.appendFileExists(filename.c_str()) != 0 {
        var error_str:string = filename+" already exists";
        halt(error_str);
    }
}

proc write_log(filename: string, out_str: string, silent: bool = false) {
    if silent == false {
      writeln(out_str);
    }
    var ret_int: int;
    ret_int = C.appendToFile(filename.c_str(),out_str.c_str());
    if ret_int != 0 {
        var err_str:string = "error in appending to "+filename;
        halt(err_str);
    }
}

proc write_params() {
    var filename:string = "params.dat";
    try {
        var paramsfile = open(filename, ioMode.cw);
        var myFileWriter = paramsfile.writer();
        myFileWriter.writeln("np\t",np.type:string,"\t",np);
        myFileWriter.writeln("nworms\t",nworms.type:string,"\t",nworms);
        myFileWriter.writeln("nsteps\t",nsteps.type:string,"\t",nsteps);
        myFileWriter.writeln("fdogic\t",fdogic.type:string,"\t",fdogic);
        myFileWriter.writeln("walldrive\t",walldrive.type:string,"\t",walldrive);
        myFileWriter.writeln("fdogicwall\t",fdogicwall.type:string,"\t",fdogicwall);
        myFileWriter.writeln("fdep\t",fdep.type:string,"\t",fdep);
        myFileWriter.writeln("fdepwall\t",fdepwall.type:string,"\t",fdepwall);
        myFileWriter.writeln("diss\t",diss.type:string,"\t",diss);
        myFileWriter.writeln("dt\t",dt.type:string,"\t",dt);
        myFileWriter.writeln("kspring\t",kspring.type:string,"\t",kspring);
        myFileWriter.writeln("kbend\t",kbend.type:string,"\t",kbend);
        myFileWriter.writeln("length0\t",length0.type:string,"\t",length0);
        myFileWriter.writeln("rcut\t",rcut.type:string,"\t",rcut);
        myFileWriter.writeln("save_interval\t",save_interval.type:string,"\t",save_interval);
        myFileWriter.writeln("boundary\t",boundary.type:string,"\t",boundary);
        myFileWriter.writeln("fluid_cpl\t",fluid_cpl.type:string,"\t",fluid_cpl);
        myFileWriter.writeln("debug\t",debug.type:string,"\t",debug);
        myFileWriter.writeln("thermo\t",thermo.type:string,"\t",thermo);
        myFileWriter.writeln("kbt\t",kbt.type:string,"\t",kbt);
        myFileWriter.writeln("sigma\t",sigma.type:string,"\t",sigma);
        myFileWriter.writeln("gamma\t",gamma.type:string,"\t",gamma);
        myFileWriter.writeln("numPoints\t",numPoints.type:string,"\t",numPoints);
        myFileWriter.writeln("numSol\t",numSol.type:string,"\t",numSol);
        myFileWriter.writeln("fluid_offset\t",fluid_offset.type:string,"\t",fluid_offset);
        myFileWriter.writeln("worm_particl_mass\t",worm_particle_mass.type:string,"\t",worm_particle_mass);
        myFileWriter.writeln("L\t",L.type:string,"\t",L);
        paramsfile.fsync();
        write_log(logfile,"params.dat written");
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

proc sudden_halt(istep:int) {
    write_xyzv(istep);
    //write_macro(nsteps);
}
