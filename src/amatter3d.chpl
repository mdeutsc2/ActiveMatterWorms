module Amatter3d{
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
use Output;
use BoundaryTypes;
use BinList;
use InitFluid;
use InitWorms;
use Forces;
use Helper;


const numTasks = here.numPUs();
// configuration
config const np = 80,//16,
            nworms = 800,//625,
            nsteps = 12000000    ,//00,
            fdogic = 0.06,
            walldrive = false,
            fdogicwall = 0.001,
            fdep = 0.25,
            dogic_fdep = 0.25, // extra attractive force when dogic shearing is present (originall 0.7)
            fdepwall = 6.0,
            rwall = 164,//125.0*rcutsmall*sqrt(2.0),
            dt = 0.015,
            kspring = 57.146436,
            k2spring = 50.0*kspring, //100
            k3spring = 75.0*kspring, //10.0 25.0 50.0
            kbend = 40.0,
            length0 = 0.8, //particle spacing on worms
            rcut = 2.5,
            save_interval = 4000,
            boundary = 2, // 1 = circle, 2 = cardioid, 3 = epicycloid, 4 = epigraph
            fluid_cpl = true,
            debug = false,
            thermo = true, // turn thermostat on? for solvent only
            thermow = false, // thermostat flag for worms
            kbt = 1.5,
            gamma = 6.0, // frictional constant for dissipative force (~1/damp)
            //numSol = 7000, // cardiod number of solution particles
            fluid_rho = 0.2,//8000, // disk number of solution particles
            sigma = 2.0,
            worm_particle_mass = 1.0,
            sw_epsilon = 2.0, // solvent-worm interaction epsilon
            L = 3.2, // thickness of cell
            fdrag = 0.01, // pairwise drag force for solvent-worm interaction
            restart_filename = "";

// parameters 
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1.0/6.0),
      r2cutsmall = rcutsmall*rcutsmall,
      diss = 0.002, // unused
      //pi = 4.0*atan(1.0), // defined in Math?
      twopi = 2*pi,
      pio4 = pi*0.25,
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
      inv_sqrt_dt = 1.0/sqrt(dt),
      length2 = 2.0*length0,
      lengthmax = (length0)*((np - 1):real),
      r2inside = (rwall - rcutsmall) * (rwall-rcutsmall),
      a = 0.13, // 0.12,0.18,0.48 layer spacing of worms in init_worms?
      dpd_ratio = 1.0,
      sqrt_gamma_term = sqrt(2.0*kbt*gamma),
      //numPoints = 5000,//5000,1200//589, //number of boundary points (for circle w/ r-75)
      numPoints = 9025, // given for epitrochoid of k=2, r = 41
      bdSpacing = 0.1, // spacing for epitrochoid
      fluid_offset = 2.0,//rcutsmall*sigma,//3.0; // z-offset of fluid
      io_interval = 500,
      restart_interval = save_interval*2,
      ww_epsilon = 0.5, // worm-worm interaction epsilon
      random_init = true; // use RSA to place particles




var wormsDomain: domain(2) = {1..nworms,1..np};
var worms: [wormsDomain] Structs.Particle;
Structs.ptc_init_counter = 1;
var bd = boundary_type_init(boundary);

var randStream = new randomStream(real); // creating random number generator

var numSol:int;
var fluid_area:real;
if random_init {
    if (bd.t == BD_TYPE.CARDIOID) {
        var ca = 1.5*(rwall/2);
        writeln("ca ",ca);
        writeln("cardioid area ",(6 * pi * ca ** 2));
        numSol = ceil(fluid_rho * (6 * pi * ca ** 2)):int; // Total number of particles to generate
    }
    if (bd.t == BD_TYPE.CIRCLE) {
        fluid_area = estimate_area(bd);
        writeln("disk area ",fluid_area);
        numSol = ceil(fluid_rho * fluid_area):int;
    }
    if (bd.t == BD_TYPE.EPICYCLOID) {
        var k = 2;
        var cylc_area = 12*pi*((rwall/4)+1)**2;
        writeln("epicycloid area ",cylc_area);
        numSol = ceil(fluid_rho * cylc_area):int;
    }
    if (bd.t == BD_TYPE.EPITROCHOID) {
        var cylc_area = 12*pi*((rwall/4)+1)**2;//estimate_area(bd);
        writeln("epicycloid area ",cylc_area);
        numSol = ceil(fluid_rho * cylc_area):int;
    }
   writeln("RSA numSol ",numSol);

} else {
   numSol = init_fluid_count();
   writeln("Array numSol ",numSol);
}


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

var restart_timestep = 0; // timestep read in from restart file
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
      if random_init {
        if (bd.t == BD_TYPE.CARDIOID) {
            solvent = init_fluid_rsa2(solvent,numSol);
        }
        if (bd.t == BD_TYPE.CIRCLE) {
            solvent = init_fluid_rsa1(solvent,numSol);
        }
        if (bd.t == BD_TYPE.EPICYCLOID) {
            solvent = init_fluid_rsa3(solvent,numSol);
        }
        if (bd.t == BD_TYPE.EPITROCHOID) { // epitrochoid
            solvent = init_fluid_rsa4(solvent,numSol);
        }
      } else {
        solvent = init_fluid(solvent, numSol);
      }
    }
    //populate the bins with lists of atoms
    update_cells(0); //again after fluid // somehow segfaults here?
    write_log(logfile,"numSol\t"+numSol:string);
    //write_log(logfile,"fluid equilibrated...5000dt");

    var macro_filename:string = "energies.dat";
    init_macro(macro_filename);
    if (restart_filename.isEmpty() == false) {
        //restart_read("amatter.restart")
        restart_timestep = restart_read(restart_filename);
        write_xyzv(0+restart_timestep);
        restart_write(0+restart_timestep);
        write_log(logfile,"simulation loaded from "+restart_filename+" at itime="+restart_timestep:string);
    } else {
        write_xyzv(0);
        restart_write(0);
    }
    
    //setting up stopwatch
    var log_str:string;
    xt.start();
    for itime in 1..nsteps {
        //writeln(itime);
        t = (itime:real) *dt;

        if (itime % io_interval == 0) {
            xt.stop();
            total_time += xt.elapsed();
            var out_str:string = "Step: "+(itime+restart_timestep):string+"\t"+
                      (io_interval/xt.elapsed()):string+"iter/s\tCalc:"+
                      ((ct.elapsed()/total_time)*100):string+"%\tElapsed:"+
                      total_time:string+" s\t Est:"+
                      convertSeconds(((nsteps-itime)*(total_time/itime)):int)+" s";
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
        wt.start();
        if (itime % nstepso1e5 == 0){
            write_macro(macro_filename,itime);
        }
        if (itime % save_interval == 0){
            write_xyzv(itime+restart_timestep);
            if (itime % restart_interval == 0) {
                restart_write(itime+restart_timestep);
            }
        }
        wt.stop();
    }
    //finalize();
    xt.stop();
    //write_macro(nsteps);
    write_log(logfile,"Total Time:"+total_time:string+" s");
}
} // end module Amatter3d
