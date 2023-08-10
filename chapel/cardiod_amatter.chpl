use Math;
use Random;
use IO;
use Time; // for stopwatch
// configuration
config const np = 40,
             nworms = 1125,
             nsteps = 25000,
             fdogic = 0.06,
             fdogicwall = 0.0,
             fdep = 1.0,
             fdepwall = 0.0,
             diss = 0.08,
             dt = 0.02, //0.02
             kspring = 57.146436,
             kbend = 40.0,
             length0 = 0.8,
             rcut = 2.5,
             save_interval = 100,
             fluid_cpl = false;

// variables
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1.0/6.0),
      r2cutsmall = rcutsmall*rcutsmall,
      rwall = 125.0*rcutsmall*sqrt(2.0),
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
      a = 0.24,
      iwalldrive = 1;

const numTasks = here.numPUs();

var wormsDomain: domain(2) = {1..nworms,1..np};
var x : [wormsDomain] real(64);
var y : [wormsDomain] real(64);
var vx: [wormsDomain] real(64);
var vy: [wormsDomain] real(64);
var vxave: [wormsDomain] real(64);
var vyave: [wormsDomain] real(64);
var fx: [wormsDomain] real(64);
var fy: [wormsDomain] real(64);
var fxold: [wormsDomain] real(64);
var fyold: [wormsDomain] real(64);
var savex: [1..np] real(64);
var savey: [1..np] real(64);
var ireverse: [1..nworms] int;
var ddx: [1..9] int;
var ddy: [1..9] int;
var hhead: [1..504*504] int; //
var ipointto: [1..nworms*np] int; // linked list, every particle poitns to another particle
var nnab: [wormsDomain] int;

var boundx: [1..150] real(64);
var boundy: [1..150] real(64);

proc main() {
    var numPoints = 150;
    var a = 1.0;
    equidistantCardioidPoints(a, numPoints);

    write_xyz(0);
}

proc equidistantCardioidPoints(a: real, numPoints: int) {
    var equidistantArcLengths: [1..numPoints] real;
    var thetaValues: [1..numPoints] real;
    
    var totalArcLength = 8 * a;
    for i in 1..numPoints {
        equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
        thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
    }
    
    for i in 1..numPoints {
        boundx[i] = a * (1 - cos(thetaValues[i])) * cos(thetaValues[i]);
        boundy[i] = a * (1 - cos(thetaValues[i])) * sin(thetaValues[i]);
    }
}

proc write_xyz(istep:int) {
    var ic:int;
    var filename = "bounds" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer();
    myFileWriter.writeln("150");
    myFileWriter.writeln("# 0");
    ic = 2;
    for iw in 1..150 {
        myFileWriter.writeln("A ",boundx[iw]," ",boundy[iw]," ", 0.0);
        //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
        ic += 1;
    }
    } catch e: Error {
        writeln(e);
    }
}
