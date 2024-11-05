use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
use CTypes; // for chpl string -> c string pointer conversion
use DynamicIters;
import Structs;


config const np=2,
            nworms = 10,
            numSol = 5,
            numPoints=5,
            fluid_cpl = true;
var wormsDomain: domain(2) = {1..nworms,1..np};
var worms: [wormsDomain] Structs.Particle;
Structs.ptc_init_counter = 1;
var solvent: [1..numSol] Structs.Particle;
Structs.ptc_init_counter = 1;
var bound: [1..numPoints] Structs.Particle;
Structs.ptc_init_counter=1;
var restartDomain: domain(1) = {1..nworms*np+numSol+numPoints};
var restart_struct: [restartDomain] Structs.Particle;
Structs.ptc_init_counter=1;

var randStream = new RandomStream(real);
for iw in 1..nworms {
    for ip in 1..np {
        worms[iw,ip].x = randStream.getNext();
        worms[iw,ip].y = randStream.getNext();
    }
}
/*
File structures
all lines starting with '#' are comment lines and contain data about the rest of the simulation
First # line is the time step of the restart
First, all worms will be written out, starting with 2 comment lines, the first containing nworms, 2nd containing np
The data will then be as following with nworms*np lines
iworm:int ip:int x:real y:real z:real vx:real vy:real vz:real
The next comment line will denote the beginning of the boundary section and contains numPoints
The data will then be as follows with numPoints lines
i:int, x:real, y:real
The next comment line will denote the beginning of the fluid section (if applicable) and contains numSol
The data will the be as follows with numSol lines
i:int x:real, y:real vx:real vy:real
*/

var filename:string = "amatter.restart";
var istep = 10;
try {
    var restart_file = open(filename, ioMode.cw);
    var restartWriter = restart_file.writer();
    // total number of particles = number of active particles + boundary + solvent (optional)
    restartWriter.writeln("# ",istep);
    restartWriter.writeln("#nworms:",nworms);
    restartWriter.writeln("#np:",np);
    for iw in 1..nworms {
        for i in 1..np {
            restartWriter.writeln(iw," ",i," ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z," ",worms[iw,i].vx," ",worms[iw,i].vy," ",worms[iw,i].vz);
        }
    }
    restartWriter.writeln("#numPoints:",numPoints);
    for i in 1..numPoints {
        restartWriter.writeln(i," ",bound[i].x," ",bound[i].y);
    }
    if fluid_cpl {
        restartWriter.writeln("#numSol:",numSol);
        for i in 1..numSol {
            restartWriter.writeln(i," ",solvent[i].x," ",solvent[i].y," ",solvent[i].vx," ",solvent[i].vy);
        }
    }
} catch e: Error {
    writeln(e);
}

try {
    var restart_file = open(filename, ioMode.r).reader();
    var read_sucess :bool;
    var linestring :string;
    //var datastring = restart_file.readString(4);
    //writeln(datastring);
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    writeln("Reading restart from step ", linestring);
    var istep_read = linestring.partition(" ")[2]:int;
    writeln(istep_read);
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    var nworms_read = linestring.partition(":")[2]:int;
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    var np_read = linestring.partition(":")[2]:int;
    writeln("nworms: ",nworms," ",nworms_read);
    writeln("np: ",np," ",np_read);
    if worms.shape != (nworms,np) {
        writeln("shape of input wrong!");
    }
    var wormsDomain2: domain(2) = {1..nworms_read,1..np_read};
    var worms2: [wormsDomain2] Structs.Particle;
    Structs.ptc_init_counter = 1;
    for iline in 1..(nworms_read*np_read) {
        read_sucess = restart_file.readLine(linestring,stripNewline=true);
        var istring : [1..8] string;
        for (i,isep) in zip(1..8,linestring.split(" ",-1)) {
            istring[i] = isep;
        }
        var iw = istring[1]:int;
        var ip = istring[2]:int;
        worms2[iw,ip].x = istring[3]:real;
        worms2[iw,ip].y = istring[4]:real;
        worms2[iw,ip].z = istring[5]:real;
        worms2[iw,ip].vx = istring[6]:real;
        worms2[iw,ip].vy = istring[7]:real;
        worms2[iw,ip].vz = istring[8]:real;
    }
    for iw in 1..nworms {
        for ip in 1..np {
            writeln(iw," ",ip," ",worms[iw,ip].x," ",worms2[iw,ip].x," ",worms[iw,ip].x==worms2[iw,ip].x," ",isClose(worms[iw,ip].x, worms2[iw,ip].x, relTol = 1e-5, absTol = 0.0));
        }
    }
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    var numPoints_read = linestring.partition(":")[2]:int;
    var bound2: [1..numPoints_read] Structs.Particle;
    Structs.ptc_init_counter=1;
    for iline in 1..numPoints_read {
        read_sucess = restart_file.readLine(linestring,stripNewline=true);
        var istring : [1..3] string;
        for (i,isep) in zip(1..3,linestring.split(" ",-1)) {
            istring[i] = isep;
        }
        var ib = istring[1]:int;
        bound2[ib].x = istring[2]:real;
        bound2[ib].y = istring[3]:real;
    }
    //assert(bound==bound2);
    if (fluid_cpl) {
        read_sucess = restart_file.readLine(linestring,stripNewline=true);
        var numSol_read = linestring.partition(":")[2]:int;
        var solvent2: [1..numSol_read] Structs.Particle;
        Structs.ptc_init_counter=1;
        for iline in 1..numSol_read {
            read_sucess = restart_file.readLine(linestring,stripNewline=true);
            var istring : [1..5] string;
            for (i,isep) in zip(1..5,linestring.split(" ",-1)) {
                istring[i] = isep;
            }
            var ip = istring[1]:int;
            solvent2[ip].x = istring[2]:real;
            solvent2[ip].y = istring[3]:real;
            solvent2[ip].vx = istring[4]:real;
            solvent2[ip].vy = istring[5]:real;
        }
        //assert(solvent==solvent2);
    }
} catch e: Error {
    writeln(e);
}
writeln("done");