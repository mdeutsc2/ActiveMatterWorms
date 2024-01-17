use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
use CTypes; // for chpl string -> c string pointer conversion
use YAML;
use RecordParser;
use DynamicIters;
import Structs;


config const np=2,
            nworms = 10,
            numSol = 5,
            numPoints=5;
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

// copy data into restart data array
// copy worms
var restart_id = 1;
for id in 1..nworms*np {
    var iw,ip : int;
    iw = 1 + ((id - 1)/np):int;
    ip = id - np*(iw - 1);
    restart_struct[id].id = worms[iw,ip].id;
    restart_struct[id].x = worms[iw,ip].x;// position
    restart_struct[id].y = worms[iw,ip].y;
    restart_struct[id].z = worms[iw,ip].z;
    restart_struct[id].vx = worms[iw,ip].vx;//velocity
    restart_struct[id].vy = worms[iw,ip].vy;
    restart_struct[id].vz = worms[iw,ip].vz;
    restart_struct[id].ptype = 1;//worms[iw,ip].ptype;
}
restart_id = nworms*np;
// copy solvent
for id in 1..numSol {
    restart_struct[restart_id+id].id = solvent[id].id;
    restart_struct[restart_id+id].x = solvent[id].x; //position
    restart_struct[restart_id+id].y = solvent[id].y;
    restart_struct[restart_id+id].z = solvent[id].z;
    restart_struct[restart_id+id].vx = solvent[id].vx;// velocitiy
    restart_struct[restart_id+id].vy = solvent[id].vy;
    restart_struct[restart_id+id].vz = solvent[id].vz;
    restart_struct[restart_id+id].ptype = 2;//solvent[id].ptype;
}
restart_id += numSol;
// copy numPoints
for id in 1..numPoints {
    restart_struct[restart_id+id].id = bound[id].id;
    restart_struct[restart_id+id].x = bound[id].x; //position
    restart_struct[restart_id+id].y = bound[id].y;
    restart_struct[restart_id+id].z = bound[id].z;
    restart_struct[restart_id+id].vx = bound[id].vx;// velocitiy
    restart_struct[restart_id+id].vy = bound[id].vy;
    restart_struct[restart_id+id].vz = bound[id].vz;
    restart_struct[restart_id+id].ptype = 3;//bound[id].ptype;
}
// write data into file
try {
    var restart_file = open("restart.yaml", ioMode.cw);
    for id in restartDomain {
        restart_file.writer().withSerializer(new yamlSerializer()).write(restart_struct);   
    }
} catch e: Error {
    writeln(e);
}
writeln("written to restart.yaml");

for id in restartDomain{
    restart_struct[id].id = 0;
    restart_struct[id].ptype=-1;
}
//read restart data
// try {
//     var restart_file = open("restart.yaml", ioMode.r);
//     for id in restartDomain {
//         restart_struct[id] = restart_file.reader().withDeserializer(new yamlDeserializer()).read(Structs.Particle);
//     }
//     //var restart_struct2 = restart_file.reader().withDeserializer(new yamlDeserializer()).read(Structs.Particle);
// } catch e: Error {
//     writeln(e);
// }

// check that restart data is the same size

// separate restart array into worms, solvent, and boundary particle arrays


// recordParser
var 
