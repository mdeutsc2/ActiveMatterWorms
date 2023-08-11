var icell:int,jcell:int,scell:int;
var icell_new:int,jcell_new:int;
var nxcells = 10;
var nycells = 10;
var ncells = nxcells*nycells;

icell = 3;
jcell = 4;
scell = icell + (jcell-1)*nxcells;
writeln("scell:",scell);
icell_new = scell%nxcells;
jcell_new = (scell-icell_new)/nxcells + 1;
writeln("icell_new:",icell_new);
writeln("jcell_new:",jcell_new);


