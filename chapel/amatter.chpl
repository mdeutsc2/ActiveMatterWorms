//


proc initialize(numParticles: int, numSolvent: int) {
    var numParticles: int 10;
    var numSolvent: int 10;
    var N: int = numParticles + numSolvent;
    var numGridCells: int;
    var pos: [1..N*4] real;
    var vel: [1..N*4] real;
    var force: [1..N*4] real;
    var force_old: [1..N*4] real;
    var sorted_pos: [1..N*4] real;
    var sorted_vel: [1..N*4] real;
    var tangent: [1..N*4] real;
    var cellStart: [1..numGridCells] int;
    var cellEnd: [1..numGridCells] int;
}
