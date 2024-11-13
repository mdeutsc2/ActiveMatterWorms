module Helper {
    use Amatter3d;
    use Output;
    use Math;
    use Random;
    // HELPER FUNCTIONS
inline proc gaussRand(mean: real, stddev: real): real {
    var u1 = randStream.next();
    var u2 = randStream.next();
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
   writeln("halted step # ",istep);
    write_xyzv(istep);
    //restart_write(istep);
    //write_macro(nsteps);
}
}