use GPU;
use Time;
config const nSteps = 10;
config const num = 1000;

var t1,t2: stopwatch;
var total_time=0.0;
writeln("on GPU:");
t1.start();
jacobi(here.gpus[0],num);
t1.stop();
total_time = t1.elapsed();
writeln(total_time," s");
writeln("on CPU:");
t2.start();
jacobi(here,num);
t2.stop();
total_time = t2.elapsed();
writeln(total_time," s");

proc jacobi(loc,n) {
  on loc {
    var A, B: [0..n+1] real;

    A[0] = 1; A[n+1] = 1;
    forall i in 1..n { A[i] = i:real; }

    for step in 1..nSteps {
      forall i in 1..n { B[i] = 0.33333 * (A[i-1] + A[i] + A[i+1]); }
      forall i in 1..n { A[i] = 0.33333 * (B[i-1] + B[i] + B[i+1]); }
    }
    //writeln(A);
  }
}
