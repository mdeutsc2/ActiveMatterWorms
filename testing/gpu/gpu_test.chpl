use Math;
use Random;
use GPU;

const device = here.gpus[0];
const host = here;
config const  n = 100;
proc main() {
   var A_h : [1..n] real;
   var B_h : [1..n] real;

   foreach i in 1..n{
      B_h[i] = 2.0;
   }

   on device {
      writeln("on gpu");
      var A_d : [1..n] real;
      var B_d : [1..n] real;
      var C_d : [1..n] real;
      B_d = B_h;
      foreach i in 1..n {
         C_d[i] = 2.0;
      }

      // foreach i in 1..n {
      //    A_d[i] = B_d[i] + C_d[i];
      // }
      test(B_d,C_d,A_d); 
      on host {
         A_h = A_d;
         var A_sum = (+ reduce A_h);
         writeln(A_sum/n);
      }
   }
}

proc test(ref A:[] real, ref B:[] real,ref C:[] real) {
   @assertOnGpu foreach i in 1..n {
      C[i] = A[i] + B[i];
   }
}