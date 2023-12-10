use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
use CTypes; // for chpl string -> c string pointer conversion

extern {
        #include <stdlib.h>
        #include <math.h>

        // static int test(const char *teststring) {
        //     printf("%s\n",teststring);
        //     return 0; // File does not exist
        // }

        static int acctest(int n, double A[n][n], double Anew[n][n]) {
            printf("in acctest\n");
            double tol = 0.0005;
            int iter_max = 10000;
            double error = 1.0;
            int iter = 0;
            while (error > tol && iter < iter_max) {
                error = 0.0;
                #pragma omp parallel for reduction(max:error)
                for(int j=1; j<n-1; j++) {
                    for (int i=1; i<n-1; i++) {
                        Anew[j][i] = 0.25*(A[j][i+1] + A[j][i-1] + A[j-1][i] + A[j+1][i]);
                        error = fmax(error,fabs(Anew[j][i]-A[j][i]));
                    }
                }
                #pragma omp parallel for
                for(int j=1; j<n-1; j++) {
                    for (int i=1; i<n-1; i++) {
                        A[j][i] = Anew[j][i];
                    }
                }
                if(iter % 100 == 0) printf("%5d, %0.6f\n", iter, error);
                //printf("%5d, %0.6f\n", iter, error);
                iter++;
            }
            return 0;
        }
    }

extern proc acctest(n:int, A: [] real, Anew: [] real): int;

proc main() {
    var ct: stopwatch;
    var teststring:string;
    var n:int;
    n = 1000;
    var A: [1..n,1..n] real;
    var Anew: [1..n,1..n] real;
    fillRandom(A);
    // teststring = "hello from c";
    // test(teststring.c_str());
    ct.start();
    acctest(n,A,Anew);
    ct.stop();
    writeln(ct.elapsed():string+" s");
}