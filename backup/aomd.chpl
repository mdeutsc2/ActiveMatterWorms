use Random;

config const steplimit = 10,
             deltaT = 0.001,
             NDIM = 2,
             nAtom = 100;

const pi = 4.0*atan(1.0);
var density:real,kinEnergy:real,potEnergy:real,pressure:real,rCut:real,
    sKinEnergy:real, sPressure:real,sTotEnergy:real, ssKinEnergy:real,
    ssPressure:real,ssTotEnergy:real,temperature:real,totEnergy:real,
    uSum:real,virSum:real,vMag:real,vSum:real,vvSum:real;
var region: [0..NDIM+1] real;
var regionH: [0..NDIM+1] real;
var r:[1..NDIM,1..nAtom] real;
var rv:[1..NDIM,1..nAtom] real;
var ra:[1..NDIM,1..nAtom] real;
var initUcell: [0..NDIM+1] int;
var stepCount:int, moreCycles,:int, runId:int, stepAvg: int, stepEquil:int, steplimit:int;
var timeNow = 0.00 : real(64);

proc singlestep() {
    stepCount = stepCount + 1;
    timeNow = stepCount * deltaT;
    writeln("Step:",stepCount,"\t",timeNow);
    ComputeForces();
    LeapfrogStep();
    ApplyBoundaryCond();
    EvalProps();
    AccumProps(1);
    if (stepCount % setAvg == 0) {
        AccumProps(2);
        PrintSummary();
        AccumProps(0);
    }
}

proc setparams() {
    rcut = pow(2.0,1.0/6.0);
    for k in 1..NDIM {
        region[k] = initUcell[k] / sqrt(density);
        regionH[k] = 0.5*region[k];
    }
    natom = initUcell[1] * initUcell[2];
    vMag = sqrt(NDIM*(1.0-1.0/nAtom)*temperature);
}
proc setupJob() {
    InitCoords();
    InitVels();
    AccumProps(0);
    stepCount = 0;
}

proc ComputeForces() {
    var dr[NDIM+1]: real;
    var f:real,fcVal:real,rrCut:real,rr:real,rri:real,rri3:real;
    var j1:int,j2:int,k:int,n:int;
    rrCut = rcut*rcut;
    for n in 1..nAtom {
        for k in 1..NDIM) {
            ra[k,n] = 0.0;
        }
    }
    uSum = 0.0;
    virSum = 0.0;
    for j1 in 1..nAtom-1 {
        for j2 in 1..nAtom {
            rr = 0.0;
            for k in 1..NDIM {
                dr[k] = r[k,j1] - r[k,j2];
                if (abs(dr[k]) > regionH[r]) {
                    dr[k] = dr[k] - signr(region[k],dr[k]);
                }
                rr = rr + dr[k]*dr[k]
            }
            if (rr < rrCut) {
                rri - 1.0/rr;
                rri3 = rri*rri*rri;
                fcVal = 48.0*rri3*(rri3-0.5)*rri;
                for k in 1..NDIM {
                    f = fcVal * dr[k];
                    ra[k,j1] = ra[k,j1] + f;
                    ra[k,j2] = ra[k,j2] - f;
                }
                uSum = uSum + 4.0 * rri3 * (rri3-1.0) + 1.0;
                virSum = virSum + fcVal * rr;
            }
        }
    }
}

proc LeapfrogStep(){
    var k:int,n:int;
    for n in 1..nAtom {
        for k in 1..NDIM {
            rv[k,n] = rv[k,n] + deltaT * ra[k,n];
            r[k,n] = r[k,n] + deltaT * rv[k,n];
        }
    }
}

proc ApplyBoundaryCond() {
    // takes care of periodic boundary conditions
    for n in 1..nAtom {
        for k in 1..NDIM {
            if (r[k,n] >= regionH[k]) {
                r[k,n] = r[k,n] - region[k];
            } else if (r[k,n] < -regionH[k]) {
                r[k,n] = r[k,n] + region[k];
            }
        }
    }
}

proc InitCoords() {
    var c: [0..NDIM+1] real;
    var gap: [0..NDIM+1] real;
    var n:int;
    for k in 1..NDIM {
        gap[k] = region[k] /initUcell[k];
    }
    n = 0;
    for nY in 1..unitUcell[2] {
        c[2] = (nY-0.5)*gap[2] - regionH[2];
        for nX in 1.. initUcell[1] {
            c[1] = (nX - 0.5) * gap[1] - region[H];
            n += 1;
            for k in 1..NDIM {
                r[k,n] = c[k];
            }
        }
    }
}

proc InitVels() {
    var vSum: [0..NDIM+1] real;
    var randr: [1..nAtom] real;
    fillRandom(randr);
    var ang:real;
    for k in 1..NDIM {
        vSum[k] = 0;
    }
    for n in 1..nAtom {
        ang = 2.0*pi*randr[n];
        rv[1,n] = vMag * cos(ang);
        rv[2,n] = vMag * sin(ang);
        for k in 1..NDIM {
            vSum[k] = vSum[k] + rv[k,n];
        }
    }
    for k in 1..NDIM {
        vSum[k] = vSum[k]/nAtom;
    }
    for n in 1..nAtom {
        for k in 1..NDIM {
            rv[k,n] = rv[k,n] - vSum[k];
        }
    }


}

proc main() {
    writeln("starting...");
    //getnamelist();
    //printnamelist();
    setparams();
    setupjob();
    var morecycles = 1;
    while (morecycles) {
        singlestep();
        if (stepCount >= steplimit) {
            morecycles = 0;
        }
    }
}
main();