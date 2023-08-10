use Math;
use Random;
use IO;
// configuration
config const np = 80,
             nworms = 1125,
             nsteps = 250,
             fdogic = 0.06,
             fdogicwall = 0.0,
             fdep = 1.0,
             fdepwall = 0.0,
             diss = 0.08,
             dt = 0.02,
             kspring = 57.146436,
             kbend = 40.0,
             length0 = 0.8,
             rcut = 2.5,
             save_interval = 10,
             fluid_cpl = false;

// variables
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1/6),
      rwall = 125.0*rcutsmall*sqrt(2.0),
      pi = 4.0*atan(1.0),
      twopi = 2*pi,
      pio4 = pi*0.25,
      density = nworms*np/(pi*rwall**2), //density for a circle
      hx = 2.0*rwall + 1.0,
      hy = hx,
      hyo2 = hy/2,
      hxo2 = hx/2,
      nxcell = ((hx/rcutsmall) - 1):int,
      nycell = ((hy/rcutsmall) - 1):int,
      dcells = hx/(nxcell:real),
      ncells = nxcell*nycell,
      nstepso500 = nsteps/500,
      gnoise = .80/sqrt(10.0)*0.8,
      dt2o2 = dt*dt*0.50,
      dto2 = dt*0.50,
      length2 = 2.0*length0,
      lengthmax = (length0)*((np - 1):real),
      r2inside = (rwall - rcutsmall) * (rwall-rcutsmall),
      a = 0.24,
      iwalldrive = 1;

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
var hhead: [1..504*504] int;
var ipointto: [1..nworms*np] int;
var nnab: [wormsDomain] int;

var randStream = new RandomStream(real);

var t = 0.0;

//main
proc main() {
    writeln("starting...");
    init_worms();
    //write_xyz(0);

    for itime in 1..nsteps {
        t = (itime:real) *dt;

        if (itime % 100 == 0) {
            writeln("Step: ",itime);
        }

        // first update positions and store old forces
        update_pos();

        //
        calc_forces();
        
        //worm_wall();

        for i in 1..ncells:int {
            hhead(i) = -1;
        }
        
        //cell_sort();

        //fluid_step
        
        //update_vel();

        if (itime % save_interval == 0){
           // write_xyz(itime);
        }
    }
    //finalize();
}

//functions
proc init_worms() {
    //array for locating neighbor cells
    var rand1:real; // temp fix for randSteam
    var r:real, dth:real, xangle:real;
    ddx[1] = 1;
    ddy[1] = 0;
    ddx[2] = 1;
    ddy[2] = 1;
    ddx[3] = 0;
    ddy[3] = 1;
    ddx[4] = -1;
    ddy[4] = 1;
    ddx[5] = -1;
    ddy[5] = 0;
    ddx[6] = -1;
    ddy[6] = -1;
    ddx[7] = 0;
    ddy[7] = -1;
    ddx[8] = 1;
    ddy[8] = -1;
    ddx[9] = 0;
    ddy[9] = 0;
    var thetanow = 5.0*pi :real;
    var rmin = a*thetanow;
    writeln("nworms\t", nworms);
    writeln("np\t", np);
    writeln("rwall\t", rwall);
    writeln("nxcells\t", nxcell, "\t nycells\t", nycell);
    writeln("density\t", density);
    //setting up the worms
    for iw in 1..nworms {
        ireverse[iw] = 0;
        rand1 = randStream.getNext();
        if (rand1 <= 0.5) {
            ireverse(iw) = 1;
        }
        for i in 1..np {
            r = a*thetanow;
            dth = length0/r;
            thetanow += dth;
            x[iw,i] = hxo2 + r*cos(thetanow);
            y[iw,i] = hyo2 + r*sin(thetanow);
            xangle = atan2(y[iw,i] - hyo2, x[iw,i] - hxo2);
            //TODO give them an initial velocity going around the circle
            vx[iw,i] = 0.0;
            vy[iw,i] = 0.0;
            fx[iw,i] = 0.0;
            fy[iw,i] = 0.0;
        }
        thetanow += 4.0*dth;
    }
    //reverse some of the worms and give them crazy colors
    for iw in 1..nworms {
        if (ireverse[iw] == 1) {
            for i in 1..np {
                savex[i] = x[iw,i];
                savey[i] = y[iw,i];
            }
            for ip in 1..np {
                x[iw,ip] = savex[np+1-ip];
                y[iw,ip] = savey[np+1-ip];
            }
        }
    }
    writeln("init done");
}

proc update_pos() {
    for iw in 1..nworms {
        for i in 1..np {
            x[iw, i] += vx[iw, i]*dt + fx[iw, i]*dt2o2;
            y[iw, i] += vy[iw, i]*dt + fy[iw, i]*dt2o2;
            fxold[iw, i] = fx[iw, i];
            fyold[iw, i] = fy[iw, i];
        }
    }
}

proc calc_forces() {
    var v1:real,v2:real,rand1:real,rand2:real,r:real,
        g1:real,fac:real,th:real;
    //zero out the force arrays and add Gaussian noise
    var rsq = 0.0;
    for iw in 1..nworms {
        for ip in 1..np {
            while(rsq >= 0.999 || rsq <= 0.001) {
                rand1 = randStream.getNext();
                rand2 = randStream.getNext();
                v1 = 2.0*rand1 - 1.0;
                v2 = 2.0*rand2 - 1.0;
                rsq = v1*v1 + v2*v2;
            }
            rand1 = randStream.getNext();
            fac = sqrt(-2.0*log(rsq)/rsq);
            g1 = v1*fac*gnoise;
            th = rand1*2*pi;
            fx[iw,ip] = g1*cos(th);
            fy[iw,ip] = g1*sin(th);
        }
    }

    //first set of springs nearest neighbor springs
    forall iw in 1..nworms {
        var dx:real,dy:real,r:real,ff:real,ffx:real,ffy:real;
        for i in 1..np-1 {
            dx = x[iw,i+1] - x[iw,i];
            dy = y[iw,i+1] - y[iw,i];
            r = sqrt(dx*dx + dy*dy);
            ff = -kspring*(r-length0)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            fx[iw,i+1] += ffx;
            fx[iw,i] -= ffx;
            fy[iw,i+1] += ffx;
            fy[iw,i] -= ffy;
        }
    }
    // bond bending terms
    forall iw in 1..nworms {
        var dx:real,dy:real,ff:real,
        ffx:real,ffy:real,r23:real,r34:real,cosvalue:real,sinvalue:real,
        f2x:real,f2y:real,f3x:real,f3y:real,f4x:real,f4y:real,dot:real,fac:real;
        var i3:int,i4:int;
        for i2 in 1..(np-2) {
            i3 = i2 + 1;
            i4 = i2 + 2;
            r23 = sqrt((x[iw,i3]-x[iw,i2])*(x[iw,i3]-x[iw,i2]) + (y[iw,i3]-y[iw,i2])*(y[iw,i3]-y[iw,i2]));
            r34 = sqrt((x[iw,i4]-x[iw,i3])*(x[iw,i4]-x[iw,i3]) + (y[iw,i4]-y[iw,i3])*(y[iw,i4]-y[iw,i3]));
            cosvalue = abs(((x[iw,i3]-x[iw,i2])*(x[iw,i4]-x[iw,i3])) + ((y[iw,i3]-y[iw,i2])*(y[iw,i4]-y[iw,i3])))/(r23*r34);
            if (cosvalue < 0.0) {
                cosvalue = 0.0;
            }
            if (cosvalue > 1.0) {
                cosvalue = 1.0;
            }
            sinvalue = sqrt(1-cosvalue*cosvalue);
            ff = -kbend*sinvalue/(r23*r34);
            dot = ((x[iw,i3]-x[iw,i2])*(x[iw,i4]-x[iw,i3])) + ((y[iw,i3]-y[iw,i2])*(y[iw,i4]-y[iw,i3]));
            
            fac = dot/(r23*r23);
            f2x = ff*((x[iw,i4]-x[iw,i3]) - fac*(x[iw,i3]-x[iw,i2]));
            f2y = ff*((y[iw,i4]-y[iw,i3]) - fac*(y[iw,i3]-y[iw,i2]));
            
            fac = dot/(r34*r34);
            f4x = ff*(fac*(x[iw,i4]-x[iw,i3]) - (x[iw,i3]-x[iw,i2]));
            f4y = ff*(fac*(y[iw,i4]-y[iw,i3]) - (y[iw,i3]-y[iw,i2]));
            f3x = -f2x - f4x;
            f3y = -f2y - f4y;

            fx[iw, i2] += f2x;
            fy[iw, i2] += f2y;

            fx[iw, i3] += f3x;
            fy[iw, i3] += f3y;

            fx[iw, i4] += f4x;
            fy[iw, i4] += f4y;
            
        }
    }
}

proc write_xyz(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename = "amatter" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer();
    myFileWriter.writeln(nworms * np + 4);
    myFileWriter.writeln("# 0");
    ic = 2;
    for iw in 1..nworms {
        dx = x[iw,1] - hxo2;
        dy = y[iw,1] - hyo2;
        xang = atan2(dy,dx);
        rx = -sin(xang);
        ry = cos(xang);
        dot = (x[iw,1] - x[iw,np])*rx + (y[iw,1] - y[iw,np])*ry;
        if (dot >= 0.0) {
            for i in 1..np {
                myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0);
                ic += 1;
            }
        }
    }
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0);
    writeln(filename,"\t",ic," lines written");
    } catch e: Error {
        writeln(e);
    }
}