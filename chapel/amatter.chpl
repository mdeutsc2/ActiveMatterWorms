use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
// configuration
config const np = 40,
             nworms = 250,
             nsteps = 250000,
             fdogic = 0.06,
             fdogicwall = 0.0,
             fdep = 1.0, // TODO: change to 4.0?
             fdepwall = 0.0,
             diss = 0.08,
             dt = 0.02, //0.02
             kspring = 57.146436,
             kbend = 40.0,
             length0 = 0.8, //particle spacing on worms
             rcut = 2.5,
             save_interval = 100,
             boundary = 1, // 1 = circle, 2 = cardiod, 3 = channel
             fluid_cpl = false;

// variables
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1.0/6.0),
      r2cutsmall = rcutsmall*rcutsmall,
      rwall = 75,//125.0*rcutsmall*sqrt(2.0),
      pi = 4.0*atan(1.0),
      twopi = 2*pi,
      pio4 = pi*0.25,
      density = nworms*np/(pi*rwall**2), //density for a circle
      hx = 2.0*rwall + 1.0,
      hy = hx,
      hyo2 = hy/2,
      hxo2 = hx/2,
      nxcell = ((hx/rcutsmall):int - 1),
      nycell = ((hy/rcutsmall):int - 1),
      dcell = hx/(nxcell:real),
      ncells = nxcell*nycell,
      nstepso500 = nsteps/500,
      gnoise = 0.80/sqrt(10.0)*0.8,
      dt2o2 = dt*dt*0.50,
      dto2 = dt*0.50,
      length2 = 2.0*length0,
      lengthmax = (length0)*((np - 1):real),
      r2inside = (rwall - rcutsmall) * (rwall-rcutsmall),
      a = 0.24, // layer spacing of worms in init_worms?
      iwalldrive = 1,
      numPoints = 840; // number of boundary points

const numTasks = here.numPUs();

var wormsDomain: domain(2) = {1..nworms,1..np};
var x : [wormsDomain] real(64);
var y : [wormsDomain] real(64);
var boundX : [1..numPoints] real(64); // x positions of boundary particles
var boundY : [1..numPoints] real(64); // y positions of boundary particles
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
var hhead: [1..ncells] int; //
var ipointto: [1..nworms*np+numPoints] int; // linked list, every particle points to another particle
var nnab: [wormsDomain] int;

var randStream = new RandomStream(real); // creating random number generator

var t = 0.0;
var ct: stopwatch, wt:stopwatch, xt:stopwatch; //calc time, io time, totaltime
//main
proc main() {
    writeln("starting...",numTasks);
    init_worms();
    write_xyz(0);
    //setting up stopwatch
    xt.start();
    for itime in 1..nsteps {
        //writeln(itime);
        t = (itime:real) *dt;

        if (itime % 100 == 0) {
	    xt.stop();
	    var total_time = xt.elapsed();
            writeln("Step: ",itime,"\t",itime/total_time,"iter/s\tCalc:",(ct.elapsed()/total_time)*100,"%\tIO:",(wt.elapsed()/total_time)*100," %\tElapsed:",total_time," s\tEst:",(nsteps/itime)*total_time," s");
	    xt.start();
        }
	ct.start();
        // first update positions and store old forces
        update_pos(itime);
        calc_forces();
        worm_wall();
        for i in 1..ncells {
            hhead[i] = -1; // 
        }
        for iw in 1..nworms*np{
            ipointto[iw] = 0;
        }
        cell_sort(itime);
        //


        //fluid_step
        
        update_vel();
	ct.stop();
        //write_xyz(itime);
        if (itime % save_interval == 0){
           wt.start();
	   write_xyz(itime);
	   wt.stop();
        }
    }
    //finalize();
    xt.stop();
    writeln("Total Time:",xt.elapsed()," s");
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
    writeln("dt\t",dt);
    //setting up the worms
    if (boundary == 1){
        // circular bound
        var equidistantThetaValues: [1..numPoints] real = 0.0;
        var deltaTheta = 2 * pi / numPoints;

        for i in 1..numPoints {
            equidistantThetaValues[i] = i * deltaTheta;
        }
        
        for i in 1..numPoints {
            boundX[i] = rwall * cos(equidistantThetaValues[i])+hxo2;
            boundY[i] = rwall * sin(equidistantThetaValues[i])+hyo2;
        }
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
                vxave[iw,i] = 0.0;
                vyave[iw,i] = 0.0;
                fx[iw,i] = 0.0;
                fy[iw,i] = 0.0;
                fxold[iw,i] = 0.0;
                fyold[iw,i] = 0.0;
            }
            thetanow += 4.0*dth;
        }
    } else if (boundary == 2) {
        // cardioid boundary
        var equidistantArcLengths: [1..numPoints] real;
        var thetaValues: [1..numPoints] real;
        var ca = rwall/2;
        var totalArcLength = 8 * ca;
        for i in 1..numPoints {
            equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
            thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
        }
        
        for i in 1..numPoints {
            boundX[i] = ca * (1 - cos(thetaValues[i])) * cos(thetaValues[i]) + hxo2;
            boundY[i] = ca * (1 - cos(thetaValues[i])) * sin(thetaValues[i]) + hyo2;
        }
        //now place worm particles
        for iw in 1..nworms {
            ireverse[iw] = 0;
            rand1 = randStream.getNext();
            if (rand1 <= 0.5) {
                ireverse(iw) = 1;
            }
            for i in 1..np {
                r = (0.9*a)*thetanow;
                dth = length0/r;
                thetanow += dth;
                x[iw,i] = hxo2 + r*cos(thetanow);
                y[iw,i] = hyo2 + r*sin(thetanow);
                xangle = atan2(y[iw,i] - hyo2, x[iw,i] - hxo2);
                //TODO give them an initial velocity going around the circle
                vx[iw,i] = 0.0;
                vy[iw,i] = 0.0;
                vxave[iw,i] = 0.0;
                vyave[iw,i] = 0.0;
                fx[iw,i] = 0.0;
                fy[iw,i] = 0.0;
                fxold[iw,i] = 0.0;
                fyold[iw,i] = 0.0;
            }
            thetanow += 4.0*dth;
        }

    } else if (boundary == 3) {
        // channel boundary
        writeln("channel boundary not implemented");
        halt();
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

proc update_pos(itime:int) {
    //forall iw in 1..nworms {
    forall iw in 1..nworms {
       foreach i in 1..np {
            x[iw, i] = x[iw,i] + vx[iw, i]*dt + fx[iw, i]*dt2o2;
            y[iw, i] = y[iw,i] + vy[iw, i]*dt + fy[iw, i]*dt2o2;
            fxold[iw, i] = fx[iw, i];
            fyold[iw, i] = fy[iw, i];
            fx[iw,i] = 0.0;
            fy[iw,i] = 0.0;
        }
    }
}

proc calc_forces() {
    var rsq:real,rand1:real,rand2:real,v1:real,v2:real,fac:real,g1:real,th:real;
    //zero out the force arrays and add Gaussian noise
    rsq = 0.0;
    for iw in 1..nworms {
        for i in 1..np {
            while ((rsq >= 0.999) || (rsq <= 0.001)) {
                rand1 = randStream.getNext();
                rand2 = randStream.getNext();
                v1 = 2.0*rand1 - 1.0;
                v2 = 2.0*rand2 - 1.0;
                rsq = v1*v1 + v2*v2;
            }
            rand1 = randStream.getNext();
            fac = sqrt(-2.0*log(rsq)/rsq);
            g1 = v1*fac*gnoise;
            th = rand1*twopi;
            fx[iw, i] = g1*cos(th) - diss*vx[iw,i];
            fy[iw, i] = g1*sin(th) - diss*vy[iw,i];
        }
    }
    //first set of springs nearest neighbor springs
    forall iw in 1..nworms {
        var ip1:int,r:real,ff:real,ffx:real,ffy:real,dx:real,dy:real;
        for i in 1..np-1 {
            ip1 = i + 1;
            dx = x[iw, ip1] - x[iw, i];
            dy = y[iw, ip1] - y[iw, i];
            r = sqrt(dx*dx + dy*dy);

            ff = -kspring*(r - length0)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            fx[iw, ip1] = fx[iw, ip1] + ffx;
            fx[iw, i] = fx[iw, i] - ffx;
            fy[iw, ip1] = fy[iw, ip1] + ffy;
            fy[iw, i] = fy[iw, i] - ffy;
        }
    }
    //bond bending terms
    forall iw in 1..nworms {
        var i3:int,i4:int,x2:real,x3:real,x4:real,y2:real,y3:real,y4:real,y23:real,y34:real,x23:real,x34:real,r23:real,r34:real,cosvalue:real;
	    var sinvalue:real,ff:real,dot:real,fac:real,f2x:real,f2y:real,f3x:real,f3y:real,f4x:real,f4y:real;
        for i2 in 1..(np-2) {
            i3 = i2 + 1;
            i4 = i2 + 2;
            //print*, i2,i3,i4
            x2 = x[iw, i2];
            y2 = y[iw, i2];
            x3 = x[iw, i3];
            y3 = y[iw, i3];
            x4 = x[iw, i4];
            y4 = y[iw, i4];
            y23 = y3 - y2;
            y34 = y4 - y3;
            x23 = x3 - x2;
            x34 = x4 - x3;
            r23 = sqrt(x23*x23 + y23*y23);
            r34 = sqrt(x34*x34 + y34*y34);

            cosvalue = abs(x23*x34 + y23*y34)/(r23*r34);
            //cosvalue = (x23*x34 + y23*y34)/(r23*r34)
            //print*, x23,x34,y23,y34,r23,r34
            if (cosvalue < 0.0) {
               cosvalue = 0.0;
               //print *, x23, x34, y23, y34, r23, r34
               //print *, "cosvalue .lt. 0.0d0"
               //stop
            }
            if (cosvalue > 1.0) {
               cosvalue = 1.0;
            }
            sinvalue = sqrt(1.0 - cosvalue*cosvalue);

            ff = -kbend*sinvalue/(r23*r34);

            dot = x23*x34 + y23*y34;
            fac = dot/(r23*r23);

            f2x = ff*(x34 - fac*x23);
            f2y = ff*(y34 - fac*y23);

            fac = dot/(r34*r34);
            f4x = ff*(fac*x34 - x23);
            f4y = ff*(fac*y34 - y23);
            f3x = -f2x - f4x;
            f3y = -f2y - f4y;

            fx[iw, i2] = fx[iw, i2] + f2x;
            fy[iw, i2] = fy[iw, i2] + f2y;

            fx[iw, i3] = fx[iw, i3] + f3x;
            fy[iw, i3] = fy[iw, i3] + f3y;

            fx[iw, i4] = fx[iw, i4] + f4x;
            fy[iw, i4] = fy[iw, i4] + f4y;
        }
    }
}

proc worm_wall_old() {
    var dx:real, dy:real, r:real, r2:real, th:real, 
           xwall:real, ywall:real, rr2:real, ffor:real, 
           dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
    var ip1:int;
    //put worm-wall interactions here, and dissipation force proportional to velocity
    for iw in 1..nworms{
        for i in 1..np{
            //dissipation proportional to v relative to local average
            // TODO swap vxave with intvx
            fx[iw, i] = fx[iw, i] - diss*(vx[iw, i] - vxave[iw, i]);
            fy[iw, i] = fy[iw, i] - diss*(vy[iw, i] - vyave[iw, i]);
            //now that we have used them, zero out vxave and vyave, recalculate below
            vxave[iw, i] = vx[iw, i];
            vyave[iw, i] = vy[iw, i];
            nnab[iw, i] = 1;
            //calculate distance to the center
            dx = x[iw, i] - hxo2;
            dy = y[iw, i] - hyo2;
            r2 = (dx*dx + dy*dy);
            //if close enough to the wall, calculate wall forces
            //use the short cut-off
            if (r2 >= r2inside) {
                //find the nearest spot on the wall
                th = atan2(dy, dx);
                xwall = hxo2 + rwall*cos(th);
                ywall = hyo2 + rwall*sin(th);
                dx = xwall - x[iw, i];
                dy = ywall - y[iw, i];
                rr2 = dx*dx + dy*dy;
                r = sqrt(rr2);
                //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
                ffor = -48.0*rr2**(-7.0) + 24.0*rr2**(-4.0);
                fx[iw, i] = fx[iw, i] + ffor*dx;
                fy[iw, i] = fy[iw, i] + ffor*dy;
               
                //iwalldrive = 1; // made this a const parameter
                //Turning on dogic drive with the wall!!!
                if (iwalldrive == 1) {
                    //first calculate unit vector along the worm
                    ip1 = i + 1;
                    if (ip1 <= np) {
                        dxi = x[iw, ip1] - x[iw, i];
                        dyi = y[iw, ip1] - y[iw, i];
                    } else {
                        dxi = x[iw, i] - x[iw, i - 1];
                        dyi = y[iw, i] - y[iw, i - 1];
                    }
                    //make it a unit vector
                    ri = sqrt(dxi*dxi + dyi*dyi);
                    dxi = dxi/ri;
                    dyi = dyi/ri;
                    //calculate the unit vector along the wall
                    dxj = -sin(th); //for cardioid, make this a vector pointing to the next wall particle
                    dyj = cos(th);

                    //if the vectors are not antiparallel, reverse the vector along the wall
                    if (dxi*dxj + dyi*dyj > 0.0) {
                        dxj = -dxj;
                        dyj = -dyj;
                    }

                    //if the two vectors have any component pointing in opposite directions
                    if (dxi*dxj + dyi*dyj < 0.0) {
                        //Find the direction for the force...
                        dx = (dxi - dxj)/2.0;
                        dy = (dyi - dyj)/2.0;
                        //normalize the direction vector
                        ri = sqrt(dx*dx + dy*dy);
                        dx = dx/ri;
                        dy = dy/ri;

                        //turn on extra-strong driving force
                        ffx = fdogicwall*dx;
                        ffy = fdogicwall*dy;
                        fx[iw, i] = fx[iw, i] + ffx;
                        fy[iw, i] = fy[iw, i] + ffy;
                    }
                }
            }
        }
    }
}

proc worm_wall() {
    //put worm-wall interactions here, and dissipation force proportional to velocity
    forall iw in 1..nworms{
        var dx:real, dy:real, r:real, r2:real, th:real, 
           xwall:real, ywall:real, rr2:real, ffor:real, 
           dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
        var ip1:int;
        for i in 1..np{
            //dissipation proportional to v relative to local average
            // TODO swap vxave with intvx
            fx[iw, i] = fx[iw, i] - diss*(vx[iw, i] - vxave[iw, i]);
            fy[iw, i] = fy[iw, i] - diss*(vy[iw, i] - vyave[iw, i]);
            //now that we have used them, zero out vxave and vyave, recalculate below
            vxave[iw, i] = vx[iw, i];
            vyave[iw, i] = vy[iw, i];
            nnab[iw, i] = 1;
            
            // calculate the force on the boundaries.
            for ib in 1..numPoints  {
                //calculate distance to the wall
                dx = x[iw,i]-boundX[ib];
                dy = y[iw,i]-boundY[ib];
                r2 = (dx*dx + dy*dy);
                //if close enough to the wall, calculate wall forces
                //use the short cut-off
                if (r2 <= r2cut) {
                    //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
                    ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                    fx[iw, i] = fx[iw, i] + ffor*dx;
                    fy[iw, i] = fy[iw, i] + ffor*dy;
                }
            }
        }
    }
}

proc cell_sort_old(itime:int) {
    var dddx:real,dddy:real,r2:real,riijj:real,ffor:real,ffx:real,ffy:real,dxi:real,dxj:real,
        ri:real,rj:real,r:real,dx:real,dy:real,dyi:real,dyj:real;
    var iworm:int,jworm:int,ip:int,jp:int,ii:int,jj:int,kk:int,i:int,ip1:int,
        jp1:int,scell:int,scnab:int,inogo:int,icnab:int,jcnab:int,icell:int,jcell:int;
    
    for iworm in 1..nworms{
        for ip in 1..np {
            ii = (iworm-1)*np+ip; //unique particle id 1<=ii<=nworms*np
            icell = 1+floor(x[iworm,ip]/dcell):int; //finds where particle is in grid
            jcell = 1+floor(y[iworm,ip]/dcell):int;
            if ((icell > nxcell) || (icell < 1)) {
                writeln("nxcell=",nxcell," icell=",icell);
                writeln(x[iworm,ip]/dcell);
                writeln(floor(x[iworm,ip]/dcell):int);
                writeln(1+floor(x[iworm,ip]/dcell):int);
                writeln("icell out of bounds\t",iworm," ",ip," ",x[iworm,ip]," ",vx[iworm,ip]," ",fx[iworm,ip]);
                write_xyz(itime);
                halt();
            }
            if ((jcell > nycell) || (jcell < 1)) {
                writeln("nxcell=",nycell," icell=",jcell);
                writeln(y[iworm,ip]/dcell);
                writeln(floor(y[iworm,ip]/dcell));
                writeln("jcell out of bounds\t",iworm," ",ip," ",x[iworm,ip]," ",vx[iworm,ip]," ",fx[iworm,ip]);
                write_xyz(itime);
                halt();
            }
            scell = icell + (jcell-1)*nxcell; // 1d-ndexing for 2d cells
            if ((scell > ncells) || (scell < 1)) {
                writeln("scell out of bounds",scell,"\t",icell,"\t",jcell);
            }
            ipointto[ii] = hhead[scell]; //
            hhead[scell] = ii;
        }
    }
    for icell in 1..nxcell {
        for jcell in 1..nycell {
            scell = icell + (jcell-1)*nxcell;
            if (hhead[scell] != -1){
                //there are particles in the cell called scell so
                //lets check all the neighbor cells
                for idir in 1..9 {
                    icnab = icell + ddx[idir];
                    if (icnab > nxcell) {break;}
                    if (icnab == 0) {break;}

                    jcnab = jcell + ddy[idir];
                    if (jcnab > nycell) {break;}
                    if (jcnab == 0) {break;}
                    scnab = icnab + (jcnab-1)*nxcell; //1d neighbor
                    if (hhead[scnab] != -1) {
                        //there are particles in the cell called scnab
                        ii = hhead[scell]; // ii is the # of the head particle
                        while (ii > 0) {

                            iworm = 1 + ((ii - 1)/np):int; // find which worm ii is in
                            ip = ii - np*(iworm - 1); // which particle in the worm is ii?
                            jj = hhead[scnab];//  head particle of neighboring cell

                            while (jj > 0) {
                                jworm = 1 + ((jj - 1)/np):int;
                                jp = jj - np*(jworm - 1);
                                inogo = 0;
                                if ((iworm == jworm) && (abs(ip-jp) <= 2)) {
                                    // on the same worm and close means no interaction calculated here
                                    inogo = 1;
                                }

                                if ((ii < jj) && (inogo == 0)) {
                                    dddx = x[jworm, jp] - x[iworm, ip];
                                    dddy = y[jworm, jp] - y[iworm, ip];
                                    r2 = dddx**2 + dddy**2;
                                    riijj = sqrt(r2);
                                    //add attractive force fdep between all pairs
                                    if (r2 <= r2cutsmall) {
                                        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
                                        ffx = ffor*dddx;
                                        ffy = ffor*dddy;
                                        fx[iworm, ip] = fx[iworm, ip] + ffx;
                                        fx[jworm, jp] = fx[jworm, jp] - ffx;
                                        fy[iworm, ip] = fy[iworm, ip] + ffy;
                                        fy[jworm, jp] = fy[jworm, jp] - ffy;

                                        //take these neighbors into account in calculating vxave and vyave
                                        vxave[iworm, ip] = vxave[iworm, ip] + vx[jworm, jp];
                                        vyave[iworm, ip] = vyave[iworm, ip] + vy[jworm, jp];
                                        nnab[iworm, ip] = nnab[iworm, ip] + 1;
                                        vxave[jworm, jp] = vxave[jworm, jp] + vx[iworm, ip];
                                        vyave[jworm, jp] = vyave[jworm, jp] + vy[iworm, ip];
                                        nnab[jworm, jp] = nnab[jworm, jp] + 1;

                                        //add 'dogic drive' to interacting pairs
                                        //first calculate unit vectors along each worm
                                        ip1 = ip + 1;
                                        if (ip1 <= np) {
                                            dxi = x[iworm, ip1] - x[iworm, ip];
                                            dyi = y[iworm, ip1] - y[iworm, ip];
                                        } else {
                                            dxi = x[iworm, ip] - x[iworm, ip - 1];
                                            dyi = y[iworm, ip] - y[iworm, ip - 1];
                                        }

                                        jp1 = jp + 1;
                                        if (jp1 <= np) {
                                            dxj = x[jworm, jp1] - x[jworm, jp];
                                            dyj = y[jworm, jp1] - y[jworm, jp];
                                        } else {
                                            dxj = x[jworm, jp] - x[jworm, jp - 1];
                                            dyj = y[jworm, jp] - y[jworm, jp - 1];
                                        }


                                        //if the two vectors have any component pointing in opposite directions
                                        if (dxi*dxj + dyi*dyj <= 0.0) {
                                            //normalize those vectors to make them unit vectors
                                            ri = sqrt(dxi*dxi + dyi*dyi);
                                            dxi = dxi/ri;
                                            dyi = dyi/ri;

                                            rj = sqrt(dxj*dxj + dyj*dyj);
                                            dxj = dxj/rj;
                                            dyj = dyj/rj;
                                            //now they are both unit vectors. Find the direction for the force...

                                            dx = (dxi - dxj)/2.0;
                                            dy = (dyi - dyj)/2.0;

                                            //normalize

                                            r = sqrt(dx*dx + dy*dy);
                                            dx = dx/r;
                                            dy = dy/r;

                                            //add an extra attractive component where kinesin drive is present

                                            ffx = fdogic*(dx) + 0.7*dddx/riijj;
                                            ffy = fdogic*(dy) + 0.7*dddy/riijj;

                                            //ffx=fdogic*(dx)
                                            //ffy=fdogic*(dy)

                                            fx[iworm, ip] = fx[iworm, ip] + ffx;
                                            fx[jworm, jp] = fx[jworm, jp] - ffx;
                                            fy[iworm, ip] = fy[iworm, ip] + ffy;
                                            fy[jworm, jp] = fy[jworm, jp] - ffy;

                                        }
                                    }
                                }
                                jj = ipointto[jj];
                                }
                            ii = ipointto[ii];
                        }
                    }
                }
            }
        }
    }
}


proc cell_sort_new(itime:int) {
    var dddx:real,dddy:real,r2:real,riijj:real,ffor:real,ffx:real,ffy:real,dxi:real,dxj:real,
        ri:real,rj:real,r:real,dx:real,dy:real,dyi:real,dyj:real;
    var iworm:int,jworm:int,ip:int,jp:int,ii:int,jj:int,kk:int,i:int,ip1:int,
        jp1:int,scell:int,scnab:int,inogo:int,icnab:int,jcnab:int,icell:int,jcell:int,
        iiboolworm:int,jjboolworm:int,ib:int,jb:int;
    
    for iworm in 1..nworms{ //sorting of worm particles into cells
        for ip in 1..np {
            ii = (iworm-1)*np+ip; //unique particle id 1<=ii<=nworms*np
            icell = 1+floor(x[iworm,ip]/dcell):int; //finds where particle is in grid
            jcell = 1+floor(y[iworm,ip]/dcell):int;
            if ((icell > nxcell) || (icell < 1)) {
                writeln("nxcell=",nxcell," icell=",icell);
                writeln(x[iworm,ip]/dcell);
                writeln(floor(x[iworm,ip]/dcell):int);
                writeln(1+floor(x[iworm,ip]/dcell):int);
                writeln("icell out of bounds\t",iworm," ",ip," ",x[iworm,ip]," ",vx[iworm,ip]," ",fx[iworm,ip]);
                write_xyz(itime);
                halt();
            }
            if ((jcell > nycell) || (jcell < 1)) {
                writeln("nycell=",nycell," jcell=",jcell);
                writeln(y[iworm,ip]/dcell);
                writeln(floor(y[iworm,ip]/dcell));
                writeln("jcell out of bounds\t",iworm," ",ip," ",y[iworm,ip]," ",vy[iworm,ip]," ",fy[iworm,ip]);
                write_xyz(itime);
                halt();
            }
            scell = icell + (jcell-1)*nxcell; // 1d-ndexing for 2d cells
            if ((scell > ncells) || (scell < 1)) {
                writeln("scell out of bounds",scell,"\t",icell,"\t",jcell);
            }
            ipointto[ii] = hhead[scell]; //
            hhead[scell] = ii;
        }
    }
    var total_worm_particles = np*nworms;
    for ib in 1..numPoints { // sorting wall particles into cells
        ii = total_worm_particles + ib; // unique id
        icell = 1 + floor(boundX[ib]/dcell):int; // finds which cell particle is in
        jcell = 1 + floor(boundY[ib]/dcell):int;
        if ((icell > nxcell) || (icell < 1)) {
                writeln("boundary particle out of cell bounds");
                writeln("nxcell=",nxcell," icell=",icell);
                halt();
            }
        if ((jcell > nycell) || (jcell < 1)) {
            writeln("boundary particle out of cell bounds");
            writeln("nycell=",nycell," jcell=",jcell);
            halt();
        }
        scell = icell + (jcell-1)*nxcell; // 1-d "street-address" for 2d cells
        if ((scell > ncells) || (scell < 1)) {
            writeln("scell out of bounds",scell,"\t",icell,"\t",jcell);
        }
        ipointto[ii] = hhead[scell];
        hhead[scell] = ii;
    }
    /*
    for icell in 1..nxcell {
       for jcell in 1..nycell {
           scell = icell + (jcell-1)*nxcell;*/

   for scell in 1..ncells {
       	icell = scell%nxcell;
        jcell = (scell-icell)/nxcell + 1;
            if (hhead[scell] != -1){
                //there are particles in the cell called scell so
                //lets check all the neighbor cells
                for idir in 1..9 {
                    icnab = icell + ddx[idir];
                    if (icnab > nxcell) {break;}
                    if (icnab == 0) {break;}

                    jcnab = jcell + ddy[idir];
                    if (jcnab > nycell) {break;}
                    if (jcnab == 0) {break;}
                    scnab = icnab + (jcnab-1)*nxcell; //1d neighbor
                    if (hhead[scnab] != -1) {
                        //there are particles in the cell called scnab
                        ii = hhead[scell]; // ii is the # of the head particle
                        while (ii > 0) {
                            if (ii > total_worm_particles) {
                                iiboolworm = 0; // is ii a worm?
                                ib = ii - total_worm_particles;
                            } else {
                                iiboolworm = 1;
                                iworm = 1 + ((ii - 1)/np):int; // find which worm ii is in
                                ip = ii - np*(iworm - 1); // which particle in the worm is ii?
                            }
                            jj = hhead[scnab];//  head particle of neighboring cell
                            while (jj > 0) {
                                if (jj >  total_worm_particles) {
                                    jjboolworm = 0; // is jj a worm?
                                    jb = jj - total_worm_particles;
                                } else {
                                    jjboolworm = 1;
                                    jworm = 1 + ((jj - 1)/np):int;
                                    jp = jj - np*(jworm - 1);
                                }
                                // worm-worm, worm-wall, wall-wall
                                if((jjboolworm == 1) && ( iiboolworm == 1)) {
                                    inogo = 0;
                                    if ((iworm == jworm) && (abs(ip-jp) <= 2)) {
                                        // on the same worm and close means no interaction calculated here
                                        inogo = 1;
                                    }

                                    if ((ii < jj) && (inogo == 0)) {
                                        dddx = x[jworm, jp] - x[iworm, ip];
                                        dddy = y[jworm, jp] - y[iworm, ip];
                                        r2 = dddx**2 + dddy**2;
                                        riijj = sqrt(r2);
                                        //add attractive force fdep between all pairs
                                        if (r2 <= r2cutsmall) {
                                            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
                                            ffx = ffor*dddx;
                                            ffy = ffor*dddy;
                                            fx[iworm, ip] = fx[iworm, ip] + ffx;
                                            fx[jworm, jp] = fx[jworm, jp] - ffx;
                                            fy[iworm, ip] = fy[iworm, ip] + ffy;
                                            fy[jworm, jp] = fy[jworm, jp] - ffy;

                                            //take these neighbors into account in calculating vxave and vyave
                                            //vxave[iworm, ip] = vxave[iworm, ip] + vx[jworm, jp];
                                            //vyave[iworm, ip] = vyave[iworm, ip] + vy[jworm, jp];
                                            //nnab[iworm, ip] = nnab[iworm, ip] + 1;
                                            //vxave[jworm, jp] = vxave[jworm, jp] + vx[iworm, ip];
                                            //vyave[jworm, jp] = vyave[jworm, jp] + vy[iworm, ip];
                                            //nnab[jworm, jp] = nnab[jworm, jp] + 1;

                                            //add 'dogic drive' to interacting pairs
                                            //first calculate unit vectors along each worm
                                            ip1 = ip + 1;
                                            if (ip1 <= np) {
                                                dxi = x[iworm, ip1] - x[iworm, ip];
                                                dyi = y[iworm, ip1] - y[iworm, ip];
                                            } else {
                                                dxi = x[iworm, ip] - x[iworm, ip - 1];
                                                dyi = y[iworm, ip] - y[iworm, ip - 1];
                                            }

                                            jp1 = jp + 1;
                                            if (jp1 <= np) {
                                                dxj = x[jworm, jp1] - x[jworm, jp];
                                                dyj = y[jworm, jp1] - y[jworm, jp];
                                            } else {
                                                dxj = x[jworm, jp] - x[jworm, jp - 1];
                                                dyj = y[jworm, jp] - y[jworm, jp - 1];
                                            }


                                            //if the two vectors have any component pointing in opposite directions
                                            if (dxi*dxj + dyi*dyj <= 0.0) {
                                                //normalize those vectors to make them unit vectors
                                                ri = sqrt(dxi*dxi + dyi*dyi);
                                                dxi = dxi/ri;
                                                dyi = dyi/ri;

                                                rj = sqrt(dxj*dxj + dyj*dyj);
                                                dxj = dxj/rj;
                                                dyj = dyj/rj;
                                                //now they are both unit vectors. Find the direction for the force...

                                                dx = (dxi - dxj)/2.0;
                                                dy = (dyi - dyj)/2.0;

                                                //normalize

                                                r = sqrt(dx*dx + dy*dy);
                                                dx = dx/r;
                                                dy = dy/r;

                                                //add an extra attractive component where kinesin drive is present

                                                ffx = fdogic*(dx) + 0.7*dddx/riijj;
                                                ffy = fdogic*(dy) + 0.7*dddy/riijj;

                                                //ffx=fdogic*(dx)
                                                //ffy=fdogic*(dy)

                                                fx[iworm, ip] = fx[iworm, ip] + ffx;
                                                fx[jworm, jp] = fx[jworm, jp] - ffx;
                                                fy[iworm, ip] = fy[iworm, ip] + ffy;
                                                fy[jworm, jp] = fy[jworm, jp] - ffy;

                                            }
                                        }
                                    }
                                } else if ((iiboolworm == 1)&&(jjboolworm == 0)) {
                                    // worm-wall (ii = worm, jj = wall)
                                    dddx = x[iworm,ip] - boundX[jb];
                                    dddy = y[iworm,ip] - boundY[jb];
                                    r2 = dddx**2 + dddy**2;
                                    riijj = sqrt(r2);
                                    //add attractive force fdep between all pairs
                                    if (r2 <= r2cutsmall) {
                                        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //fdep is attractive forceTODO: shoudl fdep = -?
                                        ffx = ffor*dddx;
                                        ffy = ffor*dddy;
                                        fx[iworm, ip] = fx[iworm, ip] + ffx;
                                        fy[iworm, ip] = fy[iworm, ip] + ffy;
                                        writeln("worm-wall");
                                    }
                                } else if ((iiboolworm == 0)&&(jjboolworm == 1)){
                                    // wall-worm (ii = wall, jj = worm)
                                    dddx = boundX[ib] - x[jworm,jp];
                                    dddy = boundY[ib] - y[jworm,jp];
                                    r2 = dddx**2 + dddy**2;
                                    riijj = sqrt(r2);
                                    if (r2 <= r2cutsmall) {
                                        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //fdep is attractive forceTODO: shoudl fdep = -?
                                        ffx = ffor*dddx;
                                        ffy = ffor*dddy;
                                        fx[jworm, jp] = fx[jworm, jp] - ffx;
                                        fy[jworm, jp] = fy[jworm, jp] - ffy;
                                        writeln("wall-worm");
                                    }
                                } else {
                                    // wall-wall interaction
                                }
                                jj = ipointto[jj];
                                }
                            ii = ipointto[ii];
                        }
                    }
                }
            }
        //}
    }
}

proc update_vel() {
    forall iw in 1..nworms {
        foreach i in 1..np {
            vx[iw, i] = vx[iw, i] + dto2*(fx[iw, i] + fxold[iw, i]);
            vy[iw, i] = vy[iw, i] + dto2*(fy[iw, i] + fyold[iw, i]);
            //vxave[iw, i] = vxave[iw, i]/nnab[iw, i];
            //vyave[iw, i] = vyave[iw, i]/nnab[iw, i];
        }
    }
}

proc write_xyz(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename:string = "amatter%{07u}.xyz".format(istep);
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer();
    myFileWriter.writeln(nworms * np + 4 + numPoints);
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
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0);
                //myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        }
    }
    for i in 1..numPoints {
        myFileWriter.writeln("I ",boundX[i]," ",boundY[i]," ",0.0);
        ic += 1;
    }
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0);
    //ic += 4;
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    writeln(filename,"\t",ic," lines written");
    } catch e: Error {
        writeln(e);
    }
}
