module Worms {
    use Parameters;
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
                bound[i].x = rwall * cos(equidistantThetaValues[i])+hxo2;
                bound[i].y = rwall * sin(equidistantThetaValues[i])+hyo2;
                bound[i].z = 0.0;
                bound[i].ptype = 3;
            }
            for iw in 1..nworms {
                ireverse[iw] = 0;
                rand1 = randStream.getNext();
                if (rand1 <= 0.5) {
                    ireverse[iw] = 1;
                }
                for i in 1..np {
                    r = a*thetanow;
                    dth = length0/r;
                    thetanow += dth;
                    worms[iw,i].x = hxo2 + r*cos(thetanow);
                    worms[iw,i].y = hyo2 + r*sin(thetanow);
                    xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                    //TODO give them an initial velocity going around the circle
                    worms[iw,i].ptype = 1;
                    // vx[iw,i] = 0.0; // NOTE: this is all set by initializer
                    // vy[iw,i] = 0.0;
                    // vxave[iw,i] = 0.0;
                    // vyave[iw,i] = 0.0;
                    // fx[iw,i] = 0.0;
                    // fy[iw,i] = 0.0;
                    // fxold[iw,i] = 0.0;
                    // fyold[iw,i] = 0.0;
                }
                thetanow += 4.0*dth;
            }
        } else if (boundary == 2) {
            // cardioid boundary
            var equidistantArcLengths: [1..numPoints] real;
            var thetaValues: [1..numPoints] real;
            var ca = 1.5*(rwall/2);
            var totalArcLength = 8 * ca;
            for i in 1..numPoints {
                equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
                thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
            }

            for i in 1..numPoints {
                bound[i].x = ca * (1 - cos(thetaValues[i])) * cos(thetaValues[i]) + hxo2 + ca;
                bound[i].y = ca * (1 - cos(thetaValues[i])) * sin(thetaValues[i]) + hyo2;
            }
            //now place worm particles
            for iw in 1..nworms {
                ireverse[iw] = 0;
                rand1 = randStream.getNext();
                if (rand1 <= 0.5) {
                    ireverse[iw] = 1;
                }
                for i in 1..np {
                    r = (0.65*a)*thetanow;
                    dth = length0/r;
                    thetanow += dth;
                    worms[iw,i].x = hxo2 + r*cos(thetanow);
                    worms[iw,i].y = hyo2 + r*sin(thetanow);
                    xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                    //TODO give them an initial velocity going around the circle
                    worms[iw,i].ptype = 1;
                    // vx[iw,i] = 0.0;
                    // vy[iw,i] = 0.0;
                    // vxave[iw,i] = 0.0;
                    // vyave[iw,i] = 0.0;
                    // fx[iw,i] = 0.0;
                    // fy[iw,i] = 0.0;
                    // fxold[iw,i] = 0.0;
                    // fyold[iw,i] = 0.0;
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
                    savex[i] = worms[iw,i].x;
                    savey[i] = worms[iw,i].y;
                }
                for ip in 1..np {
                    worms[iw,ip].x = savex[np+1-ip];
                    worms[iw,ip].y = savey[np+1-ip];
                }
            }
        }

        //init macro variable arrays
        for i in 1..nworms*np{
            KEworm[i] = 0.0;
        }
        writeln("init done");
    }

    proc update_pos(itime:int) {
        //forall iw in 1..nworms {
        forall iw in 1..nworms {
        foreach i in 1..np {
                worms[iw,i].x = worms[iw,i].x + worms[iw,i].vx*dt + worms[iw, i].fx*dt2o2;
                worms[iw,i].y = worms[iw,i].y + worms[iw,i].vy*dt + worms[iw, i].fy*dt2o2;
                worms[iw,i].fxold = worms[iw,i].fx;
                worms[iw,i].fyold = worms[iw,i].fy;
                worms[iw,i].fx = 0.0;
                worms[iw,i].fy = 0.0;
                worms[iw,i].fz = 0.0;
            }
        }
    }

    proc calc_forces() {
        var rsq:real,rand1:real,rand2:real,v1:real,v2:real,fac:real,g1:real,th:real;
        //zero out the force arrays and add Gaussian noise
        rsq = 0.0;
        /* for iw in 1..nworms {
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
                worms[iw,i].fx = g1*cos(th) - diss*worms[iw,i].vx;
                worms[iw,i].fy = g1*sin(th) - diss*worms[iw,i].vy;
            }
        } */
        //first set of springs nearest neighbor springs
        forall iw in 1..nworms {
            var ip1:int,r:real,ff:real,ffx:real,ffy:real,dx:real,dy:real;
            for i in 1..np-1 {
                ip1 = i + 1;
                dx = worms[iw,ip1].x - worms[iw,i].x;
                dy = worms[iw,ip1].y - worms[iw,i].y;
                r = sqrt(dx*dx + dy*dy);

                ff = -kspring*(r - length0)/r;
                ffx = ff*dx;
                ffy = ff*dy;
                worms[iw,ip1].fx += ffx;
                worms[iw,i].fx -= ffx;
                worms[iw,ip1].fy += ffy;
                worms[iw,i].fy -= ffy;
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
                x2 = worms[iw,i2].x;
                y2 = worms[iw,i2].y;
                x3 = worms[iw, i3].x;
                y3 = worms[iw, i3].y;
                x4 = worms[iw, i4].x;
                y4 = worms[iw, i4].y;
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

                worms[iw, i2].fx += f2x;
                worms[iw, i2].fy += f2y;

                worms[iw, i3].fx += f3x;
                worms[iw, i3].fy += f3y;

                worms[iw, i4].fx += f4x;
                worms[iw, i4].fy += f4y;
            }
        }
    }

    proc worm_wall() {
        var dx:real, dy:real, r:real, r2:real, th:real,
            xwall:real, ywall:real, rr2:real, ffor:real,
            dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
        var ip1:int;
        //put worm-wall interactions here, and dissipation force proportional to velocity
        for iw in 1..nworms{
            for i in 1..np{
                //dissipation proportional to v relative to local average
                // TODO swap vxave with intvx
                worms[iw,i].fx = worms[iw,i].fx - diss*(worms[iw,i].vx - worms[iw,i].vxave);
                worms[iw,i].fy = worms[iw,i].fy - diss*(worms[iw,i].vy - worms[iw,i].vyave);
                //now that we have used them, zero out vxave and vyave, recalculate below
                worms[iw,i].vxave = worms[iw,i].vx;
                worms[iw,i].vyave = worms[iw,i].vy;
                nnab[iw, i] = 1;
                //calculate distance to the center
                dx = worms[iw, i].x - hxo2;
                dy = worms[iw, i].y - hyo2;
                r2 = (dx*dx + dy*dy);
                //if close enough to the wall, calculate wall forces
                //use the short cut-off
                if (r2 >= r2inside) {
                    //find the nearest spot on the wall
                    th = atan2(dy, dx);
                    xwall = hxo2 + rwall*cos(th);
                    ywall = hyo2 + rwall*sin(th);
                    dx = xwall - worms[iw,i].x;
                    dy = ywall - worms[iw,i].y;
                    rr2 = dx*dx + dy*dy;
                    r = sqrt(rr2);
                    //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
                    ffor = -48.0*rr2**(-7.0) + 24.0*rr2**(-4.0);
                    worms[iw, i].fx += ffor*dx;
                    worms[iw, i].fy += ffor*dy;

                    //iwalldrive = 1; // made this a const parameter
                    //Turning on dogic drive with the wall!!!
                    if (iwalldrive == 1) {
                        //first calculate unit vector along the worm
                        ip1 = i + 1;
                        if (ip1 <= np) {
                            dxi = worms[iw, ip1].x - worms[iw, i].x;
                            dyi = worms[iw, ip1].y - worms[iw, i].y;
                        } else {
                            dxi = worms[iw, i].x - worms[iw, i - 1].x;
                            dyi = worms[iw, i].y - worms[iw, i - 1].y;
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
                            worms[iw,i].fx += ffx;
                            worms[iw,i].fy += ffy;
                        }
                    }
                }
            }
        }
    }

    proc worm_wall_new() {
        //put worm-wall interactions here, and dissipation force proportional to velocity
        forall iw in 1..nworms{
            var dx:real, dy:real, r:real, r2:real, th:real,
            xwall:real, ywall:real, rr2:real, ffor:real,
            dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
            var ip1:int;
            for i in 1..np{
                //dissipation proportional to v relative to local average
                // TODO swap vxave with intvx
                //worms[iw,i].fx = worms[iw,i].fx - diss*(worms[iw,i].vx - worms[iw,i].vxave);
                //worms[iw,i].fy = worms[iw,i].fy - diss*(worms[iw,i].vy - worms[iw,i].vyave);
                //now that we have used them, zero out vxave and vyave, recalculate below
                //worms[iw,i].vxave = worms[iw,i].vx;
                //worms[iw,i].vyave = worms[iw,i].vy;
                nnab[iw, i] = 1;

                // calculate the force on the boundaries.
                for ib in 1..numPoints  {
                    //calculate distance to the wall
                    dx = worms[iw,i].x-bound[ib].x;
                    dy = worms[iw,i].y-bound[ib].y;
                    r2 = (dx*dx + dy*dy);
                    //if close enough to the wall, calculate wall forces
                    //use the short cut-off
                    if (r2 <= r2cut) {
                        //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
                        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                        worms[iw,i].fx += ffor*dx;
                        worms[iw,i].fy += ffor*dy;
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
                icell = 1+floor(worms[iworm,ip].x/dcell):int; //finds where particle is in grid
                jcell = 1+floor(worms[iworm,ip].y/dcell):int;
                if ((icell > nxcell) || (icell < 1)) {
                    writeln("nxcell=",nxcell," icell=",icell);
                    writeln(worms[iworm,ip].x/dcell);
                    writeln(floor(worms[iworm,ip].x/dcell):int);
                    writeln(1+floor(worms[iworm,ip].x/dcell):int);
                    writeln("icell out of bounds\t",iworm," ",ip," ",worms[iworm,ip].x," ",worms[iworm,ip].vx," ",worms[iworm,ip].fx);
                    write_xyz(itime);
                    halt();
                }
                if ((jcell > nycell) || (jcell < 1)) {
                    writeln("nxcell=",nycell," icell=",jcell);
                    writeln(worms[iworm,ip].y/dcell);
                    writeln(floor(worms[iworm,ip].y/dcell));
                    writeln("jcell out of bounds\t",iworm," ",ip," ",worms[iworm,ip].y," ",worms[iworm,ip].vx," ",worms[iworm,ip].fx);
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
                                        dddx = worms[jworm, jp].x - worms[iworm, ip].x;
                                        dddy = worms[jworm, jp].y - worms[iworm, ip].y;
                                        r2 = dddx**2 + dddy**2;
                                        riijj = sqrt(r2);
                                        //add attractive force fdep between all pairs
                                        if (r2 <= r2cutsmall) {
                                            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
                                            ffx = ffor*dddx;
                                            ffy = ffor*dddy;
                                            worms[iworm,ip].fx += ffx;
                                            worms[jworm,jp].fx -= ffx;
                                            worms[iworm,ip].fy += ffy;
                                            worms[jworm,jp].fy -= ffy;

                                            //take these neighbors into account in calculating vxave and vyave
                                            //worms[iworm,ip].vxave += worms[jworm,jp].vx;
                                            //worms[iworm,ip].vyave += worms[jworm,jp].vy;
                                            nnab[iworm, ip] = nnab[iworm, ip] + 1;
                                            //worms[jworm,jp].vxave += worms[iworm,ip].vx;
                                            //worms[jworm,jp].vyave += worms[iworm,ip].vy;
                                            nnab[jworm, jp] = nnab[jworm, jp] + 1;

                                            //add 'dogic drive' to interacting pairs
                                            //first calculate unit vectors along each worm
                                            ip1 = ip + 1;
                                            if (ip1 <= np) {
                                                dxi = worms[iworm,ip1].x - worms[iworm, ip].x;
                                                dyi = worms[iworm,ip1].y - worms[iworm, ip].y;
                                            } else {
                                                dxi = worms[iworm, ip].x - worms[iworm, ip - 1].x;
                                                dyi = worms[iworm, ip].y - worms[iworm, ip - 1].y;
                                            }

                                            jp1 = jp + 1;
                                            if (jp1 <= np) {
                                                dxj = worms[jworm, jp1].x - worms[jworm, jp].x;
                                                dyj = worms[jworm, jp1].y - worms[jworm, jp].x;
                                            } else {
                                                dxj = worms[jworm, jp].x - worms[jworm, jp - 1].x;
                                                dyj = worms[jworm, jp].y - worms[jworm, jp - 1].y;
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

                                                worms[iworm,ip].fx += ffx;
                                                worms[jworm,jp].fx -= ffx;
                                                worms[iworm,ip].fy += ffy;
                                                worms[jworm,jp].fy -= ffy;
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
                icell = 1+floor(worms[iworm,ip].x/dcell):int; //finds where particle is in grid
                jcell = 1+floor(worms[iworm,ip].y/dcell):int;
                if ((icell > nxcell) || (icell < 1)) {
                    writeln("nxcell=",nxcell," icell=",icell);
                    writeln(worms[iworm,ip].x/dcell);
                    writeln(floor(worms[iworm,ip].x/dcell):int);
                    writeln(1+floor(worms[iworm,ip].x/dcell):int);
                    writeln("icell out of bounds\t",iworm," ",ip," ",worms[iworm,ip].x," ",worms[iworm,ip].vx," ",worms[iworm,ip].fx);
                    write_xyz(itime);
                    halt();
                }
                if ((jcell > nycell) || (jcell < 1)) {
                    writeln("nycell=",nycell," jcell=",jcell);
                    writeln(worms[iworm,ip].y/dcell);
                    writeln(floor(worms[iworm,ip].y/dcell));
                    writeln("jcell out of bounds\t",iworm," ",ip," ",worms[iworm,ip].y," ",worms[iworm,ip].vy," ",worms[iworm,ip].fy);
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
            icell = 1 + floor(bound[ib].x/dcell):int; // finds which cell particle is in
            jcell = 1 + floor(bound[ib].y/dcell):int;
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
                                            dddx = worms[jworm,jp].x - worms[iworm, ip].x;
                                            dddy = worms[jworm,jp].y - worms[iworm, ip].y;
                                            r2 = dddx**2 + dddy**2;
                                            riijj = sqrt(r2);
                                            //add attractive force fdep between all pairs
                                            if (r2 <= r2cutsmall) {
                                                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
                                                ffx = ffor*dddx;
                                                ffy = ffor*dddy;
                                                worms[iworm, ip].fx += ffx;
                                                worms[jworm, jp].fx -= ffx;
                                                worms[iworm, ip].fy += ffy;
                                                worms[jworm, jp].fy -= ffy;

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
                                                    dxi = worms[iworm, ip1].x - worms[iworm, ip].x;
                                                    dyi = worms[iworm, ip1].y - worms[iworm, ip].y;
                                                } else {
                                                    dxi = worms[iworm, ip].x - worms[iworm, ip - 1].x;
                                                    dyi = worms[iworm, ip].y - worms[iworm, ip - 1].y;
                                                }

                                                jp1 = jp + 1;
                                                if (jp1 <= np) {
                                                    dxj = worms[jworm, jp1].x - worms[jworm, jp].x;
                                                    dyj = worms[jworm, jp1].y - worms[jworm, jp].y;
                                                } else {
                                                    dxj = worms[jworm, jp].x - worms[jworm, jp - 1].x;
                                                    dyj = worms[jworm, jp].y - worms[jworm, jp - 1].y;
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
                                                    worms[iworm, ip].fx += ffx;
                                                    worms[jworm, jp].fx -= ffx;
                                                    worms[iworm, ip].fy += ffy;
                                                    worms[jworm, jp].fy -= ffy;

                                                }
                                            }
                                        }
                                    } else if ((iiboolworm == 1)&&(jjboolworm == 0)) {
                                        // worm-wall (ii = worm, jj = wall)
                                        dddx = worms[iworm,ip].x - bound[jb].x;
                                        dddy = worms[iworm,ip].y - bound[jb].y;
                                        r2 = dddx**2 + dddy**2;
                                        riijj = sqrt(r2);
                                        //add attractive force fdep between all pairs
                                        if (r2 <= r2cutsmall) {
                                            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //fdep is attractive forceTODO: shoudl fdep = -?
                                            ffx = ffor*dddx;
                                            ffy = ffor*dddy;
                                            worms[iworm,ip].fx += ffx;
                                            worms[iworm,ip].fy += ffy;
                                            //writeln("worm-wall");
                                        }
                                    } else if ((iiboolworm == 0)&&(jjboolworm == 1)){
                                        // wall-worm (ii = wall, jj = worm)
                                        dddx = bound[ib].x - worms[jworm,jp].x;
                                        dddy = bound[ib].y - worms[jworm,jp].y;
                                        r2 = dddx**2 + dddy**2;
                                        riijj = sqrt(r2);
                                        if (r2 <= r2cutsmall) {
                                            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj; //fdep is attractive forceTODO: shoudl fdep = -?
                                            ffx = ffor*dddx;
                                            ffy = ffor*dddy;
                                            worms[jworm,jp].fx += ffx;
                                            worms[jworm,jp].fy += ffy;
                                            //writeln("wall-worm");
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
                worms[iw, i].vx += dto2*(worms[iw,i].fx + worms[iw,i].fxold);
                worms[iw, i].vy += dto2*(worms[iw,i].fy + worms[iw,i].fyold);
                //vxave[iw, i] = vxave[iw, i]/nnab[iw, i];
                //vyave[iw, i] = vyave[iw, i]/nnab[iw, i];
                KEworm[iw*i] = 0.5*(worms[iw,i].vx * worms[iw,i].vx + worms[iw,i].vy * worms[iw,i].vy);
            }
        }
    }
}