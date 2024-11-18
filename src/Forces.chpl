module Forces {
    use Amatter3d;
    use Random;
    use DynamicIters;
    // user-defined modules
    import Structs;
    use Helper;

proc update_pos(itime:int) {
    //forall iw in 1..nworms {
    forall iw in 1..nworms {
       forall i in 1..np {
            worms[iw,i].fx = worms[iw,i].fx/worms[iw,i].m;
            worms[iw,i].fy = worms[iw,i].fy/worms[iw,i].m;
            worms[iw,i].fz = worms[iw,i].fz/worms[iw,i].m;

            //dissipation proportional to v relative to local average
            //worms[iw,i].fx = worms[iw,i].fx - diss*(worms[iw,i].vx - worms[iw,i].vxave);
            //worms[iw,i].fy = worms[iw,i].fy - diss*(worms[iw,i].vy - worms[iw,i].vyave);

            worms[iw,i].x = worms[iw,i].x + worms[iw,i].vx*dt + worms[iw, i].fx*dt2o2;
            worms[iw,i].y = worms[iw,i].y + worms[iw,i].vy*dt + worms[iw, i].fy*dt2o2;
            worms[iw,i].z = worms[iw,i].z + worms[iw,i].vz*dt + worms[iw, i].fz*dt2o2;

            worms[iw,i].fxold = worms[iw,i].fx;
            worms[iw,i].fyold = worms[iw,i].fy;
            worms[iw,i].fzold = worms[iw,i].fz;

            worms[iw,i].fx = 0.0;
            worms[iw,i].fy = 0.0;
            worms[iw,i].fz = 0.0;

            // periodic boundary conditions top and bottom for worms
            //periodic boundary conditions
            if worms[iw,i].z > L {
                //worms[iw,i].z = worms[iw,i].z - L;
                worms[iw,i].z = L;
            } else if worms[iw,i].z < 0.0 {
                //worms[iw,i].z = worms[iw,i].z + L;
                worms[iw,i].z = 0.0;
            }
        }
    }
    fluid_pos(dt);
}

proc intraworm_forces() {
    //first set of springs nearest neighbor springs
    forall iw in 1..nworms {
        var ip1:int,r:real,ff:real,ffx:real,ffy:real,ffz:real,dx:real,dy:real,dz:real;
        for i in 1..np-1 {
            ip1 = i + 1;
            dx = worms[iw,ip1].x - worms[iw,i].x;
            dy = worms[iw,ip1].y - worms[iw,i].y;
            dz = worms[iw,ip1].z - worms[iw,i].z;
            r = sqrt(dx*dx + dy*dy + dz*dz);

            ff = -kspring*(r - length0)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            ffz = ff*dz;
            worms[iw,ip1].fx += ffx;
            worms[iw,i].fx -= ffx;
            worms[iw,ip1].fy += ffy;
            worms[iw,i].fy -= ffy;
            worms[iw,ip1].fz += ffz;
            worms[iw,i].fz -= ffz;
        }
    }
    //bond bending terms
    forall iw in 1..nworms{
        var ip2:int,r:real,ff:real,ffx:real,ffy:real,ffz:real,dx:real,dy:real,dz:real;
        for i in 1..(np-2) { 
            ip2 = i + 2;
            dx = worms[iw,ip2].x - worms[iw,i].x;
            dy = worms[iw,ip2].y - worms[iw,i].y;
            dz = worms[iw,ip2].z - worms[iw,i].z;
            r = sqrt(dx*dx + dy*dy + dz*dz);

            ff = -k2spring*(r - length2)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            ffz = ff*dz;

            worms[iw,ip2].fx += ffx;
            worms[iw,i].fx -= ffx;

            worms[iw,ip2].fy += ffy;
            worms[iw,i].fy -= ffy;

            worms[iw,ip2].fz += ffz;
            worms[iw,i].fz -= ffz;
        }
    }
    // 3-spring bond-bending
    var length3 = 3.0*length0;
    forall iw in 1..nworms {
        var ip3:int,r:real,ff:real,ffx:real,ffy:real,ffz:real,dx:real,dy:real,dz:real;
        for i in 1..(np-3) { 
            ip3 = i + 3;
            dx = worms[iw,ip3].x - worms[iw,i].x;
            dy = worms[iw,ip3].y - worms[iw,i].y;
            dz = worms[iw,ip3].z - worms[iw,i].z;
            r = sqrt(dx*dx + dy*dy + dz*dz);

            ff = -k3spring*(r - length3)/r;
            ffx = ff*dx;
            ffy = ff*dy;
            ffz = ff*dz;

            worms[iw,ip3].fx += ffx;
            worms[iw,i].fx -= ffx;

            worms[iw,ip3].fy += ffy;
            worms[iw,i].fy -= ffy;

            worms[iw,ip3].fz += ffz;
            worms[iw,i].fz -= ffz;
        }
    }
    /*
    forall iw in 1..nworms {
        var i3:int,i4:int,x2:real,x3:real,x4:real,y2:real,y3:real,y4:real,z2:real,z3:real,z4:real,z23:real,z34:real,y23:real,y34:real,x23:real,x34:real,r23:real,r34:real,cosvalue:real;
	    var sinvalue:real,ff:real,dot:real,fac:real,f2x:real,f2y:real,f2z:real,f3x:real,f3y:real,f3z:real,f4x:real,f4y:real,f4z:real;
        for i2 in 1..(np-2) {
            i3 = i2 + 1;
            i4 = i2 + 2;
            //print*, i2,i3,i4
            x2 = worms[iw, i2].x;
            y2 = worms[iw, i2].y;
            z2 = worms[iw, i2].z;

            x3 = worms[iw, i3].x;
            y3 = worms[iw, i3].y;
            z3 = worms[iw, i3].z;

            x4 = worms[iw, i4].x;
            y4 = worms[iw, i4].y;
            z4 = worms[iw, i4].z;

            z23 = z3-z2;
            z34 = z4-z2;

            y23 = y3 - y2;
            y34 = y4 - y3;

            x23 = x3 - x2;
            x34 = x4 - x3;

            r23 = sqrt(x23*x23 + y23*y23 + z23*z23);
            r34 = sqrt(x34*x34 + y34*y34 + z34*z34);

            cosvalue = abs(x23*x34 + y23*y34 + z23*z34)/(r23*r34);
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

            dot = x23*x34 + y23*y34 + z23*z34;
            fac = dot/(r23*r23);

            f2x = ff*(x34 - fac*x23);
            f2y = ff*(y34 - fac*y23);
            f2z = ff*(z34 - fac*z23);

            fac = dot/(r34*r34);
            f4x = ff*(fac*x34 - x23);
            f4y = ff*(fac*y34 - y23);
            f4z = ff*(fac*z34 - z23);

            f3x = -f2x - f4x;
            f3y = -f2y - f4y;
            f3z = -f2z - f4z;

            worms[iw, i2].fx += f2x;
            worms[iw, i2].fy += f2y;
            worms[iw, i2].fz += f2z;

            worms[iw, i3].fx += f3x;
            worms[iw, i3].fy += f3y;
            worms[iw, i3].fz += f3z;

            worms[iw, i4].fx += f4x;
            worms[iw, i4].fy += f4y;
            worms[iw, i4].fz += f3z;
        }
    }
    */
}

proc calc_forces(istep:int,dt:real) {
    var binDomain: domain(1) = {1..(numBins*numBins)/2};
    forall binid in adaptive(binSpace) {
        if (bins[binid].ncount > 1) {
            // calculate the forces between atoms inside each bin
            for icount in 0..bins[binid].ncount-2 {
                var i = bins[binid].atoms[icount];
                var itype = bins[binid].types[icount];
                for jcount in (icount+1)..bins[binid].ncount-1 {
                    var j = bins[binid].atoms[jcount];
                    var jtype = bins[binid].types[jcount];
                    cell_forces(i,j,itype,jtype);
                }
            }
        }
    }
    // odd neighbors to the east (1)
    forall i in adaptive(binDomain) {
        var binid = binSpaceiodd[i];
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }
    // even neighbors to the east (1)
    forall i in adaptive(binDomain) {
        var binid = binSpaceieven[i];
        var binidnbor = bins[binid].neighbors[1];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // odd neighbors to the NE (2) (i+1,j+1)
    forall i in adaptive(binDomain) {
        var binid = binSpaceiodd[i];
        var binidnbor = bins[binid].neighbors[2];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // even neighbors to the NE (2)
    forall i in adaptive(binDomain) {
        var binid = binSpaceieven[i];
        var binidnbor = bins[binid].neighbors[2];
        if (binidnbor != -1) { // check if neighbor is valid
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) { // check if neighbor has any atoms at all
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // odd neighbors to the N (3)
    forall i in adaptive(binDomain) {
        var binid = binSpacejodd[i];
        var binidnbor = bins[binid].neighbors[3];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // even neighbors to the N (3)
    forall i in adaptive(binDomain) {
        var binid = binSpacejeven[i];
        var binidnbor = bins[binid].neighbors[3];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }

                }
            }
        }
    }

    // odd neighbors to the NW (4)
    forall i in adaptive(binDomain) {
        var binid = binSpaceiodd[i];
        var binidnbor = bins[binid].neighbors[4];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }

    // even neighbors to the NW (4)
    forall i in adaptive(binDomain) {
        var binid = binSpaceieven[i];
        var binidnbor = bins[binid].neighbors[4];
        if (binidnbor != -1) {
            if (bins[binid].ncount > 0) && (bins[binidnbor].ncount > 0) {
                for icount in 0..bins[binid].ncount-1 {
                    var i = bins[binid].atoms[icount];
                    var itype = bins[binid].types[icount];
                    for jcount in 0..bins[binidnbor].ncount-1 {
                        var j = bins[binidnbor].atoms[jcount];
                        var jtype = bins[binidnbor].types[jcount];
                        cell_forces(i,j,itype,jtype);
                    }
                }
            }
        }
    }
}

inline proc cell_forces(i:int,j:int,itype:int,jtype:int) {
    //worm-worm interaction 3D
    if (itype == 1) && (jtype == 1) {
        worm_cell_forces(i,j);
    }
    //worm-boundary interaction 2D
    if ((itype == 1) && (jtype == 3)) {
        var iw,ip:int;
        iw = 1 + ((i - 1)/np):int; // find which worm j is in
        ip = i - np*(iw - 1); // which particle in the worm is j?
        dogic_wall(iw,ip,j);
    } else if ((itype == 3) && (jtype == 1)) {
        var iw,ip:int;
        iw = 1 + ((j - 1)/np):int; // find which worm j is in
        ip = j - np*(iw - 1); // which particle in the worm is j?
        dogic_wall(iw,ip,i);
    }
    //solvent-worm interaction 2D
    if ((itype == 2) && (jtype == 1)) {
        var dx,dy,dz,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,fdrag,dv_dot_rhat:real;
        var iw,ip:int;
        iw = 1 + ((j - 1)/np):int; // find which worm j is in
        ip = j - np*(iw - 1); // which particle in the worm is j?
        dz = fluid_offset;
        dx = solvent[i].x - worms[iw,ip].x;
        dy = solvent[i].y - worms[iw,ip].y;
        r2 = (dx*dx + dy*dy + dz*dz);
        if (r2 <= r2cut) {
            //writeln("solvent worm ",sqrt(r2),"\t",sqrt(r2cut));
            //writeln("solvent worm ",r2,"\t",r2cut);
            ffor = -48.0*sw_epsilon*r2**(-7.0) + 24.0*sw_epsilon*r2**(-4.0); 
            ffx = ffor*dx;
            ffy = ffor*dy;

            // pairwise drag force - March 22
            if (fdrag > 0) {
                dvx = solvent[i].vx - worms[iw,ip].vx;
                dvy = solvent[i].vy - worms[iw,ip].vy;
                rhatx = dx/sqrt(r2);
                rhaty = dy/sqrt(r2);
                dv_dot_rhat = dvx*rhatx + dvy*rhaty;
                ffx += fdrag*dv_dot_rhat*dx;
                ffy += fdrag*dv_dot_rhat*dy;
            }
            solvent[i].fx -= ffx;
            solvent[i].fy -= ffy;
            worms[iw,ip].fx += ffx;
            worms[iw,ip].fy += ffy;
        }
    } else if ((itype == 1) && (jtype == 2)) {
        var dx,dy,dz,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,fdrag,dv_dot_rhat:real;
        var iw,ip:int;
        iw = 1 + ((i - 1)/np):int; // find which worm i is in
        ip = i - np*(iw - 1); // which particle in the worm is i?
        dz = fluid_offset;
        dx = solvent[j].x - worms[iw,ip].x;
        dy = solvent[j].y - worms[iw,ip].y;
        r2 = (dx*dx + dy*dy + dz*dz);
        if (r2 <= r2cut) {
            //writeln("worm solvent ",sqrt(r2),"\t",sqrt(r2cut));
            //writeln("worm solvent",r2,"\t",r2cut);
            ffor = -48.0*sw_epsilon*r2**(-7.0) + 24.0*sw_epsilon*r2**(-4.0);
            ffx = ffor*dx;
            ffy = ffor*dy;

            // pairwise drag force - March 22
            if (fdrag > 0) {
                dvx = solvent[i].vx - worms[iw,ip].vx;
                dvy = solvent[i].vy - worms[iw,ip].vy;
                rhatx = dx/sqrt(r2);
                rhaty = dy/sqrt(r2);
                dv_dot_rhat = dvx*rhatx + dvy*rhaty;
                ffx += fdrag*dv_dot_rhat*dx;
                ffy += fdrag*dv_dot_rhat*dy;
            }
            
            solvent[j].fx -= ffx;
            solvent[j].fy -= ffy;
            worms[iw,ip].fx += ffx;
            worms[iw,ip].fy += ffy;
        }
    }
    //solvent-solvent interaction 2D
    if (itype == 2) && (jtype == 2) {
        lj_thermo(i,j,r2cutsmall);
    }
    //solvent-boundary interaction 2D
    if ((itype == 2) && (jtype == 3)) {
        lj(j,i,sigma*r2cutsmall);
    } else if ((itype == 3) && (jtype == 2)) {
        lj(i,j,sigma*r2cutsmall);
    }
}

inline proc worm_cell_forces(i:int,j:int) {
   var dddx:real,dddy:real,dddz:real,r2:real,riijj:real,ffor:real,ffx:real,ffy:real,ffz:real,dxi:real,dxj:real,
       ri:real,rj:real,r:real,dx:real,dy:real,dz:real,dyi:real,dyj:real,dzi:real,dzj:real,r2shift:real;
   var dvx:real,dvy:real,dvz:real,dv_dot_rhat:real,rhatx:real,rhaty:real,rhatz:real,omega:real,fdissx:real,fdissy:real,fdissz:real,gauss:real,frand:real;
   var iworm:int,jworm:int,ip:int,jp:int,ip1:int,jp1:int,inogo:int;
   iworm = 1 + ((i - 1)/np):int;
   ip = i - np*(iworm - 1);
   jworm = 1 + ((j - 1)/np):int;
   jp = j - np*(jworm - 1);

   inogo = 0;
   if ((iworm == jworm) && (abs(ip-jp) <= 2)) {
      // on the same worm and close means no interaction calculated here
      inogo = 1;
   }
   if (inogo == 0) {
      dddx = worms[jworm, jp].x - worms[iworm, ip].x;
      dddy = worms[jworm, jp].y - worms[iworm, ip].y;
      dddz = worms[jworm, jp].z - worms[iworm, ip].z;
      r2 = dddx**2 + dddy**2 + dddz**2;
      riijj = sqrt(r2);
      //add attractive force fdep between all pairs
      if (r2 <= r2cutsmall) {
            ffor = -48.0*ww_epsilon*r2**(-7.0) + 24.0*ww_epsilon*r2**(-4.0) + fdep/riijj; //TODO: shoudl fdep = -?
            //ffor = -48.0*(r2-r2shift)**(-7.0)  + 24.0*(r2-r2shift)**(-4.0) + fdep/riijj;
            ffx = ffor*dddx;
            ffy = ffor*dddy;
            ffz = ffor*dddz;
            worms[iworm,ip].fx += ffx;
            worms[jworm,jp].fx -= ffx;
            worms[iworm,ip].fy += ffy;
            worms[jworm,jp].fy -= ffy;
            worms[iworm,ip].fz += ffz;
            worms[jworm,jp].fz -= ffz;
            

            // DPD thermostat https://docs.lammps.org/pair_dpd.html
            //adding dissipative force
            if (thermow) {
            dvx = worms[jworm,jp].vx - worms[iworm,ip].vx;
            dvy = worms[jworm,jp].vy - worms[iworm,ip].vy;
            dvz = worms[jworm,jp].vz - worms[iworm,ip].vz;

            rhatx = dddx/riijj;
            rhaty = dddy/riijj;
            rhatz = dddz/riijj;

            omega = (1.0-riijj/rcutsmall);
            dv_dot_rhat = dvx*rhatx + dvy*rhaty + dvz*rhatz;
            fdissx = -1.0*gamma*(omega*omega)*dv_dot_rhat*rhatx; //gamma = 1/damp (proportional to friction force)
            fdissy = -1.0*gamma*(omega*omega)*dv_dot_rhat*rhaty;
            fdissz = -1.0*gamma*(omega*omega)*dv_dot_rhat*rhatz;

            worms[iworm,ip].fx -= fdissx;
            worms[iworm,ip].fy -= fdissy;
            worms[iworm,ip].fz -= fdissz;

            worms[jworm,jp].fx += fdissx;
            worms[jworm,jp].fy += fdissy;
            worms[jworm,jp].fz += fdissz;

            // adding random forces
            gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
            frand = dpd_ratio*(inv_sqrt_dt)*omega*gauss*sqrt_gamma_term;

            worms[iworm,ip].fx += frand*rhatx;
            worms[iworm,ip].fy += frand*rhaty;
            worms[iworm,ip].fz += frand*rhatz;

            worms[jworm,jp].fx -= frand*rhatx;
            worms[jworm,jp].fy -= frand*rhaty;
            worms[jworm,jp].fz -= frand*rhatz;
            }
            //vxave,vyave
            worms[iworm,ip].vxave += worms[jworm,jp].vx;
            worms[iworm,ip].vyave += worms[jworm,jp].vy;
            worms[iworm,ip].vzave += worms[jworm,jp].vz;

            worms[jworm,jp].vxave += worms[iworm,ip].vx;
            worms[jworm,jp].vyave += worms[iworm,ip].vy;
            worms[jworm,jp].vzave += worms[iworm,ip].vz;

            nnab[iworm,ip] += 1;
            nnab[jworm,jp] += 1;

            //add 'dogic drive' to interacting pairs
            //first calculate unit vectors along each worm
            ip1 = ip + 1;
            if (ip1 <= np) {
               dxi = worms[iworm,ip1].x - worms[iworm, ip].x;
               dyi = worms[iworm,ip1].y - worms[iworm, ip].y;
               dzi = worms[iworm,ip1].z - worms[iworm, ip].z;
            } else {
               dxi = worms[iworm, ip].x - worms[iworm, ip - 1].x;
               dyi = worms[iworm, ip].y - worms[iworm, ip - 1].y;
               dzi = worms[iworm, ip].z - worms[iworm, ip - 1].z;
            }

            jp1 = jp + 1;
            if (jp1 <= np) {
               dxj = worms[jworm, jp1].x - worms[jworm, jp].x;
               dyj = worms[jworm, jp1].y - worms[jworm, jp].y;
               dzj = worms[jworm, jp1].z - worms[jworm, jp].z;
            } else {
               dxj = worms[jworm, jp].x - worms[jworm, jp - 1].x;
               dyj = worms[jworm, jp].y - worms[jworm, jp - 1].y;
               dzj = worms[jworm, jp].z - worms[jworm, jp - 1].z;
            }


            //if the two vectors have any component pointing in opposite directions
            if (dxi*dxj + dyi*dyj + dzi*dzj <= -0.5) {
            //if (dxi*dxj + dyi*dyj + dzi*dzj<= 0.0) {
               //normalize those vectors to make them unit vectors
               ri = sqrt(dxi*dxi + dyi*dyi *dzi*dzi);
               dxi = dxi/ri;
               dyi = dyi/ri;
               dzi = dzi/ri;

               rj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);
               dxj = dxj/rj;
               dyj = dyj/rj;
               dzj = dzj/rj;
               //now they are both unit vectors. Find the direction for the force...

               dx = (dxi - dxj)/2.0;
               dy = (dyi - dyj)/2.0;
               dz = (dzi - dzj)/2.0;

               //normalize

               r = sqrt(dx*dx + dy*dy + dz*dz);
               dx = dx/r;
               dy = dy/r;
               dz = dz/r;

               //add an extra attractive component where kinesin drive is present

               ffx = fdogic*(dx) + dogic_fdep*dddx/riijj;
               ffy = fdogic*(dy) + dogic_fdep*dddy/riijj;
               ffz = fdogic*(dz) + dogic_fdep*dddz/riijj;

               worms[iworm,ip].fx += ffx;
               worms[jworm,jp].fx -= ffx;
               worms[iworm,ip].fy += ffy;
               worms[jworm,jp].fy -= ffy;
               worms[iworm,ip].fz += ffz;
               worms[jworm,jp].fz -= ffz;
            }
      }
   }
}

inline proc dogic_wall(iw:int,ip:int,ib:int){
    var dx:real, dy:real, r:real, r2:real, th:real,
        xwall:real, ywall:real, rr2:real, ffor:real,
        dxi:real, dyi:real, ri:real, dxj:real, dyj:real, ffx:real, ffy:real;
    var ip1,ib1:int;
    //calculate distance to the wall
    dx = worms[iw,ip].x-bound[ib].x;
    dy = worms[iw,ip].y-bound[ib].y;
    r2 = (dx*dx + dy*dy);
    //if close enough to the wall, calculate wall forces
    //use the short cut-off
    if (r2 <= r2cutsmall) {
        r = sqrt(r2);
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdepwall/r;
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        //ffor = (1/r)*6.0 + fdepwall/r; //TODO raise this to a higher power to get the worms closer to the wall? try ^6 or ^8
        ffor = (1/r)**6.0 + fdepwall/r;
        worms[iw,ip].fx += ffor*dx;
        worms[iw,ip].fy += ffor*dy;
        if (walldrive) {
            //first calculate unit vector along the worm
            ip1 = ip + 1;
            if (ip1 <= np) {
                dxi = worms[iw, ip1].x - worms[iw, ip].x;
                dyi = worms[iw, ip1].y - worms[iw, ip].y;
            } else {
                dxi = worms[iw, ip].x - worms[iw, ip - 1].x;
                dyi = worms[iw, ip].y - worms[iw, ip - 1].y;
            }
            //make it a unit vector
            ri = sqrt(dxi*dxi + dyi*dyi);
            dxi = dxi/ri;
            dyi = dyi/ri;
            //calculate the unit vector along the wall
            ib1 = ib + 1;
            if (ib1 <= numPoints) {
                dxj = bound[ib1].x - bound[ib].x;
                dyj = bound[ib1].y - bound[ib].y;
            } else {
                dxj = bound[ib].x - bound[ib-1].x;
                dyj = bound[ib].y - bound[ib-1].y;
            }
            ri = sqrt(dxj*dxj + dyj*dyj);
            dxj = dxj/ri;
            dyj = dyj/ri;

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
                worms[iw,ip].fx += ffx;
                worms[iw,ip].fy += ffy;
            }
        }
    }
}

proc update_vel() {
    forall iw in 1..nworms {
        forall i in 1..np {
            worms[iw, i].vx += dto2*(worms[iw,i].fx + worms[iw,i].fxold);
            worms[iw, i].vy += dto2*(worms[iw,i].fy + worms[iw,i].fyold);
            //worms[iw, i].vz += dto2*(worms[iw,i].fz + worms[iw,i].fzold);

            worms[iw, i].vxave = worms[iw, i].vxave/nnab[iw, i]:real;
            worms[iw, i].vyave = worms[iw, i].vyave/nnab[iw, i]:real;
            worms[iw, i].vzave = worms[iw, i].vzave/nnab[iw, i]:real;

            KEworm[iw*i] = 0.5*worms[iw,i].m*(worms[iw,i].vx * worms[iw,i].vx + worms[iw,i].vy * worms[iw,i].vy);
            KEworm_local[iw*i] = 0.5*worms[iw,i].m*((worms[iw,i].vx-worms[iw,i].vxave) * (worms[iw,i].vx-worms[iw,i].vxave) + (worms[iw,i].vy-worms[iw,i].vyave) * (worms[iw,i].vy-worms[iw,i].vyave));
            AMworm[iw*i] = (worms[iw,i].x-hxo2)*worms[iw,i].vy - (worms[iw,i].y-hyo2)*worms[iw,i].vx; // angular momentum calculation

            nnab[iw,i] = 1;
            worms[iw,i].vxave = 0.0;
            worms[iw,i].vyave = 0.0;
            worms[iw,i].vzave = 0.0;
        }
    }
    fluid_vel(dt);
}

//FLUID FUNCTIONS

inline proc lj_thermo(i:int,j:int,r2cut_local:real) {
    var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,r,rhatx,rhaty,omega,fdissx,fdissy,gauss,frand,dv_dot_rhat :real;
    dx = solvent[j].x - solvent[i].x;
    dy = solvent[j].y - solvent[i].y;
    r2 = (dx*dx + dy*dy);
    //if (debug) {writeln("icount: ",icount,"\t",i,"\tjcount: ",jcount,"\t",j,"\t",r2,"\t",r2cut);}
    if (r2 <= r2cut_local) {
        // LJ force
        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        ffx = ffor*dx;
        ffy = ffor*dy;
        solvent[i].fx += ffx;
        solvent[i].fy += ffy;
        solvent[j].fx -= ffx;
        solvent[j].fy -= ffy;
        if (thermo) {
            // DPD thermostat
            //adding dissipative force
            dvx = solvent[j].vx - solvent[i].vx;
            dvy = solvent[j].vy - solvent[i].vy;
            r = sqrt(r2);
            rhatx = dx/r;
            rhaty = dy/r;
            omega = (1.0-r/sqrt(r2cut_local));
            dv_dot_rhat = dvx*rhatx + dvy*rhaty;
            fdissx = -1.0*gamma*(omega*omega)*dv_dot_rhat*rhatx; //gamma = 1/damp (proportional to friction force)
            fdissy = -1.0*gamma*(omega*omega)*dv_dot_rhat*rhaty;
            solvent[i].fx -= fdissx;
            solvent[i].fy -= fdissy;
            solvent[j].fx += fdissx;
            solvent[j].fy += fdissy;
            // adding random forces
            gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
            //gauss = randStream.next();
            frand = dpd_ratio*(inv_sqrt_dt)*omega*gauss*sqrt_gamma_term;
            solvent[i].fx += frand*rhatx;
            solvent[i].fy += frand*rhaty;
            solvent[j].fx -= frand*rhatx;
            solvent[j].fy -= frand*rhaty;
        }
    }
}

inline proc lj(i:int,j:int,r2cut_local:real) {
    var dx,dy,r2,ffor,ffx,ffy,sigma12,sigma6:real;
    //calculate distance to the wall
    dx = solvent[j].x - bound[i].x;
    dy = solvent[j].y - bound[i].y;
    r2 = (dx*dx + dy*dy);
    //if close enough to the wall, calculate wall forces
    //use the short cut-off
    if (r2 <= r2cut_local) {
        //ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
        //ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        //sigma12 = sigma**12;
        //sigma6 = sigma**6;
        //ffor = 48.0*sigma12*r2**(-7.0) -24*sigma6*r2**(-4.0);
        ffor = (1/r2)**2.0;
        solvent[j].fx += ffor*dx;
        solvent[j].fy += ffor*dy;
    }
}

inline proc fluid_pos(dt_fluid:real) {
    forall i in 1..numSol {
        solvent[i].x += solvent[i].vx*dt_fluid + (solvent[i].fx/2)*(dt_fluid**2);
        solvent[i].y += solvent[i].vy*dt_fluid + (solvent[i].fy/2)*(dt_fluid**2);
        // updating velocity and position
        // solvx[i] += solfx[i] * dt;
        // solvy[i] += solfy[i] * dt;
        // solx[i] += solvx[i] * dt;
        // soly[i] += solvy[i] * dt;

        // doing a velocity half-step here
        solvent[i].fx = solvent[i].fx/solvent[i].m;
        solvent[i].fy = solvent[i].fy/solvent[i].m;
        solvent[i].vx += 0.5*dt_fluid*solvent[i].fx;
        solvent[i].vy += 0.5*dt_fluid*solvent[i].fy;
        solvent[i].vz = 0.0; // just in case
        solvent[i].fx = 0.0;
        solvent[i].fy = 0.0;
        solvent[i].fz = 0.0; // just in case
    }
}

proc fluid_force_old() {
         // fluid-fluid
       var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,r,frand,gauss,fdissx,fdissy,omega:real;
         for i in 1..numSol {
             // calculating the force
             for j in (i+1)..numSol {
                 lj_thermo(i,j,r2cut);
             }
         }
         //fluid-boundary
         // var r2cutsol = r2cut
         var r2cutsol = sigma*r2cutsmall;
         forall i in 1..numSol {
             // calculate the force on the boundaries.
             for ib in 1..numPoints  {
                 lj(ib,i,r2cutsol);

             }
            //lj(1,i,r2cutsol);
         }

    //     // fluid-worms
    //     var dz = fluid_offset;
    //     for i in 1..numSol{
    //         var dx,dy,r2,ffor,ffx,ffy:real;
    //         for iw in 1..nworms {
    //             for ip in 1..np {
    //                 dx = solvent[i].x - worms[iw,i].x;
    //                 dy = solvent[i].y - worms[iw,i].y;
    //                 r2 = (dx*dx + dy*dy + dz*dz);
    //                 if (r2 <= r2cutsmall) {
    //                     ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
    //                     ffx = ffor*dx;
    //                     ffy = ffor*dy;
    //                     solvent[i].fx += ffx;
    //                     solvent[i].fy += ffy;
    //                     worms[iw,ip].fx -= ffx;
    //                     worms[iw,ip].fy -= ffy;
    //                 }
    //             }
    //         }
    //     }
}

inline proc fluid_vel(dt_fluid:real) {
    //update velocities
    forall i in 1..numSol {
        //solvent[i].vx += 0.5*dt_fluid*(solvent[i].fxold + solvent[i].fx);
        //solvent[i].vy += 0.5*dt_fluid*(solvent[i].fyold + solvent[i].fy);
        //solvent[i].fxold = solvent[i].fx;
        //solvent[i].fyold = solvent[i].fy;
        solvent[i].fx = solvent[i].fx/solvent[i].m;
        solvent[i].fy = solvent[i].fy/solvent[i].m;
        solvent[i].vx += 0.5*dt_fluid*solvent[i].fx;
        solvent[i].vy += 0.5*dt_fluid*solvent[i].fy;
        // calculating kinetic energy here too
        KEsol[i] = 0.5*solvent[i].m*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
        AMsol[i] = (solvent[i].x-hxo2)*solvent[i].vy - (solvent[i].y-hyo2)*solvent[i].vx; // angular momentum calculation
    }
}

proc fluid_step(istep:int,dt_fluid:real) {
    // update positions
    fluid_pos(dt_fluid);
    // update_fluid_cells(istep);
    fluid_force_old();

    fluid_vel(dt_fluid);
}
}