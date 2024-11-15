import warp as wp

@wp.func
def LJ_force(Xij: wp.vec3, r2:float,eps:float,sigma:float,fdep:float):
    return (48.0*eps*r2**(-7.0) - 24.0*eps*r2**(-4.0) + fdep)*Xij

@wp.func
def dist_sq(X: wp.vec3):
    return wp.dot(X,X)

@wp.func
def dist_sq(X: wp.vec2):
    return wp.dot(X,X)

@wp.kernel
def update_positions(grid: wp.uint64,
                     P: wp.array(dtype=wp.vec3),
                     V: wp.array(dtype=wp.vec3),
                     F: wp.array(dtype=wp.vec3),
                     Fold: wp.array(dtype=wp.vec3),
                     dt: float, L: float, nactive:wp.int32,nmobile:wp.int32):
    tid = wp.tid()
    i = wp.hash_grid_point_id(grid,tid)
    if i <= nactive: #active particle update
        P[i] += V[i]*dt + F[i]*dt*dt*0.5
        if P[i][2] > L:
            P[i][2] = L
        elif P[i][2] < 0.0:
            P[i][2] = 0.0
    elif (i > nactive) and (i <= nmobile): # fluid particle update
        #solvent[i].x += solvent[i].vx*dt_fluid + (solvent[i].fx/2)*(dt_fluid**2);
        #solvent[i].y += solvent[i].vy*dt_fluid + (solvent[i].fy/2)*(dt_fluid**2);
        P[i] += V[i]*dt + F[i]*dt*dt*0.5
        P[i][2] = 0.0
        # solvent[i].vx += 0.5*dt_fluid*solvent[i].fx; //velocity half-step here
        # solvent[i].vy += 0.5*dt_fluid*solvent[i].fy;
        V[i] += 0.5*dt*F[i]
    Fold[i] = F[i]
    F[i] = wp.vec3(0.0,0.0,0.0)
    
@wp.kernel
def update_velocities(grid: wp.uint64,
                      V: wp.array(dtype=wp.vec3),
                      F: wp.array(dtype=wp.vec3),
                      Fold: wp.array(dtype=wp.vec3),
                      dt: float,
                      nactive:wp.int32,nmobile:wp.int32):
    tid = wp.tid()
    i = wp.hash_grid_point_id(grid,tid)
    if i <= nactive: #active particle update
        V[i] += dt * 0.5 * (F[i] + Fold[i])
    elif (i > nactive) and (i <= nmobile): #fluid particle update
        V[i] += 0.5*dt*F[i]

    
@wp.kernel
def apply_forces(grid: wp.uint64,
                 F: wp.array(dtype=wp.vec3),
                 X: wp.array(dtype=wp.vec3),
                 V: wp.array(dtype=wp.vec3),
                 T: wp.array(dtype=wp.int32),
                 rcut: float, rcutsmall:float, nworms:int, np:int, fdep:float,kbt:float,dt:float,fdogic:float,dogic_fdep:float,fdepwall:float,fluid_offset:float,fdrag:float):
    tid = wp.tid()
    eps = 1.0
    sigma = 1.0
    ww_epsilon = 0.5
    sw_epsilon = 2.0
    i = wp.hash_grid_point_id(grid,tid)
    Xi = X[i]
    Vi = V[i]
    Ti = T[i]
    Fi = wp.vec3(0.0,0.0,0.0)
    rng = wp.rand_init(9383421,tid)
    r2cut = rcut*rcut
    neighbors = wp.hash_grid_query(grid,Xi,rcut)
    for j in neighbors:
        Xj = X[j]
        Vj = V[j]
        Tj = T[j]
        # worm-worm interaction
        if (Ti == 1) and (Tj == 1):
            # both are worms
            iworm = int(i/np)
            ip = i - np*iworm
            jworm = int(j/np)
            jp = j - np*jworm

            r2cutsmall = rcutsmall*rcutsmall
            if ((iworm == jworm) and (wp.abs(ip-jp) <= 2)):
                # particles are too close
                pass
            else:
                Xij = Xi - Xj
                r2 = dist_sq(Xij)
                if r2 <= r2cutsmall:
                    riijj = wp.sqrt(r2)
                    ffor = -48.0*ww_epsilon*r2**(-7.0) + 24.0*ww_epsilon*r2**(-4.0) + fdep/riijj
                    Fi += ffor*Xij
                    """
                    # DPD thermostat https://docs.lammps.org/pair_dpd.html
                    #if thermow:
                    #adding dissipative force
                    dVij = Vi - dVij
                    rhat = Xij / riijj
                    omega = (1.0-riijj/rcutsmall)
                    dv_dot_rhat = wp.dot(dVij,rhat)
                    fdiss = -1.0*gamma*(omega*omega)*dv_dot_rhat*rhat # gamma = 1/damp (proportional to friction force)
                    Fi -= fdiss
                    # adding random force
                    gauss = wp.randn(r)
                    frand = ((inv_sqrt_dt)*omega*gauss*sqrt_gamma_term)*rhat
                    Fi += frand
                    """
                    #add 'dogic drive' to interacting pairs
                    #first calculate unit vectors along each worm
                    ipp1 = ip + 1
                    if (ipp1 <= np):
                        ip1 = (iworm)*np+ipp1
                        dXi = X[ip1] - Xi
                    else:
                        im1 = (iworm)*np+(ip-1)
                        dXi = Xi - X[im1]
                    jpp1 = jp + 1
                    if jpp1 <= np:
                        jp1 = (jworm)*np+jpp1
                        dXj = X[jp1] - Xj
                    else:
                        jm1 = (jworm)*np+(jp-1)
                        dXj = Xi - X[jm1]
                    # if the two vectors have any component pointing in opposite directions
                    if (wp.dot(dXi,dXj) <= -0.5):
                        # normalize it
                        ri = wp.sqrt(wp.dot(dXi,dXi))
                        dXi = dXi/ri
                        rj = wp.sqrt(wp.dot(dXj,dXj))
                        dXj = dXj/rj
                        #now they are both unit vectors. Find the direction for the force...
                        dF = (dXi-dXj)/2.0
                        # and normalize
                        r = wp.sqrt(wp.dot(dF,dF))
                        dF = dF/r
                        #add an extra attractive component where kinesin drive is present
                        ff = fdogic*dF + dogic_fdep*Xij/riijj
                        Fi += ff
        # worm-boundary interaction 2D
        elif (Ti == 1) and (Tj == 3):
            #calculate distance to the wall
            dX = wp.vec3(Xi[0] - Xj[0],Xi[1] - Xj[1],0.0)
            r2 = dist_sq(dX)
            if r2 <= r2cutsmall:
                r = sqrt(r2)
                #ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdepwall/r;
                #ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                ffor = (1.0/r)*6.0 + fdepwall/r # raise this to a higher power to get the worms closer to the wall? try ^6 or ^8
                #ffor = (1/r)**6.0 - fdepwall/r;
                ff = ffor*dX
                #ff3d = wp.vec3(ff[0],ff[1],0.0)
                Fi += ff
        # worm-fluid interaction
        elif (Ti == 2) and (Tj == 1):
            dX = wp.vec3(Xi[0] - Xj[0],Xi[1] - Xj[1],fluid_offset)
            r2 = dist_sq(dX)
            if r2 <= r2cut:
                ffor = -48.0*sw_epsilon*r2**(-7.0) + 24.0*sw_epsilon*r2**(-4.0)
                ffx = ffor*dX[0]
                ffy = ffor*dX[1]
                if fdrag > 0: #pairwise drag force
                    dV = wp.vec2(Vi[0]-Vj[0],Vi[1]-Vj[1])
                    rhat = dX/r2
                    rhat2d = wp.vec2(rhat[0],rhat[1])
                    dv_dot_rhat = wp.dot(dV,rhat2d)
                    ffx += fdrag*dv_dot_rhat*dX[0]
                    ffy += fdrag*dv_dot_rhat*dX[1]
                Fi += wp.vec3(ffx,ffy,0.0)
        elif (Ti == 1) and (Tj == 2):
            dX = wp.vec3(Xj[0] - Xi[0],Xj[1] - Xi[1],fluid_offset)
            r2 = dist_sq(dX)
            if r2 <= r2cut:
                ffor = -48.0*sw_epsilon*r2**(-7.0) + 24.0*sw_epsilon*r2**(-4.0)
                ffx = ffor*dX[0]
                ffy = ffor*dX[1]
                if fdrag > 0: #pairwise drag force
                    dV = wp.vec2(Vj[0]-Vi[0],Vj[1]-Vi[1])
                    rhat = dX/r2
                    rhat2d = wp.vec2(rhat[0],rhat[1])
                    dv_dot_rhat = wp.dot(dV,rhat2d)
                    ffx += fdrag*dv_dot_rhat*dX[0]
                    ffy += fdrag*dv_dot_rhat*dX[1]
                Fi -= wp.vec3(ffx,ffy,0.0)
        # fluid-fluid interaction
        elif (Ti == 2) and (Tj == 2):
            dX = wp.vec3(Xj[0] - Xi[0],Xj[1]-Xi[1],0.0)
            r2 = dist_sq(dX)
            if r2 <= r2cut:
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0)
                #ff  = ffor*dX
                Fi += ffor*dX
        # fluid-boundary
        elif (Ti == 2) and (Tj == 3):
            dX = wp.vec3(Xi[0]-Xj[0],Xi[1]-Xj[1],0.0)
            r2 = dist_sq(dX)
            if r2 <= r2cut:
                ffor = (1.0/r2)**2.0
                #ff = ffor*dX
                Fi += ffor*dX
    if (Ti == 3):
        Fi = wp.vec3(0.0,0.0,0.0)    
    F[i] = Fi   

@wp.kernel
def bond_bending_forces(grid: wp.uint64,
                        F: wp.array(dtype=wp.vec3),
                        X: wp.array(dtype=wp.vec3),
                        bond_offset: int,np: int, kspring: float, length0: float):
    iw,ip = wp.tid()
    #wp.printf("%d %d %d\n",iw,ip,ip+bond_offset)
    id = (iw+1)*np - (ip+1)
    idp1 = (iw+1)*np - (ip+bond_offset+1)
    Xij = X[idp1] - X[id]
    r = wp.sqrt(dist_sq(Xij))
    ff = -kspring*(r - length0)/r
    F[idp1] += wp.vec3(ff*X[idp1][0],ff*X[idp1][1],ff*X[idp1][2])
    F[id] -= wp.vec3(ff*X[id][0],ff*X[id][1],ff*X[id][2])
 