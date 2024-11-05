import warp as wp

@wp.func
def LJ_force(Xij: wp.vec3, r2:float,eps:float,sigma:float,fdep:float):
    return (48.0*eps*r2**(-7.0) - 24.0*eps*r2**(-4.0) + fdep)*Xij

@wp.func
def dist_sq(X: wp.vec3):
    return wp.dot(X,X)

@wp.kernel
def update_positions(grid: wp.uint64,
                     P: wp.array(dtype=wp.vec3),
                     V: wp.array(dtype=wp.vec3),
                     F: wp.array(dtype=wp.vec3),
                     Fold: wp.array(dtype=wp.vec3),
                     dt: float,Lx: float, Ly: float, Lz: float):
    tid = wp.tid()
    i = wp.hash_grid_point_id(grid,tid)
    P[i] += V[i]*dt + F[i]*dt*dt*0.5
    #P[i] = apply_periodic_bc(P[i],Lx,Ly,Lz)
    Fold[i] = F[i]
    F[i] = wp.vec3(0.0,0.0,0.0)
    
@wp.kernel
def update_velocities(grid: wp.uint64,
                      V: wp.array(dtype=wp.vec3),
                      F: wp.array(dtype=wp.vec3),
                      Fold: wp.array(dtype=wp.vec3),
                      dt: float):
    tid = wp.tid()
    i = wp.hash_grid_point_id(grid,tid)
    V[i] += dt * 0.5 * (F[i] + Fold[i])
    
@wp.kernel
def apply_forces(grid: wp.uint64,
                 F: wp.array(dtype=wp.vec3),
                 X: wp.array(dtype=wp.vec3),
                 V: wp.array(dtype=wp.vec3),
                 rcut: float, nworms:int, np:int, fdep:float,kbt:float,dt:float):
    tid = wp.tid()
    eps = 1.0
    sigma = 1.0
    i = wp.hash_grid_point_id(grid,tid)
    Xi = X[i]
    Vi = V[i]
    Fi = wp.vec3(0.0,0.0,0.0)
    rcut2 = rcut*rcut
    iworm = int(i/np)
    ip = i - np*iworm

    neighbors = wp.hash_grid_query(grid,Xi,rcut)
    for j in neighbors:
        if j != i:
            Xj = X[j]
            jworm = int(i/np)
            jp = j - np*jworm
            if iworm != jworm:
                # calculate LJ force
                Xij = Xi - Xj
                r2 = dist_sq(Xij)
                if r2 <= rcut2:
                    Fi += LJ_force(Xij,wp.sqrt(r2),eps,sigma,fdep)
                    #Fi += langevin_thermo()
                #if r2 < rcl2:
                    # add cutoff
            elif (iworm == jworm and abs(ip-jp) > 5):
                # calculate LJ force
                Xij = Xi - Xj
                r2 = dist_sq(Xij)
                if r2 <= rcut2:
                    Fi += LJ_force(Xij,wp.sqrt(r2),eps,sigma,fdep)
                    #Fi += langevin_thermo()
                #if r2 < rcl2:
                    # add cutoff
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
 