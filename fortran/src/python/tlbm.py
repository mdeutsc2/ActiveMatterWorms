# taichi-based lbm fluid solver for just the fluid flow example
import os,sys
if sys.platform == "darwin":
    os.environ['KMP_DUPLICATE_LIB_OK']='True'
import timeit
import taichi as ti
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt


ti.init(arch=ti.cpu)
show_plot = False
save_vtk = False
save_interval = 1000

# PARAMS
nx = 801#domain size
ny = 201
niu = 0.01 #viscosity
steps = 20
tau = 3.0*niu + 0.5
inv_tau = 1.0/tau
rho = ti.field(dtype=ti.f32, shape=(nx,ny))
vel = ti.Vector.field(2,dtype=ti.f32, shape=(nx,ny))
mask = ti.field(dtype=ti.f32, shape=(nx,ny))
f_old = ti.Vector.field(9,dtype=ti.f32, shape=(nx,ny))
f_new = ti.Vector.field(9,dtype=ti.f32, shape=(nx,ny))
w = ti.field(dtype=ti.f32, shape=9)
e = ti.field(dtype=ti.f32, shape=(9,2))

w.from_numpy(np.array([ 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0], dtype=np.float32))
e.from_numpy(np.array([[0, 0],
                        [1, 0], 
                        [0, 1], 
                        [-1, 0], 
                        [0, -1], 
                        [1, 1],
                        [-1, 1], 
                        [-1, -1], 
                        [1, -1]], dtype=np.int32))

# [left,top,right,bottom] =0, Dirichlet, =1, Neumann
bc_type = ti.field(dtype=ti.f32,shape=4)
bc_type.from_numpy(np.array([0,0,1,0], dtype=np.int32))
# if bc_type = 0 we need to specify a velocity in bc_value
bc_value = ti.field(dtype=ti.f32, shape=(4,2))
bc_value.from_numpy(np.array([[0.1,0.0],
                                [0.0,0.0],
                                [0.0,0.0],
                                [0.0,0.0]], dtype=np.float32))
cy = 0 # whether to place cylindrical object
cy_para = ti.field(dtype=ti.f32, shape=3)
cy_para.from_numpy(np.array([160.0,100.0,20.0], dtype=np.float32))

@ti.func
def f_eq(i,j,k,e,vel,w,rho):
    eu = ti.cast(e[k,0],ti.f32)*vel[i,j][0] + ti.cast(e[k,1],ti.f32)*vel[i,j][1]
    uv = vel[i,j][0]**2.0 + vel[i,j][1]**2.0
    return w[k]*rho[i,j]*(1.0 + 3.0*eu + 4.5*eu**2 - 1.5*uv)

@ti.kernel
def lbm_init():
    for i,j in rho:
        vel[i,j][0] = 0.0
        vel[i,j][1] = 0.0
        rho[i,j] = 1.0
        mask[i,j] = 0.0
        for k in ti.static(range(9)):
            f_new[i,j][k] = f_eq(i,j,k,e,vel,w,rho)
            f_old[i,j][k] = f_new[i,j][k]
        if (cy==1):
            if ((ti.cast(i,ti.f32) - cy_para[0])**2.0 + (ti.cast(j,ti.f32) - cy_para[1])**2.0 <= cy_para[2]**2.0):
                mask[i,j] = 1.0

@ti.kernel
def collide_and_stream():
    for i,j in ti.ndrange((1,nx-1),(1,ny-1)):
        for k in ti.static(range(9)):
            ip = i - ti.cast(e[k,0],ti.int32)
            jp = j - ti.cast(e[k,1],ti.int32)
            f_new[i,j][k] = (1.0-inv_tau)*f_old[ip,jp][k] + f_eq(ip,jp,k,e,vel,w,rho)*inv_tau

@ti.kernel
def update_macro_vars():
    for i,j in ti.ndrange((1,nx-1),(1,ny-1)):
        rho[i,j] = 0.0
        vel[i,j][0] = 0.0
        vel[i,j][1] = 0.0
        for k in ti.static(range(9)):
            f_old[i,j][k] = f_new[i,j][k]
            rho[i,j] += f_new[i,j][k]
            vel[i,j][0] += (ti.cast(e[k,0],ti.f32)*f_new[i,j][k])
            vel[i,j][1] += (ti.cast(e[k,1],ti.f32)*f_new[i,j][k])
        
        vel[i,j][0] /= rho[i,j]
        vel[i,j][1] /= rho[i,j]

@ti.kernel
def apply_bc():
    # left and right bc
    #print("dr,ibc,jbc,inb,jnb,bc_type[dr]")
    for j in ti.ndrange(1,ny-1):
        # left: dr = 0, ibc=0, jbc = j, inb = 1, jnb = j
        apply_bc_core(1,0,0,j,1,j,vel,bc_type,bc_value,rho,e,w)

        #right: dr = 2, ibc = nx-1, jbc = j, inb = nx-2, jnb = j
        apply_bc_core(1,2,nx-1,j,nx-2,j,vel,bc_type,bc_value,rho,e,w)
    
    # top and bottom bc
    for i in ti.ndrange(nx):
        #print(i)
        #top: dr = 1, ibc=i, jbc=ny-1, inb = i, jnb = ny-2
        apply_bc_core(1,1,i,ny-1,i,ny-2,vel,bc_type,bc_value,rho,e,w)

        #bottom: dr = 3, ibc = 1, jbc = 0, inb = i, jnb = 1
        #print(3,1,0,i,1,bc_type[3])
        apply_bc_core(1,3,i,0,i,1,vel,bc_type,bc_value,rho,e,w)
    # cylindrical obstacle
    for i,j in ti.ndrange(nx,ny):
        if (cy == 1) and mask[i,j] == 1:
            vel[i,j][0] = 0.0 # velocity is zero at solid boundary
            vel[i,j][1] = 0.0
            inb = 0
            jnb = 0
            if (ti.cast(i,ti.f32) >= cy_para[0]):
                inb = i+1
            else:
                inb = i-1
            if (ti.cast(j,ti.f32) >= cy_para[1]):
                jnb = j+1
            else:
                jnb = j-1
            apply_bc_core(0,0,i,j,inb,jnb,vel,bc_type,bc_value,rho,e,w)

@ti.func
def apply_bc_core(outer,dr,ibc,jbc,inb,jnb,vel,bc_type,bc_value,rho,e,w):
    if (outer == 1): #handle outer boundary
        if (bc_type[dr] == 0):
            vel[ibc,jbc][0] = bc_value[dr,0]
            vel[ibc,jbc][1] = bc_value[dr,1]
        elif (bc_type[dr] == 1):
            vel[ibc,jbc][0] = vel[inb,jnb][0]
            vel[ibc,jbc][1] = vel[inb,jnb][1]
        rho[ibc,jbc] = rho[inb,jnb]
        for k in ti.static(range(9)):
            f_old[ibc,jbc][k] = f_eq(ibc,jbc,k,e,vel,w,rho) - f_eq(inb,jnb,k,e,vel,w,rho) + f_old[inb,jnb][k]

if show_plot == True:
    gui = ti.GUI('lbm solver', (nx, 2 * ny))

print(rho.shape)
print(tau)
print(bc_value[0,0])
print(bc_value[0,0]*cy_para[2]/niu)
tic = timeit.default_timer()
lbm_init()
for istep in range(steps):
    print(istep," f_new,f_old before collide_and_stream",np.sum(f_new.to_numpy())," ",np.sum(f_old.to_numpy()))
    print(istep," rho,vel, before collide_and_stream ",np.sum(rho.to_numpy()),np.sum(vel.to_numpy()[:,:,0]),np.sum(vel.to_numpy()[:,:,1]))
    collide_and_stream()   
    print(istep," f_new,f_old after collide_and_stream ",np.sum(f_new.to_numpy())," ",np.sum(f_old.to_numpy()))
    print(istep," rho,vel, after collide_and_stream ",np.sum(rho.to_numpy()),np.sum(vel.to_numpy()[:,:,0]),np.sum(vel.to_numpy()[:,:,1]))
    
    print(istep," f_new,f_old before update_macro_vars ",np.sum(f_new.to_numpy())," ",np.sum(f_old.to_numpy()))
    print(istep," rho,vel, before update_macro_vars ",np.sum(rho.to_numpy()),np.sum(vel.to_numpy()[:,:,0]),np.sum(vel.to_numpy()[:,:,1]))
    update_macro_vars()
    print(istep," f_new,f_old after update_macro_vars ",np.sum(f_new.to_numpy())," ",np.sum(f_old.to_numpy()))
    print(istep," rho,vel, after update_macro_vars ",np.sum(rho.to_numpy()),np.sum(vel.to_numpy()[:,:,0]),np.sum(vel.to_numpy()[:,:,1]))

    apply_bc()
    print(istep,"\t",np.amax(rho.to_numpy()),"\t",np.amax(vel.to_numpy()))
    if (istep % 1000 == 0):
        elapsed = (timeit.default_timer()-tic)/1000
        print('Step: {timestep} \t {it_time} s'.format(timestep=istep, it_time = elapsed))
        # ti.imwrite((img[:,:,0:3]*255).astype(np.uint8), 'fig/karman_'+str(i).zfill(6)+'.png')
        tic = timeit.default_timer()
    if show_plot == True:
        #displaying vorticity
        veln = vel.to_numpy()
        ugrad = np.gradient(veln[:,:,0])
        vgrad = np.gradient(veln[:,:,1])
        vor = ugrad[1]-vgrad[0]
        veln_mag = (veln[:,:,0]**2.0 + veln[:,:,1]**2.0)**0.5
        ## color map
        colors = [(1, 1, 0), (0.953, 0.490, 0.016), (0, 0, 0),
            (0.176, 0.976, 0.529), (0, 1, 1)]
        my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            'my_cmap', colors)
        vor_img = cm.ScalarMappable(norm=matplotlib.colors.Normalize(
            vmin=-0.02, vmax=0.02),cmap=my_cmap).to_rgba(vor)
        vel_img = cm.plasma(veln_mag / 0.15)
        img = np.concatenate((vor_img, vel_img), axis=1)
        gui.set_image(img)
        gui.show()
    if (save_vtk == True) and (istep % save_interval == 0):
        filename = "vel-"+str(istep).zfill(8)+".vtk"
        with open(filename,"w+") as f:
            veln = vel.to_numpy()
            f.write("# vtk DataFile Version 2.0\n")
            f.write(filename+"\n")
            f.write('ASCII\n')
            f.write('DATASET STRUCTURED_POINTS\n')
            f.write('DIMENSIONS '+str(nx)+" "+str(ny)+" 1\n") #vtk requires data in 3d
            f.write('ORIGIN 0 0 0\n')
            f.write('SPACING 1 1 1\n')
            f.write('POINT_DATA '+str(nx*ny)+"\n")
            f.write('VECTORS velocity_field float\n')
            for j in range(veln.shape[1]): #vtk's are column major!!!
                for i in range(veln.shape[0]):
                    f.write(str(veln[i,j,0])+" "+str(veln[i,j,1])+" 0.0\n")

print("Done!")
