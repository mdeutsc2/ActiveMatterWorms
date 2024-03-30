import taichi as ti
from datetime import date
import numpy
from tqdm import tqdm
import asyncio

#ti.init(arch=ti.cpu,advanced_optimization=False,debug=True,cpu_max_num_threads=1)
#ti.init(arch=ti.cpu,kernel_profiler=True)
ti.init(arch=ti.gpu)

vec = ti.math.vec3


@ti.dataclass
class Particle:
   p: vec # position
   m: ti.f32 # mass
   t: ti.int32 # type / id #, 1 = active, 2 = fluid, 3 = boundary
   v: vec    # velocity
   f: vec    # forces
   fold: vec # old forces


# parameters

dt = 0.01
dt2o2 = (dt*dt)/2
dto2 = dt/2
L = 3.2

save_interval = 500
nworms = 800
np = 80
n_active = nworms*np
n_bound = 5000
n_fluid = 0

n = n_active + n_bound + n_fluid

nsteps = 1000
kspring = 57.146436
k2spring = 50.0*kspring
k3spring = 75.0*kspring
length0 = 0.8
length2 = 2.0*length0
length3 = 3.0*length0
rwall = 150
fdep = 0.25
fdogic =0.06
dogic_fdep = 0.25
fdepwall = 6.0
hx = 2.0*rwall + 1.0
hy = hx
hyo2 = hy/2
hxo2 = hx/2

rcut = 2.5
r2cut = rcut*rcut
rcutsmall = 2**(1.0/6.0)
r2cutsmall = rcutsmall*rcutsmall

# then create a 1d StructField with n elements to gather all the particles in the computational domain
pf = Particle.field(shape=(n,)) # pf -> particle field, first n_active, then n_bound, then n_fluid

# structures for cells
grid_n = int(ti.ceil(rwall/rcut))
print(grid_n,grid_n)
ptc_count = ti.field(dtype=ti.i32,shape=(grid_n, grid_n),name="ptc_count")
column_sum = ti.field(dtype=ti.i32, shape=grid_n, name="column_sum")
prefix_sum = ti.field(dtype=ti.i32, shape=(grid_n, grid_n), name="prefix_sum")
particle_id = ti.field(dtype=ti.i32, shape=n, name="particle_id")
list_head = ti.field(dtype=ti.i32, shape=grid_n * grid_n)
list_cur = ti.field(dtype=ti.i32, shape=grid_n * grid_n)
list_tail = ti.field(dtype=ti.i32, shape=grid_n * grid_n)

#arrays for init
i_reverse = ti.field(ti.f32, shape=nworms)
save_x = ti.field(ti.f32,shape=np)
save_y = ti.field(ti.f32,shape=np)

@ti.kernel
def init_worms():
   #placing worm particles
   thetanow = 5.0*numpy.pi
   a = 0.14
   for iw in range(nworms):
      if ti.random(ti.f32) < 0.5:
         i_reverse[iw] = 1
      dth = 0.0
      for ip in range(np):
         i = (iw)*np+(ip)
         r = a*thetanow
         dth = length0/r
         thetanow += dth
         pf[i].p[0] = r*ti.cos(thetanow) + hxo2
         pf[i].p[1] = r*ti.sin(thetanow) + hyo2
         pf[i].p[2] = 0.0
         pf[i].v = 0.0
         pf[i].f = 0.0
         pf[i].fold = 0.0
         pf[i].t = 1
         pf[i].m = 1.0
      thetanow += 2.0*dth

   for iw in range(nworms):
      if i_reverse[iw] == 1:
         for ip in range(np):
            i = (iw)*np+(ip)
            save_x[ip] = pf[i].p[0]
            save_y[ip] = pf[i].p[1]
         for ip in range(np):
            i = (iw)*np+(ip)
            pf[i].p[0] = save_x[ip]
            pf[i].p[1] = save_y[ip]
   
   # placing boundary particles
   for ib in range(n_bound):
      theta = ib*(2*numpy.pi/n_bound)
      pf[n_active+ib].p[0] = rwall*ti.cos(theta) + hyo2
      pf[n_active+ib].p[1] = rwall*ti.sin(theta) + hxo2
      pf.p[2] = 0.0
      pf[ib].v = 0.0
      pf[ib].f = 0.0
      pf[ib].fold = 0.0
      pf[ib].t = 3
      pf[ib].m = 1.0
   print(pf[0].p,pf[0].p[:2])

@ti.kernel
def update_pos():
   for i in pf:
      if i < n_active: #(n_active = nworms*np)
         pf[i].f /= pf[i].m
         pf[i].p += pf[i].v*dt + pf[i].f*dt2o2
         pf[i].fold = pf[i].f
         pf[i].f = 0.0
         if pf[i].p[2] > L:
            pf[i].p[2] = L
         elif pf[i].p[2] < 0.0:
            pf[i].p[2] = L
      if i > n_active+n_bound:
         # half-velocity update for fluid particles
         pf[i].p += pf[i].v*dt + pf[i].f*dt2o2 
         pf[i].p[2] = 0.0
         pf[i].f /= pf[i].m
         pf[i].v += 0.5*dt*pf[i].f
         pf[i].v[2] = 0.0
         pf[i].f = 0.0

@ti.kernel
def update_vel():
    for i in pf:
        if i < n_active:
            pf[i].v += dto2*(pf[i].f + pf[i].fold)
        if i > n_active+n_bound:
            # other half-velocity update for fluid particles
            pf[i].f /= pf[i].m
            pf[i].v += 0.5*dt*pf[i].f

@ti.kernel
def intraworm_forces():
   #first set of springs nearest neighbor springs
   for I in ti.grouped(ti.ndrange(nworms,np-1)):
      i = (I[0])*np+(I[1])
      ip1 = i + 1
      #print(I,i,ip1)
   
      dr = pf[ip1].p - pf[i].p
      r = ti.sqrt(dr[0]*dr[0] + dr[1]*dr[1])
      ff = (-kspring*(r - length0)/r)*dr

      pf[i].f -= ff
      pf[ip1].f += ff 

   # bond bending terms
   for I in ti.grouped(ti.ndrange(nworms,np-2)):
      i = (I[0])*np+(I[1])
      ip2 = i + 2

      dr = pf[ip2].p - pf[i].p
      r = ti.sqrt(dr[0]*dr[0] + dr[1]*dr[1])
      ff = (-k2spring*(r - length2)/r)*dr

      pf[i].f -= ff
      pf[ip2].f += ff 

   # 3-spring bond-bending
   for I in ti.grouped(ti.ndrange(nworms,np-3)):
      i = (I[0])*np+(I[1])
      ip3 = i + 3

      dr = pf[ip3].p - pf[i].p
      r = ti.sqrt(dr[0]*dr[0] + dr[1]*dr[1])
      ff = (-k3spring*(r - length3)/r)*dr

      pf[i].f -= ff
      pf[ip3].f += ff 

@ti.kernel
def calc_forces(pf: ti.template()):
   ptc_count.fill(0)
   #print(-1,ptc_count.shape)
   #ti.loop_config(serialize=True)
   for i in range(n):
      grid_idx = ti.floor((pf[i].p[:2])/hx * grid_n, int)
      #print(grid_idx,pf[i].p[:2],ti.floor(pf[i].p[:2], int))
      ptc_count[grid_idx] += 1
      #print(grid_idx,pf[i].p[:2],ti.floor(pf[i].p[:2], int),ptc_count[grid_idx])
   #print(0)
   for i in range(grid_n):
      sum = 0
      for j in range(grid_n):
         sum += ptc_count[i, j]
      column_sum[i] = sum
   #print(1)
   prefix_sum[0, 0] = 0
   ti.loop_config(serialize=True)
   for i in range(1, grid_n):
      prefix_sum[i, 0] = prefix_sum[i - 1, 0] + column_sum[i - 1]
   #print(2)
   for i in range(grid_n):
      for j in range(grid_n):
         if j == 0:
            prefix_sum[i, j] += ptc_count[i, j]
         else:
            prefix_sum[i, j] = prefix_sum[i, j - 1] + ptc_count[i, j]

         linear_idx = i * grid_n + j

         list_head[linear_idx] = prefix_sum[i, j] - ptc_count[i, j]
         list_cur[linear_idx] = list_head[linear_idx]
         list_tail[linear_idx] = prefix_sum[i, j]
   #print(3)
   for i in range(n):
      grid_idx = ti.floor(pf[i].p[:2]/hx * grid_n, int)
      linear_idx = grid_idx[0] * grid_n + grid_idx[1]
      ptc_location = ti.atomic_add(list_cur[linear_idx], 1)
      particle_id[ptc_location] = i
   """
   # brute force
   for i in range(n):
      for j in range(i+1,n):
         resolve_forces(i,j)
   """
   # fast
   for i in range(n):
      grid_idx = ti.floor(pf[i].p[:2]/hx * grid_n, int)
      x_begin = max(grid_idx[0] - 1, 0)
      x_end = min(grid_idx[0] + 2, grid_n)

      y_begin = max(grid_idx[1] - 1, 0)
      y_end = min(grid_idx[1] + 2, grid_n)

      for neigh_i in range(x_begin, x_end):
         for neigh_j in range(y_begin, y_end):
            neigh_linear_idx = neigh_i * grid_n + neigh_j
            for p_idx in range(list_head[neigh_linear_idx],list_tail[neigh_linear_idx]):
               j = particle_id[p_idx]
               if i < j:
                  pass
                  #resolve_forces(i, j)

@ti.func
def resolve_forces(i,j):
   # worm-worm
   if pf[i].t == 1 and pf[j].t == 1:
      iworm = int(i/np)
      ip = i - np*iworm
      jworm = int(i/np)
      jp = j - np*jworm
      if (iworm != jworm) and (ti.abs(ip-jp) > 2):
         dr = pf[j].p - pf[i].p
         r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
         if r2 < r2cutsmall:
            riijj = ti.sqrt(r2)
            ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0) + fdep/riijj
            ff = ffor*dr
            pf[i].f += ff
            pf[j].f -= ff
            # TODO put filament thermostat here
            
            # add active shear to interacting pairs
            ip1 = ip + 1
            dri = vec(0.0,0.0,0.0)
            drj = vec(0.0,0.0,0.0)
            if (ip1 <= np):
               iip1 = iworm*np+ip1
               dri = pf[iip1].p - pf[i].p
            else:
               dri = pf[i].p - pf[i-1].p

            jp1 = jp+1
            if (jp1 <= np):
               jjp1 = jworm*np + jp1
               drj = pf[jjp1].p - pf[j].p
            else:
               drj = pf[j].p - pf[j-1].p

            # if the two vectors have any component < 60 degrees
            if (ti.math.dot(dri,drj) <= -0.5):
               # normalize the vectors to make them unit vectors
               dri_norm = ti.math.normalize(dri)
               drj_norm = ti.math.normalize(drj)

               dr_direction = (dri_norm-drj_norm)/2.0
               dr_d_norm = ti.math.normalize(dr_direction)

               # add extra attractive component where kinesin drive is present
               ff = fdogic*dr_d_norm + dogic_fdep*dr/riijj
               pf[i].f += ff
               pf[j].f -= ff
            

   # worm-solvent

   # worm-boundary
   if pf[i].t == 1 and pf[j].t == 3:
      dr = pf[i].p - pf[j].p
      r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
      if (r2 <= r2cutsmall):
         r = ti.sqrt(r2)
         ffor = (1/r)*6.0 + fdepwall/r
         ff = ffor*dr
         pf[i].f += ff
         # TODO add dogic walldrive

   if pf[i].t == 3  and pf[j].t == 1:
      dr = pf[i].p - pf[j].p
      r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
      if (r2 <= r2cutsmall):
         r = ti.sqrt(r2)
         ffor = (1/r)*6.0 + fdepwall/r
         ff = ffor*dr
         pf[j].f += ff
         #TODO add dogic walldrive

   # fluid boundary

def write_xyzv(istep):
   xyzfile = open("amatter"+str(istep).zfill(9)+".xyz", "w") 
   xyzfile.write(str(n_active+n_bound+n_fluid+4)+"\n")
   xyzfile.write("# "+str(istep)+"\n")
   for iw in range(nworms):
      i = iw*np
      xangle = ti.atan2(pf[i].p[1],pf[i].p[0])
      rx = -ti.sin(xangle)
      ry = ti.cos(xangle)
      dot = (pf[i].p[0] - pf[i+np-1].p[0])*rx + (pf[i].p[1] - pf[i+np-1].p[1])*ry
      if dot >= 0.0:
         for ip in range(np):
            i = (iw)*np+(ip)
            xyzfile.write("A "+str(pf[i].p[0])+
                          " "+str(pf[i].p[1])+
                          " "+str(pf[i].p[2])+
                          " "+str(pf[i].v[0])+
                          " "+str(pf[i].v[1])+
                          " "+str(pf[i].v[2])+"\n")
      else:
         for ip in range(np):
            i = (iw)*np+(ip)
            xyzfile.write("B "+str(pf[i].p[0])+
                          " "+str(pf[i].p[1])+
                          " "+str(pf[i].p[2])+
                          " "+str(pf[i].v[0])+
                          " "+str(pf[i].v[1])+
                          " "+str(pf[i].v[2])+"\n")

   for i in range(n_active,n_active+n_bound):
      xyzfile.write("I "+str(pf[i].p[0])+
                    " "+str(pf[i].p[1])+
                    " 0.0 0.0 0.0 0.0\n")


   for i in range(n_active+n_bound,n_active+n_bound+n_fluid):
      xyzfile.write("S "+str(pf[i].p[0])+" "+str(pf[i].p[1])+" 0.0 "+str(pf[i].v[0])+" "+str(pf[i].v[1])+" 0.0\n")
   
   xyzfile.write("E "+str(hxo2-rwall)+" "+str(hyo2-rwall)+" 0.0 0.0 0.0 0.0\n")
   xyzfile.write("E "+str(hxo2-rwall)+" "+str(hyo2+rwall)+" 0.0 0.0 0.0 0.0\n")
   xyzfile.write("E "+str(hxo2+rwall)+" "+str(hyo2-rwall)+" 0.0 0.0 0.0 0.0\n")
   xyzfile.write("E "+str(hxo2+rwall)+" "+str(hyo2+rwall)+" 0.0 0.0 0.0 0.0\n")
   xyzfile.close()

def main():
   print("=============== AMATTER3D ===============")
   print(date.today()," :Matt Deutsch, Kent State University")
   
   # save params to file
   #write_params();

   init_worms()

   #init_fluid(solvent, numSol);

   #restart_write(0);
   for itime in tqdm(range(nsteps)):
      if (itime % save_interval == 0) or itime == 0:
        write_xyzv(itime)
        #pass
      #   if (itime % io_interval == 0) {
      #       xt.stop();
      #       total_time += xt.elapsed();
      #       var out_str:string = "Step: "+(itime+restart_timestep):string+"\t"+
      #                 (io_interval/xt.elapsed()):string+"iter/s\tCalc:"+
      #                 ((ct.elapsed()/total_time)*100):string+"%\tIO:"+
      #                 ((wt.elapsed()/total_time)*100):string+" %\tElapsed:"+
      #                 total_time:string+" s\t Est:"+
      #                 ((nsteps-itime)*(total_time/itime)):string+" s";
      #       write_log(logfile,out_str);
      #       //writeln((+ reduce solvent.vx),"\t",(+ reduce solvent.vy));
      #       //writeln((+ reduce solvent.x),"\t",(+ reduce solvent.y));
      #       xt.restart();
      #       }
      update_pos()

      intraworm_forces()

      #print("calculating forces")
      calc_forces(pf)
      #   if (fluid_cpl) {
      #       KEsol_total[itime] = (+ reduce KEsol);
      #       AMsol_total[itime] = (+ reduce AMsol);
      #   } else {
      #       KEsol_total[itime] = 0.0;
      #   }

      update_vel()
      #   KEworm_total[itime] = (+ reduce KEworm);
      #   KEworm_local_total[itime] = (+ reduce KEworm_local);
      #   AMworm_total[itime] = (+ reduce AMworm);

      #   if (itime % nstepso1e5 == 0){
      #       write_macro(macro_filename,itime);
      #   }
      #   if (itime % save_interval == 0){
      #       write_xyzv(itime+restart_timestep);
      #       if (itime % restart_interval == 0) {
      #           restart_write(itime+restart_timestep);
      #       }
      #   }
   #print("Total Time:"+total_time:string+" s");
   #ti.profiler.print_kernel_profiler_info('trace')
   #ti.profiler.print_kernel_profiler_info()

if __name__ == "__main__":
   main()
