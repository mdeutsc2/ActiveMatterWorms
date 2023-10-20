#!/usr/local/anaconda3/bin/python3
#amatter.py

import os,sys
import timeit
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from dataclasses import dataclass
import math
from datetime import date

save_xyz = False
save_interval = 2000
output_steps = 1

@dataclass
class Params:
    np: int
    nworms: int
    nsteps: int
    fdogic: float
    fdogicwall: float
    fdep: float
    fdepwall: float
    diss: float
    dt: float
    kspring: float
    kbend: float
    length0: float

params = Params(400, # np
                250, #nworms
                10000, #nsteps
                0.06, #fdogic
                0.0, #fdogicwall
                1.0, #fdep
                0.0, #fdepwall,
                0.08, #diss
                0.02, #dt
                57.146436, #kspring
                40.0, #kbend
                0.8 #length0
                )

x = np.zeros((params.nworms,params.np))
y = np.zeros((params.nworms,params.np))
vx = np.zeros((params.nworms,params.np))
vy = np.zeros((params.nworms,params.np))
vxave = np.zeros((params.nworms,params.np))
vyave = np.zeros((params.nworms,params.np))
fx = np.zeros((params.nworms,params.np))
fy = np.zeros((params.nworms,params.np))
fxold = np.zeros((params.nworms,params.np))
fyold = np.zeros((params.nworms,params.np))
savex = np.zeros(params.np)
savey = np.zeros(params.np)

#integer arrays
ireverse = np.zeros(params.nworms)
ddx = np.zeros(params.nworms)
ddy = np.zeros(params.nworms)
hhead = np.zeros(504*504)
ipointto = np.zeros(params.nworms*params.np)
nnab = np.zeros((params.nworms,params.np))

rcut = 2.5
r2cut = rcut*rcut
rcutsmall = 2.0**(1.0/6.0)
rwall = 125.0*rcutsmall*math.sqrt(2.)
r2cutsmall = rcutsmall*rcutsmall
r2inside = (rwall - rcutsmall)**2
thetanow = 5.0*np.pi
a = 0.24
rmin = a*thetanow
density = (params.nworms*params.np)/(np.pi*rwall**2)

print("rwall:",rwall)
print("nworms:",params.nworms,"\t np:",params.np)
print("density:",density)
filename = "amatter_"+str(date.today().month) + "-" + str(date.today().day) + "-"+str(date.today().year)+".xyz"
print(filename)
hx = 2.0*rwall+1.0
hy = hx
hyo2 = 0.5*hy
hxo2 = 0.5*hx

nxcell = int((hx/rcutsmall)-1)
nycell = int((hy/rcutsmall)-1)
dcell = hx/float(nxcell)
ncells = nxcell*nycell
print("nxcells:",nxcell,"\tnycells:",nycell,"\tncells",ncells)

gnoise = 0.8/math.sqrt(10)*0.8
dt2o2 = params.dt*params.dt*0.5
dto2 = params.dt*0.5

length2 = 2*params.length0
lengthmax = params.length0*float(params.np-1)
iwalldrive = True

ddx[1]=1
ddy[1]=0
ddx[2]=1
ddy[2]=1
ddx[3]=0
ddy[3]=1
ddx[4]=-1
ddy[4]=1
ddx[5]=-1
ddy[5]=0
ddx[6]=-1
ddy[6]=-1
ddx[7]=0
ddy[7]=-1
ddx[8]=1
ddy[8]=-1
ddx[9]=0
ddy[9]=0

#PYTHON FUNCTIONS
def init_worms():
    thetanow = 5.0*np.pi
    for iw in range(params.nworms):
        ireverse[iw] = 0
        if np.random.rand() <= 0.5:
            ireverse[iw] = 1
        for i in range(params.np):
            r = a*thetanow
            dth = params.length0/r
            thetanow += dth
            x[iw,i] = hxo2+r*np.cos(thetanow)
            y[iw,i] = hyo2+r*np.sin(thetanow)
            xangle = math.atan2(y[iw,i]-hyo2,x[iw,i]-hxo2)
            #give initial velocity here
            vx[iw,i] = 0.0
            vy[iw,i] = 0.0
            fx[iw,i] = 0.0
            fy[iw,i] = 0.0
        thetanow += 4.0*dth
    
    for iw in range(params.nworms):
        if ireverse[iw] == 1:
            for i in range(params.np):
                savex[i] = x[iw,i]
                savey[i] = y[iw,i]
            for i in range(params.np):
                x[iw,i] = savex[params.np-1-i]
                y[iw,i] = savey[params.np-1-i]
    print("init_worms successful")
    
def write_xyz():
    if save_xyz:
        with open(filename,"a+") as fxyz:
            fxyz.write("# 0")
            xn = x.to_numpy()
            yn = y.to_numpy()
            for iw in range(params.nworms):
                dx = xn[iw,1] - hxo2
                dy = xn[iw,1] - hyo2
                rx = -np.sin(math.atan2(dy,dx))
                ry = np.cos(math.atan2(dy,dy))
                dot=(xn[iw,0]-xn[iw,params.np-1])*rx+(yn[iw,0]-yn[iw,params.np-1])*ry
                if dot >= 0.0:
                    for i in range(params.np):
                        fxyz.write("A "+str(xn[iw,i])+" "+str(yn[iw,i]) + " 0.0")
                else:
                    for i in range(params.np):
                        fxyz.write("B "+str(xn[iw,i])+" "+str(yn[iw,i]) + " 0.0")
            fxyz.write("E "+str(hxo2-rwall) + " " + str(hyo2-rwall) + " 0.0")
            fxyz.write("E "+str(hxo2-rwall) + " " + str(hyo2+rwall) + " 0.0")
            fxyz.write("E "+str(hxo2+rwall) + " " + str(hyo2-rwall) + " 0.0")
            fxyz.write("E "+str(hxo2+rwall) + " " + str(hyo2+rwall) + " 0.0")

def check_cells():
    for i in range(ncells):
        hhead[i] = -1

    for iw in range(params.nworms):
        for ip in range(params.np):
            ii = (iw-1)*params.np+ip
            icell = 1+np.floor(x[iw,ip]/dcell)
            jcell = 1+np.floor(y[iw,ip]/dcell)
            if icell > nxcell or icell < 1:
                print("nxcell=",nxcell)
                print("icell out of bounds ",iw,ip,x[iw,ip])
                print("icell=",icell,"\tdcell=",dcell)
                exit()
            if jcell > nycell or jcell < 1:
                print("nycell=",nycell)
                print("jcell out of bounds ",iw,ip,y[iw,ip])
                print("jcell=",jcell,"\tdcell=",dcell)
                exit()
            scell = int(icell + (jcell-1)*nxcell)
            if scell > ncells or scell < 0:
                print("scell out of bounds")
            
            ipointto[ii] = hhead[scell]
            hhead[scell] = ii

def assign_cells():
    for icell in range(nxcell):
        for jcell in range(nycell):
            #print("icell:",icell," jcell:",jcell)
            scell = int(icell + (jcell-1)*nxcell)
            if hhead[scell] != -1:
                #there are particles in the cell called scell so
                # lets check all the neighbor cells
                for idir in range(9):
                    icnab = icell+ddx[idir]
                    if icnab > nxcell or icnab <= 0:
                        break
                    jcnab = jcell + ddy[idir]
                    if jcnab > nycell or jcnab <= 0:
                        break

                    scnab = int(icnab + (jcnab-1)*nxcell)
                    if hhead[scnab] != -1:
                        #there are particles in the cell called scnab
                        ii = int(hhead[scell])
                        while ii > 0:
                            iworm = 1 + int((ii-1)/params.np)
                            if iworm == 0:
                                print("iworm == 0")
                                exit()
                            ip = int(ii-params.np*(iworm-1)) - 1
                            iworm = iworm - 1
                            #print(iworm,ip)
                            jj = int(hhead[scnab])

                            while jj > 0:
                                jworm = 1+int((jj-1)/params.np)
                                if jworm == 0:
                                    print("jworm == 0")
                                    exit()
                                jp = jj-params.np*(jworm-1) - 1
                                jworm = jworm - 1
                                #print(iworm,ip,jworm,jp)
                                inogo = 0
                                if (iworm == jworm) and (ip-jp < 2):
                                    inogo = 1 #on the same worm and close means no interaction
                                if (ii < jj) and inogo != 0:
                                    dddx = x[jworm,jp] - x[iworm,ip]
                                    dddy = y[jworm,jp] - y[jworm,jp]
                                    r2 = dddx*dddx + dddy*dddy
                                    riijj = np.sqrt(r2)

                                    #add attractive force between all pairs
                                    if r2 <= r2cutsmall:
                                        ffor = -48.0*r2**(-7) + 24.0*r2**(-4) + params.fdep/riijj
                                        ffx = ffor*dddx
                                        ffy = ffor*dddy
                                        fx[iworm,ip] = fx[iworm,ip] + ffx
                                        fx[jworm,jp] = fx[jworm,jp] - ffx
                                        fy[iworm,ip] = fy[iworm,ip] + ffy
                                        fy[jworm,jp] = fy[jworm,jp] - ffy

                                        # take neighbors into account in calculations
                                        vxave[iworm,ip]=vxave[iworm,ip]+vx[jworm,jp]
                                        vyave[iworm,ip]=vyave[iworm,ip]+vy[jworm,jp]
                                        nnab[iworm,ip]=nnab[iworm,ip]+1
                                        vxave[jworm,jp]=vxave[jworm,jp]+vx[iworm,ip]
                                        vyave[jworm,jp]=vyave[jworm,jp]+vy[iworm,ip]
                                        nnab[jworm,jp]=nnab[jworm,jp]+1

                                        # add dogic drive to interacting pairs
                                        # first calculate unit vectors along each worm
                                        ip1=ip+1
                                        if ip1 <= params.np-1:
                                            dxi=x[iworm,ip1]-x[iworm,ip]
                                            dyi=y[iworm,ip1]-y[iworm,ip]
                                        else:
                                            dxi=x(iworm,ip)-x[iworm,ip-1]
                                            dyi=y(iworm,ip)-y[iworm,ip-1]

                                        jp1=jp+1
                                        if jp1 <= params.np-1:
                                            dxj=x[jworm,jp1]-x[jworm,jp]
                                            dyj=y[jworm,jp1]-y[jworm,jp]
                                        else:
                                            dxj=x[jworm,jp]-x[jworm,jp-1]
                                            dyj=y[jworm,jp]-y[jworm,jp-1]

                                        if dxi*dxj+dyi*dyj <= 0:
                                            #normalize those vectors to make them unit vectors
                                            ri = np.sqrt(dxi*dxi+dyi*dyi)
                                            dxi = dxi/ri
                                            dyi = dyi/ri

                                            rj = np.sqrt(dxj*dxj + dyj*dyj)
                                            dxj = dxj/rj
                                            dyj = dyj/rj

                                            #now they are both unit vectors. Find the direction for the force...
                                            dx = (dxi-dxj)/2.0
                                            dy = (dyi-dyj)/2.0

                                            r = np.sqrt(dx*dx+dy*dy)
                                            dx = dx/r
                                            dy = dy/r

                                            #add an extra attractive component where kinesin drive is present
                                            ffx = params.fdogic*dx + 0.7*dddx/riijj
                                            ffy = params.fdogic*dy + 0.7*dddy/riijj
                                            
                                            fx[iworm,ip]=fx[iworm,ip]+ffx
                                            fx[jworm,jp]=fx[jworm,jp]-ffx
                                            fy[iworm,ip]=fy[iworm,ip]+ffy
                                            fy[jworm,jp]=fy[jworm,jp]-ffy
                                jj = int(ipointto[jj])
                            ii = int(ipointto[ii])
  
def update_pos():
    for iw in range(params.nworms):
        for i in range(params.np):
            print(iw,i,x[iw,i],vx[iw,i],fx[iw,i],fxold[iw,i])
            x[iw,i] = x[iw,i] + vx[iw,i]*params.dt + fx[iw,i]*dt2o2
            y[iw,i] = y[iw,i] + vy[iw,i]*params.dt + fx[iw,i]*dt2o2
            fxold[iw,i] = fx[iw,i]
            fyold[iw,i] = fy[iw,i]

def update_vel():
    for iw in range(params.nworms):
        for i in range(params.np):
            vx[iw,i] = vx[iw,i] + dto2*(fx[iw,i] + fxold[iw,i])
            vy[iw,i] = vy[iw,i] + dto2*(fy[iw,i] + fyold[iw,i])
            vxave[iw,i] = vxave[iw,i]/nnab[iw,i]
            vyave[iw,i] = vyave[iw,i]/nnab[iw,i]

def calc_forces():
    # zeroing out force arrays and adding gaussian noise
    for iw in range(params.nworms):
        for i in range(params.np):
            rsq = 0.0
            v1 = 1.0
            while rsq >= 0.999 or rsq <= 0.001:
                #print("rsq>=0.999 or rsq<= 0.0001")
                v1 = 2.0*np.random.rand()-1.0
                v2 = 2.0*np.random.rand()-1.0
                rsq = v1*v1+v2*v2
            fac = np.sqrt(-2.0*np.log(rsq)/rsq)
            g1 = v1*fac*gnoise
            th = np.random.rand()*np.pi*2.0
            fx[iw,i] = g1*np.cos(th)
            fy[iw,i] = g1*np.cos(th)

    #first-neighbor strings
    for iw in range(params.nworms):
        for i in range(params.np-1): #serialized
            ip1 = i+1
            dx = x[iw,ip1] - x[iw,i]
            dy = y[iw,ip1] - y[iw,i]
            r = np.sqrt(dx*dx+dy*dy)
            ff = -params.kspring*(r-params.length0)/r
            ffx = ff*dx
            ffy = ff*dy
            fx[iw,ip1] = fx[iw,ip1] + ff*dx
            fx[iw,i] = fx[iw,i] - ff*dx
            fy[iw,i] = fy[iw,i] - ff*dy
            fy[iw,ip1] = fy[iw,ip1] + ff*dy

    #bond-bending terms
    for iw in range(params.nworms):
        for i2 in range(params.np-2): #serialized
            i3 = i2+1
            i4 = i2+2
            #print(i2,i3,i4)
            x2=x[iw,i2]
            y2=y[iw,i2]
            x3=x[iw,i3]
            y3=y[iw,i3]
            x4=x[iw,i4]
            y4=y[iw,i4]
            y23=y3-y2
            y34=y4-y3
            x23=x3-x2
            x34=x4-x3
            r23=np.sqrt(x23*x23+y23*y23)
            r34=np.sqrt(x34*x34+y34*y34)

            cosvalue = (x23*x34+y23*y34)/(r23*r34)
            # print(iw," ",x23," ",x34," ",y23," ",y34," ",r23," ",r34)
            if cosvalue >= 1.0:
                cosvalue = 1.0
            cr = 1.0-cosvalue*cosvalue
            if cr < 0.0:
                if math.isclose(cr,0.0,rel_tol=1e-9,abs_tol=1e-9):
                    cr = 0.0
                else:
                    print(cr)
                    exit()
                # print(x2,x3,x23)
                # print(x3,x4,x34)
                # print(y2,y3,y23)
                # print(y3,y4,y34)
                # print(1.0- cosvalue*cosvalue)
                # print("line378, neg sqrt")
                # exit()

            sinvalue = np.sqrt(cr)

            ff = -params.kbend*sinvalue/(r23*r34)
            dot = x23*x34+y23*y34
            fac = dot/(r23*r23)

            f2x = ff*(x34-fac*x23)
            f2y = ff*(y34-fac*y23)

            fac = dot/(r34*r34)
            f4x = ff*(fac*x34-x23)
            f4y = ff*(fac*y34-y23)
            f3x = -f2x-f4x
            f3y = -f2y-f4y

            fx[iw,i2] = fx[iw,i2] + f2x
            fy[iw,i2] = fy[iw,i2] + f2y

            fx[iw,i3] = fx[iw,i3] + f3x
            fy[iw,i3] = fx[iw,i3] + f3y
            
            fx[iw,i4] = fx[iw,i4] + f4x
            fy[iw,i4] = fx[iw,i4] + f4y


def worm_wall_int():
    for iw in range(params.nworms):
        for i in range(params.np):
            #dissipation proportional to v relative to local average
            fx[iw,i] = fx[iw,i] - params.diss*(vx[iw,i] - vxave[iw,i])
            fy[iw,i] = fy[iw,i] - params.diss*(vy[iw,i] - vyave[iw,i])
            #now that we have used them, zero them
            vxave[iw,i] = vx[iw,i]
            vyave[iw,i] = vy[iw,i]
            #calculate distance to the center
            nnab[iw,i] = 1
            dx = x[iw,i]-hxo2
            dy = y[iw,i]-hyo2
            r2 = (dx*dx+dy*dy)
            if r2 >= r2inside:
                th = math.atan2(dy,dx)
                xwall = hxo2+rwall*np.cos(th)
                ywall = hyo2+rwall*np.sin(th)
                dx = xwall-x[iw,i]
                dy = ywall-y[iw,i]
                rr2 = dx*dx+dy*dy
                ffor = -48.0*rr2**(-7) + 24.0*rr2**(-4)
                fx[iw,i]=fx[iw,i] + ffor*dx
                fy[iw,i] = fy[iw,i] + ffor*dy
                if iwalldrive:
                    ip1 = i+1
                    dxi = x[iw,i] - x[iw,i-1]
                    dyi = y[iw,i] - x[iw,i-1]
                    if ip1 <= params.np-1:
                        dxi = x[iw,ip1] - x[iw,i]
                        dyi = y[iw,ip1] - y[iw,i]
                    #make it a unit vector
                    ri = np.sqrt(dxi*dxi+dyi*dyi)
                    dxi = dxi/ri
                    dyi = dyi/ri
                    dxj = -np.sin(th)
                    dyj = np.cos(th)
                    #if the vectors are not anti-parallel, reverse the vector along the wall
                    if (dxi*dxj+dyi*dyj) > 0.0:
                        dxj = -dxj
                        dyj = -dyj
                    #if the two vectors have any component pointing in opposite directions
                    if (dxi*dxj+dyi*dyj) < 0.0:
                        #find the direction for the force...
                        dx = (dxi-dxj)/2.0
                        dy = (dyi-dyj)/2.0
                        #normalize along the direction vector
                        ri = np.sqrt(dx*dx+dy*dy)
                        dx = dx/ri
                        dy = dy/ri 
                        #turn on dogic drive on the wall
                        ffx = params.fdogicwall*dx
                        ffy = params.fdogicwall*dy
                        fx[iw,i] = fx[iw,i] + ffx
                        fy[iw,i] = fy[iw,i] + ffy

#INIT    
tic = timeit.default_timer()
init_worms()
write_xyz()
# MAIN LOOP
for itime in range(params.nsteps):
    t = float(itime)*params.dt
    if (itime % output_steps == 0):
        elapsed = (timeit.default_timer()-tic)/output_steps
        print('Step: {timestep} \t Avg.Time/step{it_time} s \t Sim. Time {st}'.format(timestep=itime, it_time = elapsed,st=t))
        tic = timeit.default_timer()
    update_pos()
    calc_forces()
    worm_wall_int()
    check_cells()
    assign_cells()
    update_vel()
    if itime%save_interval == 0:
        write_xyz()
    print("Done!")