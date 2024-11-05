import numpy as np
import warp as wp
from tqdm.auto import tqdm
from kernels import *
from enum import Enum

class BDType(Enum):
    CIRCLE = 1
    CARDIOID = 2
    EPICYCLOID1 = 3
    EPICYCLOID2 = 4

class AmatterSim:
    def __init__(self):
        wp.config.verify_cuda = True
        wp.init()
        #INPUT
        self.name = "test1"
        self.chkpt = ""
        self.test = False
        self.slurm = False
        self.in_file = ""
        self.backend = "cpu" #"cuda"
        # DIMENSIONS
        self.np = 80
        self.nworms = 800
        self.nsteps = 2000
        self.rwall = 164
        self.dt = 5e-4
        self.save_interval = 20
        self.dimensions = (200,200,100)
        # FORCES
        self.rcut = 2.5
        self.rcutsmall = 2.0**(1.0/6.0)
        self.fdogic = 0.06
        self.walldrive = False
        self.fdogicweall = 0.001
        self.fdep = 0.25
        self.fdepwall = 6.0
        self.ww_epsilon = 0.5
        # THERMOSTAT
        self.kbt = 1.5
        self.gamma = 5.0
        self.thermo = True
        self.thermow = False
        # BONDS
        self.kspring = 57.146436
        self.k2spring = 50.0*self.kspring
        self.k3spring = 50.0*self.kspring
        self.kbend = 40.0
        self.length0 = 0.8
        #FLUID 
        self.fluid_cpl = True
        self.fluid_rho = 0.2
        self.fluid_offset = 2.0 #rcutsmall*sigma,//3.0; // z-offset of fluid
        self.random_init = True # RSA fluid init
        # BOUNDARY
        self.nbounds = 5000
        self.numSol = 0
        self.k = 2 # number of cusps for EPICYCLOID
        self.bd = 2 # boundary initial type
        if self.bd == 1:
            self.bd = BDType.CIRCLE
            print("disk area ",np.pi*self.rwall*self.rwall)
            self.numSol = np.ceil(self.fluid_rho * (np.pi*self.rwall**2))
        elif self.bd == 2:
            self.bd = BDType.CARDIOID
            ca = 1.5*(self.rwall/2)
            print("ca ",ca)
            cardioid_area = (6 * np.pi * ca ** 2)
            print("cardoid area ",cardioid_area)
            self.numSol = np.ceil(self.fluid_rho*cardioid_area)
        elif self.bd == 3:
            self.bd = BDType.EPICYCLOID1
            if self.k == 0:
                print("error k=",self.k)
                exit()
            cylc_area = 12*np.pi*((self.rwall/4)+1)**2
            self.numSol = np.ceil(self.fluid_rho*cylc_area)
        elif self.bd == 4:
            self.bd = BDType.EPICYCLOID2
            if self.k == 0:
                print("error k=",self.k)
                exit()
            cylc_area = 12*np.pi*((self.rwall/4)+1)**2
            self.numSol = np.ceil(self.fluid_rho*cylc_area)
            print("NOT DONE YET!")
            exit()
        else:
            if self.numSol == 0:
                print("no boundary defined, numSol=",self.numSol)
                print(self.bd)
                exit()
        #derived
        self.N_active = self.np*self.nworms
        self.N = self.np*self.nworms + self.numSol + self.nbounds
        self.hx = 2.0*self.rwall+1.0
        self.hy = self.hx
        self.hxo2 = self.hx/2
        self.hyo2 = self.hy/2
        self.inv_sqrt_dt = 1.0/np.sqrt(self.dt)
        self.a = 0.13
        self.dpd_ratio = 1.0
        self.sqrt_gamma_term = np.sqrt(2.0*self.kbt*self.gamma)
        self.restart_interval = self.save_interval*2

        # Arrays
        self.F = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend)
        self.Fold = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend)
        self.V = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend)
        self.P = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend)
        self.grid = wp.HashGrid(int(self.dimensions[0]/self.rcut),
                                int(self.dimensions[1]/self.rcut),
                                int(self.dimensions[2]/self.rcut),
                                device=self.backend) #adjust this sizing dynamically
        self.grid_cell_size = self.rcut
        

    def init_worms(self):
        tmp_pos = np.zeros((self.N,3))
        ireverse = np.zeros(self.N_active)
        #placing worm particles
        thetanow = 5.0*np.pi
        a = 0.14
        for iw in np.arange(self.nworms):
            if np.random.rand() < 0.5:
                ireverse[iw] = 1
            dth = 0.0
            for ip in np.arange(self.np):
                i = (iw)*self.np+(ip)
                r = a*thetanow
                dth = self.length0/r
                thetanow += dth
                tmp_pos[i][0] = r*np.cos(thetanow) + self.hxo2
                tmp_pos[i][1] = r*np.sin(thetanow) + self.hyo2
                tmp_pos[i][2] = 0.0
            thetanow += 2.0*dth

        for iw in np.arange(self.nworms):
            savex = np.zeros(self.np)
            savey = np.zeros(self.np)
            if ireverse[iw] == 1:
                for ip in range(self.np):
                    i = (iw)*self.np+(ip)
                    savex[ip] = tmp_pos[i][0]
                    savey[ip] = tmp_pos[i][1]
                for ip in range(self.np):
                    i = (iw)*self.np+(ip)
                    tmp_pos[i][0] = savex[ip]
                    tmp_pos[i][0] = savey[ip]
        self.P = wp.from_numpy(tmp_pos, dtype=wp.vec3, device=self.backend)

    def init_bounds(self):
        pass

    def init_fuild(self):
        pass
        
    def init_arrays(self):
        self.V = wp.from_numpy(np.random.rand(self.N,3)*0.1,dtype=wp.vec3,device=self.backend)

    def write_xyzv(self,istep):
        ic = 0
        filename = f"amatter{istep:010d}.xyz"
        with open(filename, 'w') as xyzfile:
                # number of active particles + 4 edge-defining particles + boundary + solvent (optional)
                xyzfile.write(f"{self.N}\n")
                xyzfile.write("# 0\n")
                ic = 2
                
                for iw in np.arange(0,self.nworms,self.np):
                    dx = self.P[iw][0] - self.hxo2
                    dy = self.P[iw][1] - self.hyo2
                    xang = np.atan2(dy, dx)
                    rx = -np.sin(xang)
                    ry = np.cos(xang)
                    i = (iw)*self.np+(self.np-1)
                    dot = (self.P[iw][0] - self.P[i][0]) * rx + (self.P[iw][1] - self.P[i][1]) * ry
                    
                    if dot >= 0.0:
                        for i in range(0,self.np):
                            id = (iw)*self.np+(i)
                            xyzfile.write(f"A {self.P[i][0]} {self.P[i][1]} {self.P[i][2]} {worms[iw][i][3]} {worms[iw][i][4]} {worms[iw][i][5]}\n")
                            ic += 1
                    else:
                        for i in range(np):
                            xyzfile.write(f"B {worms[iw][i][0]} {worms[iw][i][1]} {worms[iw][i][2]} {worms[iw][i][3]} {worms[iw][i][4]} {worms[iw][i][5]}\n")
                            ic += 1
                
                for i in range(numPoints):
                    xyzfile.write(f"I {bound[i][0]} {bound[i][1]} 0.0 0.0 0.0 0.0\n")
                    ic += 1
                
                if fluid_cpl:
                    for i in range(numSol):
                        xyzfile.write(f"S {solvent[i][0]} {solvent[i][1]} 0.0 {solvent[i][2]} {solvent[i][3]} 0.0\n")
                        ic += 1
                
                xyzfile.write(f"E {hxo2 - rwall} {hyo2 - rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {hxo2 - rwall} {hyo2 + rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {hxo2 + rwall} {hyo2 - rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {hxo2 + rwall} {hyo2 + rwall} 0.0 0.0 0.0 0.0\n")
                
                write_log(logfile, f"{filename}\t{ic} lines written", True)

    def run(self):
        for istep in tqdm(range(self.nsteps)):
            # REBUILD HASH GRID
            self.grid.build(self.P,self.grid_cell_size)
            # UPDATE POSITIONS
            wp.launch(kernel=update_positions,
                      dim = len(self.P),
                      inputs=[self.grid.id,self.P,self.V,self.F,self.Fold,self.dt,self.dimensions[0],self.dimensions[1],self.dimensions[2]],
                      device=self.backend)
            # LJ FORCES
            wp.launch(kernel=apply_forces,
                      dim=len(self.P),
                      inputs=[self.grid.id,self.F,self.P,self.V,self.rcut,self.nworms,self.np,self.fdep,self.kbt,self.dt],
                      device=self.backend)
            # BOND
            for bond_offset in np.arange(1,4):
                wp.launch(kernel=bond_bending_forces,
                          dim=(self.nworms,self.np-bond_offset),
                          inputs = [self.grid.id,self.F,self.P,bond_offset,self.np,self.kspring,self.length0],
                          device=self.backend)
            # wp.launch(kernel=bond_bending_forces,
            #           dim = self.nworms,
            #           inputs = [self.grid.id,self.F,self.P,2,self.np,self.k2spring,self.length0*2.0],
            #           device=self.backend)
            # wp.launch(kernel=bond_bending_forces,
            #           dim = self.nworms,
            #           inputs = [self.grid.id,self.F,self.P,3,self.np,self.k3spring,self.length0*3.0],
            #           device=self.backend)
            # UPDATE VELOCITIES
            wp.launch(kernel=update_velocities,
                      dim = len(self.P),
                      inputs=[self.grid.id,self.V,self.F,self.Fold,self.dt],
                      device=self.backend)




