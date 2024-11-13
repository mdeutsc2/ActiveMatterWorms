import numpy as np
import warp as wp
from tqdm.auto import tqdm
from kernels import *
from init_helper import InitHelper
from io_helper import IOHelper
from enum import Enum

class BDType(Enum):
    CIRCLE = 1
    CARDIOID = 2
    EPITROCHOID = 3

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
        self.L = 3.2
        self.np = 80
        self.nworms = 800
        print(self.np,self.nworms, self.np*self.nworms)
        self.nsteps = 2000
        self.rwall = 164
        self.dt = 5e-4
        self.save_interval = 200
        self.dimensions = (200,200,100)
        # FORCES
        self.rcut = 2.5
        self.rcutsmall = 2.0**(1.0/6.0)
        self.fdogic = 0.06
        self.dogic_fdep = 0.25
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
        self.fluid_rho = 0.05#0.2
        self.fluid_offset = 2.0 #rcutsmall*sigma,//3.0; // z-offset of fluid
        self.random_init = True # RSA fluid init
        self.fdrag = 0.01
        # BOUNDARY
        self.nbounds = 5000
        self.numSol = 0
        self.k = 2 # number of cusps for EPITROCHOID
        self.bd = 1 # boundary initial type
        if self.bd == 1:
            self.bd = BDType.CIRCLE
            print("disk area ",np.pi*self.rwall*self.rwall)
            self.numSol = int(np.ceil(self.fluid_rho * (np.pi*self.rwall**2)))
            print("# fluid ptc: ",self.numSol)
        elif self.bd == 2:
            self.bd = BDType.CARDIOID
            ca = 1.5*(self.rwall/2)
            print("ca ",ca)
            cardioid_area = (6 * np.pi * ca ** 2)
            print("cardoid area ",cardioid_area)
            self.numSol = int(np.ceil(self.fluid_rho*cardioid_area))
        elif self.bd == 3:
            self.bd = BDType.EPITROCHOID
            if self.k == 0:
                print("error k=",self.k)
                exit()
            cylc_area = 12*np.pi*((self.rwall/4)+1)**2
            self.numSol = int(np.ceil(self.fluid_rho*cylc_area))
        else:
            if self.numSol == 0:
                print("no boundary defined, numSol=",self.numSol)
                print(self.bd)
                exit()

        #derived parameters
        self.N_active = int(self.np*self.nworms)
        self.N_mobile = int(self.np*self.nworms + self.numSol)
        self.N = int(self.np*self.nworms + self.numSol + self.nbounds)
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
        self.F = wp.zeros(shape = self.N, dtype=wp.vec3, device=self.backend)
        self.Fold = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend)
        self.V = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend)
        self.P = wp.zeros(shape=self.N, dtype=wp.vec3, device=self.backend) # x,y,z,type
        self.T = wp.zeros(shape=self.N,dtype=wp.int32,device=self.backend )
        self.Phost = np.zeros((self.N,3)) # for initialization
        self.Thost = np.zeros(self.N,dtype=int) # particle types
        self.grid = wp.HashGrid(int(self.dimensions[0]/self.rcut),
                                int(self.dimensions[1]/self.rcut),
                                int(self.dimensions[2]/self.rcut),
                                device=self.backend) #adjust this sizing dynamically
        self.grid_cell_size = self.rcut
        self.init = InitHelper(self)
        self.io = IOHelper(self)

        
    #imported methods
    def write_xyzv(self,istep):
        return self.io.write_xyzv(istep)
    #CIRCLE
    def init_worms1(self):
        return self.init.init_worms1()
    def init_bounds1(self):
        return self.init.init_bounds1()
    def init_fluid1(self):
        return self.init.init_fluid1()
    # CARDOID
    def init_worms2(self):
        return self.init.init_worms2()
    def init_bounds2(self):
        return self.init.init_bounds2()
    def init_fluid2(self):
        return set.init.init_fluid2()
    # EPITROCHOID
    def init_worms3(self):
        return self.init.init_worms3()
    def init_bounds3(self):
        return self.init.init_bounds3()
    def init_fluid3(self):
        return self.init.init_fluid3()

    def init_arrays(self):
        if self.bd == BDType.CIRCLE:
            self.init_bounds1()
            self.init_worms1()
            self.init_fluid1()
        elif self.bd == BDType.CARDIOID:
            self.init_bounds2()
            self.init_worms2()
            self.init_fluid1()
        elif self.bd == BDType.EPITROCHOID:
            self.init_bounds3()
            self.init_worms3()
            self.init_fluid3()
        self.V = wp.from_numpy(np.random.rand(self.N,3)*0.1,dtype=wp.vec3,device=self.backend)
        self.P = wp.from_numpy(self.Phost,dtype=wp.vec3,device=self.backend)
        self.T = wp.from_numpy(self.Thost,dtype=wp.int32,device=self.backend)

    def write_xyzv_old(self,istep):
        ic = 0
        P_out = self.P.numpy()
        V_out = self.V.numpy()
        filename = f"amatter{istep:010d}.xyz"
        with open(filename, 'w') as xyzfile:
                # number of active particles + 4 edge-defining particles + boundary + solvent (optional)
                xyzfile.write(f"{int(self.np*self.nworms + self.nbounds + self.numSol + 4)}\n")
                #xyzfile.write(f"{self.N}\n")
                xyzfile.write("# 0\n")
                ic = 2
                for iw in np.arange(0,self.nworms):
                    ifirst = (iw)*self.np
                    dx = P_out[ifirst][0] - self.hxo2
                    dy = P_out[ifirst][1] - self.hyo2
                    xang = np.arctan2(dy, dx)
                    rx = -np.sin(xang)
                    ry = np.cos(xang)
                    ilast = (iw)*self.np+self.np-1
                    dot = (P_out[ifirst][0] - P_out[ilast][0]) * rx + (P_out[ifirst][1] - P_out[ilast][1]) * ry
                    if dot >= 0.0:
                        for i in range(0,self.np):
                            id = (iw)*self.np+i
                            xyzfile.write(f"A {P_out[id][0]} {P_out[id][1]} {P_out[id][2]} {0.0} {0.0} {0.0}\n")
                    else:
                        for i in range(0,self.np):
                            id = (iw)*self.np+i
                            xyzfile.write(f"B {P_out[id][0]} {P_out[id][1]} {P_out[id][2]} {0.0} {0.0} {0.0}\n")
                offset = int(self.np*self.nworms + self.numSol)
                for i in range(self.nbounds):
                    xyzfile.write(f"I {P_out[i+offset][0]} {P_out[i+offset][1]} 0.0 0.0 0.0 0.0\n")
                offset = int(self.np*self.nworms)
                for i in np.arange(self.numSol):
                    xyzfile.write(f"S {P_out[i+offset][0]} {P_out[i+offset][1]} 0.0 0.0 0.0 0.0\n")
                
                xyzfile.write(f"E {self.hxo2 - self.rwall} {self.hyo2 - self.rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {self.hxo2 - self.rwall} {self.hyo2 + self.rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {self.hxo2 + self.rwall} {self.hyo2 - self.rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {self.hxo2 + self.rwall} {self.hyo2 + self.rwall} 0.0 0.0 0.0 0.0\n")

    def run(self):
        print("Starting")
        for istep in tqdm(range(self.nsteps)):
            # REBUILD HASH GRID
            self.grid.build(self.P,self.grid_cell_size)
            # UPDATE POSITIONS
            wp.launch(kernel=update_positions,
                      dim = self.N_mobile,
                      inputs=[self.grid.id,self.P,self.V,self.F,self.Fold,self.dt,self.L,self.N_active,self.N_mobile],
                      device=self.backend)
            # LJ FORCES
            wp.launch(kernel=apply_forces,
                        dim=len(self.P),
                        inputs=[self.grid.id,self.F,self.P,self.V,self.T,self.rcut,self.rcutsmall,
                                self.nworms,self.np,self.fdep,self.kbt,self.dt,self.fdogic,self.dogic_fdep,
                                self.fdepwall,self.fluid_offset,self.fdrag],
                        device=self.backend)
            # BOND FORCES

            # UPDATE VELOCITIES
            wp.launch(kernel=update_velocities,
                      dim = self.N_mobile,
                      inputs=[self.grid.id,self.V,self.F,self.Fold,self.dt,self.N_active,self.N_mobile],
                      device=self.backend)
            if (istep % self.save_interval == 0):
                sim.write_xyzv_old(istep)

if __name__ == "__main__":
    sim = AmatterSim()
    sim.init_arrays()
    sim.write_xyzv_old(0)
    sim.run()
    # print(sim.P)