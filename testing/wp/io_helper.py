#io_helper.py
import numpy as np

class IOHelper():
    def __init__(self, parent_instance):
        # The parent instance is passed (i.e., LargeClass instance)
        # Automatically copy all attributes from the parent class to the helper class
        self.__dict__.update(parent_instance.__dict__)

    def write_xyzv(self,istep):
        ic = 0
        P_out = self.P.numpy()
        print(self.P.numpy())
        exit()
        V_out = self.V.numpy()
        filename = f"amatter{istep:010d}.xyz"
        with open(filename, 'w') as xyzfile:
                # number of active particles + 4 edge-defining particles + boundary + solvent (optional)
                xyzfile.write(f"{int(self.np*self.nworms + self.nbounds + 4)}\n")
                #xyzfile.write(f"{self.N}\n")
                xyzfile.write("# 0\n")
                ic = 2
                for iw in np.arange(self.nworms):
                    ifirst = (iw)*self.np
                    dx = P_out[ifirst][0] - self.hxo2
                    dy = P_out[ifirst][1] - self.hyo2
                    xang = np.arctan2(dy, dx)
                    rx = -np.sin(xang)
                    ry = np.cos(xang)
                    ilast = (iw)*self.np+(self.np-1)
                    dot = (P_out[ifirst][0] - P_out[ilast][0]) * rx + (P_out[ifirst][1] - P_out[ilast][1]) * ry
                    if dot >= 0.0:
                        for i in range(0,self.np):
                            id = (iw)*self.np+(i)
                            xyzfile.write(f"A {P_out[i][0]} {P_out[i][1]} {P_out[i][2]} {0.0} {0.0} {0.0}\n")
                            ic += 1
                    else:
                        for i in range(0,self.np):
                            xyzfile.write(f"B {P_out[i][0]} {P_out[i][1]} {P_out[i][2]} {0.0} {0.0} {0.0}\n")
                            ic += 1
                offset = int(self.np*self.nworms + self.numSol)
                for i in range(self.nbounds):
                    xyzfile.write(f"I {P_out[i+offset][0]} {P_out[i+offset][1]} 0.0 0.0 0.0 0.0\n")
                
                # if fluid_cpl:
                #     for i in range(numSol):
                #         xyzfile.write(f"S {solvent[i][0]} {solvent[i][1]} 0.0 {solvent[i][2]} {solvent[i][3]} 0.0\n")
                #         ic += 1
                
                xyzfile.write(f"E {self.hxo2 - self.rwall} {self.hyo2 - self.rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {self.hxo2 - self.rwall} {self.hyo2 + self.rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {self.hxo2 + self.rwall} {self.hyo2 - self.rwall} 0.0 0.0 0.0 0.0\n")
                xyzfile.write(f"E {self.hxo2 + self.rwall} {self.hyo2 + self.rwall} 0.0 0.0 0.0 0.0\n")
