#init_helper.py
import numpy as np
import warp as wp
from tqdm.auto import tqdm

class InitHelper():
    def __init__(self, parent_instance):
        # The parent instance is passed (i.e., LargeClass instance)
        # Automatically copy all attributes from the parent class to the helper class
        self.__dict__.update(parent_instance.__dict__)
    def print_sim_info(self):
        pass

    def init_bounds1(self):
        #circular bound
        equidistantThetaValues = np.zeros((self.nbounds))
        deltaTheta = 2*np.pi / self.nbounds
        for i in np.arange(self.nbounds):
            equidistantThetaValues[i] = i * deltaTheta
        #self.np*self.nworms + self.numSol + self.nbounds
        offset = int(self.np*self.nworms + self.numSol)
        for i in np.arange(self.nbounds):
            self.Phost[i+offset,0] = self.rwall * np.cos(equidistantThetaValues[i])+self.hxo2
            self.Phost[i+offset,1] = self.rwall * np.sin(equidistantThetaValues[i])+self.hyo2
            self.Phost[i+offset,2] = 0.0
            self.Thost[i+offset] = 3#self.Phost[i+offset,3] = 3.0


    def init_bounds2(self):
        pass
    def init_bounds3(self):
        pass


    def init_worms1(self):
        # CIRCLE
        thetanow = 5.0*np.pi
        dth = 0.0
        for iw in np.arange(self.nworms):
            worm_z_height = np.random.rand()*self.L
            thetanow += 2.0*dth
            for ip in np.arange(self.np):
                i = (iw)*self.np+(ip)
                r = self.a*thetanow
                dth = self.length0/r
                thetanow += dth
                x = r*np.cos(thetanow) + self.hxo2
                y = r*np.sin(thetanow) + self.hyo2
                self.Phost[i][0] = r*np.cos(thetanow) + self.hxo2
                self.Phost[i][1] = r*np.sin(thetanow) + self.hyo2
                self.Phost[i][2] = worm_z_height
                self.Thost[i] = 1 # ptype

        for iw in np.arange(self.nworms):
            if np.random.rand() < 0.5:
                savepos = np.zeros((self.np,3))
                for ip in range(self.np):
                    i = iw * self.np+ip
                    savepos[ip,:] = self.Phost[i,:]
                for ip in range(self.np):
                    i = iw*self.np + ip
                    self.Phost[i,:] = savepos[self.np-1-ip,:]
    
    def init_worms2(self):
        pass

    def init_worms3(self):
        pass

    def init_fluid1(self):
        #CIRCLE
        minDist = 1.5
        numPlaced = 0
        offset = int(self.np*self.nworms)
        pbar = tqdm(total = self.numSol)
        while numPlaced <= self.numSol:
            x = 2*self.rwall*np.random.rand()
            y = 2*self.rwall*np.random.rand()
            dx = x - self.hxo2
            dy = y - self.hyo2
            r = np.sqrt(dx*dx + dy*dy)
            if r <= 0.95*self.rwall:
                isInCardioid = True
            else:
                isInCardioid = False
            tooClose = False
            for i in np.arange(numPlaced):
                dist = np.sqrt((self.Phost[numPlaced+offset,0]-x)**2 + (self.Phost[numPlaced+offset,1]-y)**2)
                if dist < minDist:
                    tooClose = True
                    break
            if isInCardioid and not tooClose:
                #print(numPlaced,"/",int(self.numSol)," ",numPlaced/(np.pi*self.rwall**2))
                self.Phost[numPlaced + offset,0] = x
                self.Phost[numPlaced + offset,1] = y
                self.Phost[numPlaced + offset,2] = 0.0
                self.Thost[numPlaced+offset] = 2 #self.Phost[numPlaced + offset,3] = 2.0
                numPlaced += 1
                pbar.update(1)



    def init_fluid2(self):
        pass
    def init_fluid3(self):
        pass