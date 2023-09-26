// Functions for GPU Device
// 5.12.15
// Mike Varga

#ifndef _KERNELS_H
#define _KERNELS_H

#include "cuda_runtime.h"

/*
*	Global scope function to print the parameters defined in the simulation
*	based on what files have been included and if they've been memcpy-ed t0
*	the device.
*/
__global__ void CheckParametersOnDevice(){
#ifdef __SIMULATION_PARAMETERS_H__
	printf("\nSimulation Parameters:\n");
	printf("----------------------\n");
	printf("XBOX = %f\nYBOX = %f\n", dev_simParams._XBOX, dev_simParams._YBOX);
	printf("DT = %f\nNSTEPS = %i\n(inner) = %i\n", dev_simParams._DT, dev_simParams._NSTEPS, dev_simParams._NSTEPS_INNER);
	printf("PRATE = %i\nMAXFILE = %i\n", dev_simParams._FRAMERATE, dev_simParams._FRAMESPERFILE);
#endif

#ifdef __WORMS_PARAMETERS_H__
	printf("\nWorms Parameters:\n");
	printf("-----------------\n");
	printf("NP = %i\nNWORMS = %i\n", dev_Params._NP, dev_Params._NWORMS);
	printf("NPARTICLES = %i\n", dev_Params._NPARTICLES);
	printf("LISTSETGAP = %i\nNMAX = %i\n", dev_Params._LISTSETGAP, dev_Params._NMAX);
	printf("EPSILON = %f\nSIGMA = %f\n", dev_Params._EPSILON, dev_Params._SIGMA);
	printf("DRIVE = %f\n", dev_Params._DRIVE);
	printf("K1 = %f\nK2 = %f\nK3 = %f\nKa = %f\n", dev_Params._K1, dev_Params._K2, dev_Params._K3, dev_Params._Ka);
	printf("L1 = %f\nL2 = %f\nL3 = %f\n", dev_Params._L1, dev_Params._L2, dev_Params._L3);
	printf("KBT = %f\nGAMMA = %f\nDAMP = %f\n", dev_Params._KBT, dev_Params._GAMMA, dev_Params._DAMP);
	printf("SIG6 = %f\n2SIG6 = %f\nLJAMP = %f\n", dev_Params._SIGMA6, dev_Params._2SIGMA6, dev_Params._LJ_AMP);
	printf("RMIN = %f\nR2MIN = %f\n", dev_Params._RMIN, dev_Params._R2MIN);
	printf("RCUT = %f\nR2CUT = %f\n", dev_Params._RCUT, dev_Params._R2CUT);
	printf("BUFFER = %f\n", dev_Params._BUFFER);
	printf("LANDSCALE = %f\nXLINKERDENSITY = %f\n", dev_Params._LANDSCALE, dev_Params._XLINKERDENSITY);
	printf("Kx = %f\nLx = %f\n", dev_Params._Kx, dev_Params._Lx);
#endif

#ifdef __FLUID_PARAMETERS_H__
	printf("\nFluid Parameters:\n");
	printf("-----------------\n");
	printf("NFLUID = %i\n", dev_flParams._NFLUID);
	printf("LISTSETGAP = %i\nNMAX = %i\n", dev_flParams._LISTSETGAP, dev_flParams._NMAX);
	printf("EPSILON = %f\nSIGMA = %f\n", dev_flParams._EPSILON, dev_flParams._SIGMA);
	printf("KBT = %f\nGAMMA = %f\n", dev_flParams._KBT, dev_flParams._GAMMA);
	printf("SIG6 = %f\n2SIG6 = %f\nLJAMP = %f\n", dev_flParams._SIGMA6, dev_flParams._2SIGMA6, dev_flParams._LJ_AMP);
	printf("RMIN = %f\nR2MIN = %f\n", dev_flParams._RMIN, dev_flParams._R2MIN);
	printf("RCUT = %f\nR2CUT = %f\n", dev_flParams._RCUT, dev_flParams._R2CUT);
	printf("BUFFER = %f\n", dev_flParams._BUFFER);
#endif
}

#endif