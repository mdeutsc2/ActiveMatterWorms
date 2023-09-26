
#ifndef __WORMS_KERNEL__CALCULATE_THETA_H__
#define __WORMS_KERNEL__CALCULATE_THETA_H__
// 2D
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"

__global__ void CalculateThetaKernel(float *r,
									int rshift,
									float *thphi,
									int tshift)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id < dev_Params._NPARTICLES){
		int p = id % dev_Params._NP;
		if (p < (dev_Params._NP - 1)){
			int pp1 = id + 1;
			float dr[3] = { 0, 0, 0 }; // always 3d
			for_D_ dr[d] = r[pp1 + d*rshift];
			BC_dr(dr, dev_simParams._BOX);
			const float rmag = mag(dr);
			thphi[id] = atan2f(dr[1] / rmag, dr[0] / rmag);
			thphi[id + tshift] = PI / 2.0f - asinf(dr[2] / rmag);
		}
	}
}

//.. launch one for each worm
__global__ void FinishCalculateThetaKernel(float *thphi, int tshift){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < dev_Params._NWORMS){

		// get to last particle in worm tid
		int id = (tid * dev_Params._NP) + (dev_Params._NP - 1);
		thphi[id] = thphi[id - 1];
		thphi[id + tshift] = thphi[(id - 1) + tshift];
	}
}

#endif
