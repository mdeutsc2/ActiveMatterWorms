
#ifndef __FLUID_KERNEL__NOISE_H__
#define __FLUID_KERNEL__NOISE_H__

#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWFluidParameters.h"
#include "NWSimulationParameters.h"

__global__ void FluidNoiseKernel(
	float *fx, float *fy, 
	float *vx, float *vy, 
	float *x, float *y, 
	float * randx)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{
		float fxid = 0.0;
		float fyid = 0.0;

#ifdef _DRAG
		fxid -= dev_flParams._GAMMA * vx[id];
		fyid -= dev_flParams._GAMMA * vy[id];
#endif

#ifdef _NOISE
		float scalor = sqrtf(2.0f * dev_flParams._GAMMA * dev_flParams._KBT / dev_simParams._DT);
		fxid += scalor * randx[id];
		fyid += scalor * randx[id + dev_flParams._NFLUID];
#endif

		//.. assign temp fxid and fyid to memory
		fx[id] = fxid;
		fy[id] = fyid;
	}
}

#endif