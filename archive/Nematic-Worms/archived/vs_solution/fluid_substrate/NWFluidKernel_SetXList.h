
#ifndef __FLUID_KERNEL__SET_X_LIST_H__
#define __FLUID_KERNEL__SET_X_LIST_H__

#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWFluidParameters.h"
#include "NWWormsParameters.h"

//.. finds the closes worm particle i to fluid particle id
__global__ void SetXListFluidKernel(
	float *flx, float *fly, 
	float *wx, float *wy, 
	int *xlist)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id < dev_flParams._NFLUID)
	{
		//.. search through all worm particles (fix for better search)
		float _flx = flx[id], _fly = fly[id];
		float dx, dy, rr, minrr = dev_simParams._XBOX;
		int mini = -1;
		for (int i = 0; i < dev_Params._NPARTICLES; i++)
		{
			dx = wx[i] - _flx;
			dy = wy[i] - _fly;
			//DevicePBC(dx, dev_simParams._XBOX);
			//DevicePBC(dy, dev_simParams._YBOX);
			rr = dx*dx + dy*dy;

			//.. only pics over XCUT2 interaction range
			if (rr > _XCUT2) continue;
			if (rr < minrr)
			{
				minrr = rr;
				mini = i;
			}
		}
		xlist[id] = mini;
	}
}

#endif