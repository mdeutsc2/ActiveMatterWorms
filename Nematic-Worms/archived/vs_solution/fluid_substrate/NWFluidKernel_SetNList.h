
#ifndef __FLUID_KERNEL__SET_N_LIST_H__
#define __FLUID_KERNEL__SET_N_LIST_H__

#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWFluidParameters.h"
#include "NWSimulationParameters.h"

__global__ void SetFluidNList_N2Kernel(float *x, float *y, int *nlist)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{
		int found = 0;
		float dx, dy, rr;
		for (int i = 0; i < dev_flParams._NFLUID; i++)
		{
			dx = x[i] - x[id];
			dy = y[i] - y[id];
			DevicePBC(dx, dev_simParams._XBOX);
			DevicePBC(dy, dev_simParams._YBOX);
			rr = dx*dx + dy*dy;

			if (rr > dev_flParams._R2CUT + dev_flParams._BUFFER) continue;
			if (i == id) continue;

			int plcid = id * dev_flParams._NMAX + found++;
			nlist[plcid] = i;
		}
	}
}

#endif