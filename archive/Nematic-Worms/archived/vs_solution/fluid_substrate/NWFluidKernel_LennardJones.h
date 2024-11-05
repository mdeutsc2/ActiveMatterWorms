
#ifndef __FLUID_KERNEL__LENNARD_JONES_H__
#define __FLUID_KERNEL__LENNARD_JONES_H__

#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWFluidParameters.h"
#include "NWSimulationParameters.h"

__global__ void LJFluidNListKernel(float *fx, float *fy, float *x, float *y, int *nlist)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{

		int listid = id*dev_flParams._NMAX;
		float fxid = 0.0;
		float fyid = 0.0;

		//.. loop through neighbors
		for (int i = 0; i < dev_flParams._NMAX; i++)
		{
			//.. neighbor id
			int nid = nlist[listid + i];

			//.. if no more players
			if (nid == -1) break;

			float dx = x[nid] - x[id];
			float dy = y[nid] - y[id];
			DevicePBC(dx, dev_simParams._XBOX);
			DevicePBC(dy, dev_simParams._YBOX);
			float rr = dx * dx + dy * dy;

			//.. stop if too far
			if (rr > dev_flParams._R2CUT) continue;

			//.. calculate LJ force
			float ff = dev_flParams._LJ_AMP*((dev_flParams._2SIGMA6 / (rr*rr*rr*rr*rr*rr*rr)) - (1.00000f / (rr*rr*rr*rr)));
			fxid -= ff * dx;
			fyid -= ff * dy;
		}

		//.. assign tmp to memory
		fx[id] += fxid;
		fy[id] += fyid;
	}
}

#endif