
#ifndef __WORMS_KERNEL__SET_X_LIST_H__
#define __WORMS_KERNEL__SET_X_LIST_H__
// 2D only
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWSimulationParameters.h"
#include "NWWormsParameters.h"
#include "NWFluidParameters.h"

//.. tells that there is a particle near to drive it
__global__ void SetXListWormsKernel(
	float *wx, float *wy, 
	float *flx, float *fly, 
	int *xlist)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id < dev_Params._NPARTICLES)
	{
		//.. search through all fluid particles (fix for better search)
		float wxid = wx[id], wyid = wy[id];
		float dx, dy, rr, minrr = dev_simParams._XBOX;
		int mini = -1;
		for (int i = 0; i < dev_flParams._NFLUID; i++)
		{
			dx = flx[i] - wxid;
			dy = fly[i] - wyid;
			//DevicePBC(dx, dev_simParams._XBOX);
			//DevicePBC(dy, dev_simParams._YBOX);
			rr = dx*dx + dy*dy;

			//.. only sets when within XCUT2 interaction range
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
