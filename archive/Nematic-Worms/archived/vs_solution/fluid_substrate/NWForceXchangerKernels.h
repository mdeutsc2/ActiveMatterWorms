
#ifndef __FORCE_X_CHANGER__KERNELS_H__
#define __FORCE_X_CHANGER__KERNELS_H__

#include "NWDeviceFunctions.h"
#include "NWParams.h"

__global__ void UpdateXchangeListKernel(
	float *flx, float *fly,
	float *wx, float *wy,
	int *xlist, float *randx)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{
		//.. if not bound, try to bind
		float rand = randx[id];
		//float rand2 = randx[id + _NFLUID];
		if (xlist[id] == -1)
		{
			//.. loop over worm particles (need better search)
			for (int i = 0; i < dev_Params._NPARTICLES; i++)
			{
				float dx = flx[id] - wx[i];
				float dy = fly[id] - wy[i];
				DevicePBC(dx, dev_simParams._XBOX);
				DevicePBC(dy, dev_simParams._YBOX);
				float rr = dx * dx + dy * dy;

				if (rr >= _XCUT2) continue;
				if (rand > _MAKEBINDING)
				{
					xlist[id] = i; // fluid id linked to particle i
					i = dev_Params._NPARTICLES;
				}
			}
		}
		//.. if bound, try to release or move
		else
		{
			if (rand < -_LOSSBINDING) xlist[id] = -1;
			//else if (rand2 < _WALKBINDING)
			//{
			//	xlist[id] += 1;
			//	if (xlist[id] >= _NP) // walked off the end
			//		xlist[id] = -1;
			//}
		}
	}
}

__global__ void EliminateXListDoublesKernel(int *xlist)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{
		bool found = false;
		for (int i = 0; i < dev_Params._NPARTICLES; i++)
		{
			//.. if fluid id is linked to particle i
			if (xlist[id] == i)
			{
				//.. if already link exists, break others
				if (found) xlist[id] = -1;
				found = true; // flag the link exists
			}
		}
	}
}

__global__ void XchangeForceKernel(
	float *flfx, float *flfy,
	float *wfx, float *wfy,
	float *flx, float *fly,
	float *wx, float *wy,
	int *xlist)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{
		//.. if worm particle id is bound to something
		if (xlist[id] != -1)
		{
			//.. calculate distence and force
			float dx = flx[id] - wx[xlist[id]];
			float dy = fly[id] - wy[xlist[id]];
			DevicePBC(dx, dev_simParams._XBOX);
			DevicePBC(dy, dev_simParams._YBOX);
			float r = sqrtf(dx * dx + dy * dy);
			float ff = -_KBIND * r;
			float ffx = ff * dx;
			float ffy = ff * dy;

			//.. apply force to fluid and worm particles
			wfx[xlist[id]] -= ffx;
			wfy[xlist[id]] -= ffy;
			flfx[id] += ffx;
			flfy[id] += ffy;
		}
	}
}

#endif