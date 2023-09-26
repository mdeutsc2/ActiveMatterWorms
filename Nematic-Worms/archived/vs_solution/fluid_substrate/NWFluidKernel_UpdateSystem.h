
#ifndef __FLUID_KERNEL__UPDATE_SYSTEM_H__
#define __FLUID_KERNEL__UPDATE_SYSTEM_H__

#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWFluidParameters.h"

__global__ void UpdateFluidKernel(float *fx, float *fy, float *fx_old, float *fy_old, float *vx, float *vy, float *x, float *y)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_flParams._NFLUID)
	{

		//.. change in velocity
		float dvx = 0.5f * (fx[id] + fx_old[id]) * dev_simParams._DT;
		float dvy = 0.5f * (fy[id] + fy_old[id]) * dev_simParams._DT;

		//.. change in position
		float dx = vx[id] * dev_simParams._DT + 0.5f * fx_old[id] * dev_simParams._DT * dev_simParams._DT;
		float dy = vy[id] * dev_simParams._DT + 0.5f * fy_old[id] * dev_simParams._DT * dev_simParams._DT;

		//.. save forces
		fx_old[id] = fx[id];
		fy_old[id] = fy[id];

		//.. update positions and velocities
		x[id] += dx;
		y[id] += dy;
		vx[id] += dvx;
		vy[id] += dvy;

		//.. boundary conditions
		DeviceMovementPBC(x[id], dev_simParams._XBOX);
		DeviceMovementPBC(y[id], dev_simParams._YBOX);
	}
}

#endif