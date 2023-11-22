
#ifndef __WORMS_KERNEL__UPDATE_SYSTEM_H__
#define __WORMS_KERNEL__UPDATE_SYSTEM_H__
// 2D
#include "NWmain.h"
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
// -------------------------------------------------------------------------------------
//.. Update positions and velocities of particles then save forces
//	 Update list of cell positions for neighbor finding
__global__ void UpdateSystemKernel(float *f,
								   int fshift,
								   float *f_old,
								   int foshift,
								   float *v,
								   int vshift, 
								   float *r,
								   int rshift,
								   int *cell,
								   int cshift,
								   float dt)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES)
	{
		//.. local components
		float dv[_D_], dr[_D_], rid[_D_], fid[_D_];
		for_D_ rid[d] = r[id + d*rshift];
		for_D_ fid[d] = f[id + d*fshift];

		//.. boundary conditions
		BC_r(fid, rid, dev_simParams._BOX);

		//.. change in velocity
		for_D_ dv[d] = 0.5f * (fid[d] + f_old[id + d*foshift]) * dt;

		//.. change in position
		for_D_ dr[d] = v[id + d*vshift] * dt + 0.5f * f_old[id + d*foshift] * dt * dt;

		//.. save forces
		for_D_ f_old[id + d*foshift] = fid[d];

		//.. update positions
		for_D_ rid[d] += dr[d];

		//.. boundary conditions and apply new pos
		//BC_r(rid, dev_simParams._BOX);
		for_D_ r[id + d*rshift] = rid[d];

		//.. update velocities
		for_D_ v[id + d*vshift] += dv[d];

		//.. update cell address
		for_D_ cell[id + d*cshift] = (int)(rid[d] / dev_Params._DCELL);
	}
}


// ---------------------------------------------------------------------------------------
//.. Update positions and velocities of particles then save forces
//	 Update list of cell positions for neighbor finding
__global__ void FastUpdateKernel(float *f, int fshift,
								 float *f_old, int foshift,
								 float *v, int vshift,
								 float *r, int rshift,
								 int *cell, int cshift,
								 float dt){
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int bid = blockIdx.y;
	if (tid < dev_Params._NPARTICLES){

		const int fid = tid + bid * fshift;
		const int foid = tid + bid * foshift;
		const int vid = tid + bid * vshift;
		const int rid = tid + bid * rshift;
		const int cid = tid + bid * cshift;

		//.. boundary conditions
		BC_r(f[fid], r[rid], dev_simParams._BOX[bid], bid); // only applies to

		float dvx = 0.5f * (f[fid] + f_old[foid]) * dt;
		float dx = v[vid] * dt + 0.5f * f_old[foid] * dt * dt;
		f_old[foid] = f[fid];
		r[rid] += dx;
		v[vid] += dvx;

		//.. boundary conditions
		//BC_r(r[rid], dev_simParams._BOX[bid]); // only applies to

		//if (bid == 0) BC_r(r[rid], dev_simParams._XBOX);
		//else if (bid == 1) BC_r(r[rid], dev_simParams._YBOX);
		//else if (bid == 2) BC_r(r[rid], dev_simParams._ZBOX);
		
		//.. update cell list
		//cell[cid] = (int)(r[rid] / dev_Params._DCELL);

		f[fid] = 0.0f;
	}
}

#endif
