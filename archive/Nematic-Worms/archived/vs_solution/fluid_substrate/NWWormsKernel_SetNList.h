
#ifndef __WORMS_KERNEL__SET_N_LIST_H__
#define __WORMS_KERNEL__SET_N_LIST_H__
// 2D
//----------------------------------------------------------------------------------------------
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
//----------------------------------------------------------------------------------------------
__global__ void SetNeighborList_N2Kernel(float *r,
										 int rshift, 
										 int *nlist,
										 int nlshift,
										 int *cell,
										 int cshift)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){

		const int w1 = id / dev_Params._NP;
		int found = 0;
		int cid[_D_], ncid[_D_], dcid[_D_], ncell[_D_]; // cell indices occupied by particle
		for_D_ cid[d] = cell[id + d*cshift]; // assign local from global
		for_D_ ncell[d] = (int)(ceil(dev_simParams._BOX[d] / dev_Params._DCELL)); // calc number of cells

		float dr[_D_], rnab[_D_], rid[_D_], _rr; // position vectors
		for_D_ rid[d] = r[id + d*rshift]; // assign local from global

		const float cutoff = dev_Params._R2CUT + dev_Params._BUFFER;
		bool next;
		for (int p2 = 0; p2 < dev_Params._NPARTICLES; p2++){
			
			int w2 = p2 / dev_Params._NP;

			//.. stop if in same worm
			if (w2 == w1) continue;

			//next = false; // set to do calculate
			//for_D_ ncid[d] = cell[p2 + d*cshift]; // get neighbor cell indices
			//for_D_ dcid[d] = abs(cid[d] - dcid[d]); // mag of distance in each dim
			//for_D_ if (dcid[d] == (ncell[d] - 1)) dcid[d] = 1; // over the boundary fix
			//for_D_ if (dcid[d] > 1) next = true; // set flag if to far away
			//if (next) continue; // to next particle

			//.. add to nlist if within range
			for_D_ rnab[d] = r[p2 + d*rshift];
			_rr = CalculateRR(rid, rnab, dr);
			if ((_rr <= cutoff) && (found < dev_Params._NMAX)){
				int plcid = id + found++ * nlshift;
				nlist[plcid] = p2;
			}
		}
	}
}
// -------------------------------------------------------------------------------------------
__global__ void FastSetNeighborListKernel(float *r, int rshift,
	int nlist, int nlshift){
	int id = threadIdx.x + blockDim.x * blockIdx.x; // particle index
	int nid = threadIdx.y + blockDim.y * blockIdx.y; // target index
	int nparts = dev_Params._NPARTICLES;
	int np = dev_Params._NP;
	if (id < nparts && nid < nparts){
		int w1 = id / np;
		int w2 = nid / np;
		if (w1 != w2){
			//float rid[3], rnab[3];
			//.. unfinshed
		}


	}
}
#endif