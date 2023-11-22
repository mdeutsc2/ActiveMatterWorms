
#ifndef __WORMS_KERNEL__LENNARD_JONES_H__
#define __WORMS_KERNEL__LENNARD_JONES_H__
// 2D
//-------------------------------------------------------------------------------------------------
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
//------------------------------------------------------------------------------------------------
__global__ void LennardJonesNListKernel(float *f,
										int fshift,
										float *r,
										int rshift,
										int *nlist,
										int nshift )
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){

		float fid[_D_], rid[_D_], dr[_D_], _r[_D_];
		float _f, _rr;
		const int np = dev_Params._NP;
		const float drive = dev_Params._DRIVE;
		int pnp = id % np; // particle in chain
		for_D_{ fid[d] = 0.0f; rid[d] = r[id + d*rshift];}

		bool extensile = dev_Params._EXTENSILE;

		//.. loop through neighbors
		int nnp, nid;
		for (int i = 0; i < dev_Params._NMAX; i++)
		{
			//.. neighbor id
			nid = nlist[id + i*nshift];
			nnp = nid % np;

			//.. if no more players
			if (nid == -1) break;

			for_D_ _r[d] = r[nid + d*rshift];
			_rr = CalculateRR(rid, _r, dr);

			//.. stop if too far
			if (_rr > dev_Params._R2CUT) continue;

			//.. calculate LJ force
			_f = CalculateLJ(_rr);

			//.. extensile driving mechanism applied here!!! need a better place
			if (extensile){
				if ((pnp < np - 1) && (nnp < np - 1)){
					float rnxt1[_D_], rnxt2[_D_], x[_D_], u[_D_], f_ext[_D_];
					for_D_ rnxt1[d] = r[(id + 1) + d*rshift];
					for_D_ rnxt2[d] = r[(nid + 1) + d*rshift];
					CalculateRR(rid, rnxt1, x); // vector connecting id to next particle
					CalculateRR(_r, rnxt2, u); // vector connecting nid to next particle

					//..  average of x and -u (to preserve detailed balance)
					for_D_ f_ext[d] = (x[d] - u[d]) / 2.0f;
					for_D_ fid[d] += (f_ext[d] * dev_Params._DRIVE) / sqrt(_rr); // try this for now
					//float dotprod = dot(x, u); // dot product
					//if (dotprod < -75.0f){ // if anti-parallel
					//	for_D_ fid[d] += dotprod * drive * 
					//}
				}
			}

			//.. apply forces to directions
			for_D_ fid[d] -= _f * dr[d];
		}

		//.. assign tmp to memory
		for_D_ f[id + d*fshift] += fid[d];

#ifdef __PRINT_FORCES__
		if (id == __PRINT_INDEX__)
			printf("\n\tLJ Kernel:\n\tf = { %f, %f, %f }", f[id], f[id + fshift], f[id + 2*fshift]);
#endif
	}
}
//---------------------------------------------------------------------------------------------------
//.. each block handles a particle and all interactions, a block for each particle basically
__global__ void FastLennardJonesNListKernel(float *f, int fshift,
											float *r, int rshift,
											int *nlist, int nlshift){
	int pid = blockIdx.x + gridDim.x * blockIdx.x; // particle index is total block id
	int nab = threadIdx.x + blockDim.x * threadIdx.y; // neighbor count
	const int nmax = dev_Params._NMAX;
	if ((pid < dev_Params._NPARTICLES) && (nab < nmax)){
		extern __shared__ float f_shared[];
		int nid = nlist[pid + nab * nlshift];
		if (nid != -1){
			float rid[_D_], rnab[_D_], dr[_D_], _f, rr;
			for_D_ rid[d] = r[pid + d*rshift];
			for_D_ rnab[d] = r[nid + d*rshift];
			rr = CalculateRR(rid, rnab, dr);
			if (rr < dev_Params._R2CUT){
				_f = CalculateLJ(rr);
				for_D_ f_shared[nab + d*nmax] -= _f * dr[d];
			}
		}
		else{
			for_D_ f_shared[nab + d*nmax] = 0.0f; // sets to zero if no partcle theres
		}

		__syncthreads(); // sync before adding them all up
		float fid[_D_];
		for_D_ {
			fid[d] = 0.0f;
			for (int n = 0; n < nmax; n++){
				fid[d] += f_shared[n + d*nmax];
			}
		}

		//.. apply to global forces
		for_D_ f[pid + d*fshift] += fid[d];
	}
}
#endif
