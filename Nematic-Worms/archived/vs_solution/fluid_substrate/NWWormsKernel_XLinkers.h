#ifndef __WORMS_KERNEL__XLINKERS_H__
#define __WORMS_KERNEL__XLINKERS_H__
// -----------------------------------------------------------------------------------------
#include "NWmain.h"
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"

// -----------------------------------------------------------------------------------------
//	Count the current number of cross-links (launch one)
__global__ void XLinkerCountKernel(int *xlink, int *xcount){
	int new_count = 0;
	for (int i = 0; i < dev_Params._NPARTICLES; i++){
		if (xlink[i] != -1) new_count++;
	}
	(*xcount) = new_count;
}

// -----------------------------------------------------------------------------------------
// Randomly attempt to break some bonds
__global__ void XLinkerBreakKernel(int *xlink, float* randNum, float percentLoss)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){
		if (xlink[id] != -1){ // if linked
			if (randNum[id] < percentLoss){
				xlink[id] = -1;
			}
		}
	}
}

// -----------------------------------------------------------------------------------------
//	Update the cross linkers based on intended percentages
__global__ void XLinkerUpdateKernel(float *r,
							  int rshift,
							  int *xlink,
							  int *nlist,
							  int nlshift,
							  float *randNum, // one per particle needed uniform
							  float change_percent,
							  float linkCutoff)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){
		if (change_percent < 0) { // remove certain percentage
			const float target_percent = abs(change_percent);
			if (xlink[id] != -1) { // if linked
				if (randNum[id] < target_percent){ // attempt remove xlink
					xlink[id] = -1; // remove link
				}
			}
		}
		// FIX SO LOWER NUMBER NEIGHBORS AREN'T PICKED FIRST !!!!
		else { // add certain percentage 
			if (xlink[id] == -1) { // if unlinked
				if (randNum[id] < change_percent) { // if attempting adding 
					float rid[_D_], dr[_D_], _r[_D_], _rr;
					for_D_ rid[d] = r[id + d*rshift];
					for_D_ dr[d] = 0.0f;
					for (int n = 0; n < dev_Params._NMAX; n++){ // find someone to link with
						int nid = nlist[id + n * nlshift]; // neighbor id
						if (nid == -1) break; // end of neighbors
						for_D_ _r[d] = r[nid + d*rshift];
						_rr = CalculateRR(rid, _r, dr); // distane sqrd
						if (_rr > linkCutoff*linkCutoff) continue; // if close enough
						xlink[id] = nid; // assign linkage
					}
				}
			}
		}
	}
}
// ----------------------------------------------------------------------------------------
__global__ void XLinkerForceKernel(float *f,
								   int fshift,
								   float *r,
								   int rshift,
								   int *xlink,
								   bool f_on_self) 
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){
		int xid = xlink[id]; // cross linked to
		if (xid != -1){ // if linked to someone
			float k = dev_Params._Kx;
			int pid;
			if (f_on_self) {
				pid = id; // apply force to self
				k = -k; // set force direction
			}
			else pid = xid; // apply force to particle linked to
			float rid[_D_], rnab[_D_], dr[_D_]; // local position vectors
			for_D_ { // for all dimensions
				rid[d] = r[id + d*rshift]; // assign values for id
				rnab[d] = r[xid + d*rshift]; // assign values for xid
			}
			const float _r = sqrt(CalculateRR(rid, rnab, dr)); // distance
			const float _f = -k * (_r - dev_Params._Lx) / _r; // magnitude of force
			for_D_ f[pid + d*fshift] += _f * dr[d]; // apply force component
		}
	}
} // --------------------------------------------------------------------------------------
#endif