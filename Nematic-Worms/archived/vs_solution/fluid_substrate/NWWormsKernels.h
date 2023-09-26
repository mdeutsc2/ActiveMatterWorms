#ifndef __WORMS_KERNELS_H__
#define __WORMS_KERNELS_H__
// 2D
// -----------------
//.. cuda includes
// -----------------
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "curand_kernel.h"
#include "math_functions.h"
#include <stdio.h>

// ------------------------------
//.. individual kernel includes
// ------------------------------
#include "NWWormsKernel_InternalForce.h"
#include "NWWormsKernel_Noise.h"
#include "NWWormsKernel_LennardJones.h"
#include "NWWormsKernel_DriveForce.h"
#include "NWWormsKernel_LandscapeForce.h"
#include "NWWormsKernel_CalculateTheta.h"
#include "NWWormsKernel_UpdateSystem.h"
#include "NWWormsKernel_SetNList.h"
#include "NWWormsKernel_SetXList.h"
#include "NWWormsKernel_BondBending.h"
#include "NWWormsKernel_XLinkers.h"

// --------------------
// Just for debugging
// --------------------

__global__ void AddForce(float *f, int fpitch, int dim, float force){
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){
		int fshift = fpitch / sizeof(float);
		int index = id + dim * fshift;
		f[index] += force;
	}
}

//// ---------------------------------------------------------------------
//// launched one thread per worm
//__global__ void BondBendingForces(float *f, 
//								  int fshift,
//								  float *r, 
//								  int rshift){
//	int wid = threadIdx.x + blockDim.x * blockIdx.x;
//	int nworms = dev_Params._NWORMS;
//	int np = dev_Params._NP;
//	const float k_a = dev_Params._Ka;
//	if (wid < nworms){
//		
//		//.. loop through particles in worm (excluding last)
//		int p2id, p3id;
//		float r1[3], r2[3], r3[3], f1[3], f2[3], f3[3];
//		float r12[3], r23[3];
//		float BOX[] = { dev_simParams._XBOX, dev_simParams._YBOX };
//		for (int p1id = wid * np; p1id < ((wid + 1)*np) - 2; p1id++){
//			
//			//.. particle ids
//			p2id = p1id + 1;
//			p3id = p2id + 1;
//
//			//.. grab memory and calculate distances
//			for (int d = 0; d < 3; d++){
//
//				//.. memory
//				r1[d] = r[p1id + d *rshift];
//				r2[d] = r[p2id + d *rshift];
//				r3[d] = r[p3id + d *rshift];
//
//				//.. distances
//				r12[d] = r2[d] - r1[d];
//				r23[d] = r3[d] - r2[d];
//
//				//.. PBC
//				if (d < 2) {
//					DevicePBC(r12[d], BOX[d]);
//					DevicePBC(r23[d], BOX[d]);
//				}
//
//			}
//
//			//.. calculate terms
//			float dot_r12_r23 = dot(r12, r23);
//			float r12r12 = dot(r12, r12);
//			float r23r23 = dot(r23, r23);
//			float mag12inv = 1.0f / mag(r12);
//			float mag23inv = 1.0f / mag(r23);
//			float a = k_a * mag12inv * mag23inv;
//			float A[] = { a, a, a };
//	
//			//..  calculate forces x, y, z
//			for (int d = 0; d < 3; d++){
//
//				//.. particle 1
//				f1[d] = A[d] * (r23[d] - ((dot_r12_r23 / r12r12) * r12[d]));
//
//				//.. particle 2
//				f2[d] = A[d] * (((dot_r12_r23 / r12r12) * r12[d]) - ((dot_r12_r23 / r23r23) * r23[d]) + r12[d] - r23[d]);
//
//				//.. particle 3
//				f3[d] = A[d] * (((dot_r12_r23 / r23r23) * r23[d]) - r12[d]);
//
//				//.. apply to containers
//				f[p1id + d * fshift] -= f1[d];
//				f[p2id + d * fshift] -= f2[d];
//				f[p3id + d * fshift] -= f3[d];
//			}
//		}
//	}
//}

#endif