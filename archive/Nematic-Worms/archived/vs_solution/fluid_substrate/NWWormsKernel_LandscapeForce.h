
#ifndef __WORMS_KERNEL__LANDSCAPE_FORCE_H__
#define __WORMS_KERNEL__LANDSCAPE_FORCE_H__
// 2D
#include "NWmain.h"
#include "NWParams.h"
#include "NWSimulationParameters.h"
#include "NWWormsParameters.h"
//old
/*__global__ void WormsLandscapeKernel(float *fx, float *fy, float *x, float *y){
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){
		
		float _KK = 2.0f * PI * _N / dev_simParams._XBOX;

		float dx = x[id] - (dev_simParams._XBOX / 2.0f);
		float dy = y[id] - (dev_simParams._YBOX / 2.0f);
		DevicePBC(dx, dev_simParams._XBOX);
		DevicePBC(dy, dev_simParams._YBOX);
		float r = sqrtf(dx * dx + dy * dy);
		float coskr = cosf(_KK * r);
		fx[id] += -_KK * (coskr / r) * dx;
		fy[id] += -_KK * (coskr / r) * dy;
	}

}*/

__global__ void WormsLandscapeKernel(float *f,
									 int fshift,
									 float *r,
									 int rshift){

	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){

		////.. position of giant attractor
		//const float eps = 1.0f;
		//const float sig = 20.0f;
		//const float X = dev_simParams._XBOX / 2.0f;
		//const float Y = dev_simParams._YBOX / 2.0f;
		//const float Z = 0.0f;
		
		//.. position vectors
		float rid[_D_], rnab[_D_], dr[_D_];
		for_D_ {
			rid[d] = r[id + d*rshift];
			rnab[d] = dev_simParams._BOX[d] / 2.0f; // attract to center		
		}
		
		CalculateRR(rid, rnab, dr); // calculate distance
		BC_dr(dr, dev_simParams._BOX); // boundary conditions
		for_D_ f[id + d*fshift] += dev_Params._LANDSCALE * dr[d]; // calc and apply forces

		//.. harmonic potential zeroed around z = 0
		//f[id + 2 * fshift] -= dev_Params._LANDSCALE * r[id + 2 * rshift];
		//
		//.. attraction to giant attractor at { Lx/2, Ly/2, 0 }
		//const float rr = CalculateRR_3d(rid, rnab, dr);
		//const float _f = CalculateLJ_3d(rr, sig, eps);
		//for (int d = 0; d < 3; d++) 
		//	f[id + d * fshift] -= _f * dr[d];
	}
}

#endif
