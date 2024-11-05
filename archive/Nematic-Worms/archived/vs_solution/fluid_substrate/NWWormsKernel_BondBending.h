#ifndef __WORMS_KERNEL__BOND_BENDING_H__
#define __WORMS_KERNEL__BOND_BENDING_H__
// 2D
// ------------------------------------------------------------------------------------------
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
// ---------------------------------------------------------------------
// launched one thread per worm
__global__ void BondBendingForces(float *f,
	int fshift,
	float *r,
	int rshift){
	int wid = threadIdx.x + blockDim.x * blockIdx.x;
	int nworms = dev_Params._NWORMS;
	int np = dev_Params._NP;
	const float ka_ratio = dev_Params._Ka_RATIO;
	const float ka1 = dev_Params._Ka;
	const float ka2 = dev_Params._Ka2;
	if (wid < nworms){

		//.. chose ka to use (default to ka1)
		const float k_a = (wid <= nworms*ka_ratio ? ka1 : ka2);

		//.. loop through particles in worm (excluding last)
		int p2id, p3id;
		float r1[_D_], r2[_D_], r3[_D_], f1[_D_], f2[_D_], f3[_D_];
		float r12[_D_], r23[_D_];
		//float BOX[] = { dev_simParams._XBOX, dev_simParams._YBOX };
		for (int p1id = wid * np; p1id < ((wid + 1)*np) - 2; p1id++){

			//.. particle ids
			p2id = p1id + 1;
			p3id = p2id + 1;

			//.. grab memory and calculate distances
			for_D_ {

				//.. memory
				r1[d] = r[p1id + d *rshift];
				r2[d] = r[p2id + d *rshift];
				r3[d] = r[p3id + d *rshift];

				//.. distances
				r12[d] = r2[d] - r1[d];
				r23[d] = r3[d] - r2[d];
			}

			//.. boundary conditions
			BC_dr(r12, dev_simParams._BOX);
			BC_dr(r23, dev_simParams._BOX);

			//.. calculate terms
			float dot_r12_r23 = dot(r12, r23);
			float r12r12 = dot(r12, r12);
			float r23r23 = dot(r23, r23);
			float mag12inv = 1.0f / mag(r12);
			float mag23inv = 1.0f / mag(r23);
			float a = k_a * mag12inv * mag23inv;
			float A[] = { a, a, a }; // always 3d

			//..  calculate forces x, y, z
			for_D_ {

				//.. particle 1
				f1[d] = A[d] * (r23[d] - ((dot_r12_r23 / r12r12) * r12[d]));

				//.. particle 2
				f2[d] = A[d] * (((dot_r12_r23 / r12r12) * r12[d]) - ((dot_r12_r23 / r23r23) * r23[d]) + r12[d] - r23[d]);

				//.. particle 3
				f3[d] = A[d] * (((dot_r12_r23 / r23r23) * r23[d]) - r12[d]);

				//.. apply forces to all 3 particles
				f[p1id + d * fshift] -= f1[d];
				f[p2id + d * fshift] -= f2[d];
				f[p3id + d * fshift] -= f3[d];
			}
		}
	}
}

#endif
