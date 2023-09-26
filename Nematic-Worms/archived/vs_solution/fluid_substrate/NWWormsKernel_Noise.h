#ifndef __WORMS_KERNELS__NOISE_KERNEL__
#define __WORMS_KERNELS__NOISE_KERNEL__
// 2D
// ----------------------------------------------------------------------------
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
// ----------------------------------------------------------------------------
__global__ void NoiseKernel(float *f, int fshift,
							float *rn, float mag)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int bid = blockIdx.y;
	if (tid < dev_Params._NPARTICLES){
		f[tid + bid * fshift] += mag * rn[tid + bid * dev_Params._NPARTICLES];
	}
}

#endif