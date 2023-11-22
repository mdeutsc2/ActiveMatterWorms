
#ifndef __FLUID_KERNELS_H__
#define __FLUID_KERNELS_H__

// -------------------------------------------
//.. cuda includes
// -------------------------------------------
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "curand_kernel.h"
#include "math_functions.h"

// -------------------------------------------
//.. individual kernel includes and prototypes
// -------------------------------------------

__global__ void FluidNoiseKernel(float *fx, float *fy, float *vx, float *vy, float *x, float *y, float * randx);
#include "NWFluidKernel_Noise.h"

__global__ void LJFluidNListKernel(float *fx, float *fy, float *x, float *y, int *nlist);
#include "NWFluidKernel_LennardJones.h"

__global__ void UpdateFluidKernel(float *fx, float *fy, float *fx_old, float *fy_old, float *vx, float *vy, float *x, float *y);
#include "NWFluidKernel_UpdateSystem.h"

__global__ void SetFluidNList_N2Kernel(float *x, float *y, int *nlist);
#include "NWFluidKernel_SetNList.h"

__global__ void SetXListFluidKernel(float *wx, float *wy, float *flx, float *fly, int *xlist);
#include "NWFluidKernel_SetXList.h"

#endif