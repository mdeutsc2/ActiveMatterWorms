
#ifndef __BOUNDARY_FUNCTIONS_H__
#define __BOUNDARY_FUNCTIONS_H__

#include "cuda_runtime.h"

__host__ void MovementBC(float &pos, float L)
{
	if (pos > L) pos -= L;
	if (pos < 0) pos += L;
}

__host__ void PBC(float &dr, float L)
{
	if (dr > L / 2.0) dr -= L;
	if (dr < -L / 2.0) dr += L;
}

#endif