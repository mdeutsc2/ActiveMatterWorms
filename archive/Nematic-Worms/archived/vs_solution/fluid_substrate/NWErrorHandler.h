
#ifndef __ERROR_HANDLER_H__
#define __ERROR_HANDLER_H__

#include "cuda_runtime.h"
#include <stdlib.h>
#include <stdio.h>

#define _LOSS_TOLERANCE 100

void ErrorHandler(cudaError_t Status)
{
	if (cudaSuccess != Status)
	{
		fprintf(stderr, "\nfail msg:  '%s'\n", cudaGetErrorString(Status));
		//fprintf(stderr, "\npress any key to clear memory...");
		//std::cin.get();
		//ShutDownGPUDevice();
		//CleanUpHost();
		abort();
	}
}

#endif