
#ifndef __FORCE_XCHANGER_H__
#define __FORCE_XCHANGER_H__

#include "curand_kernel.h"

#include "NWmain.h"
#include "NWWorms.h"
#include "NWFluid.h"
#include "NWParams.h"
#include "NWForceXchangerKernels.h"
#include "NWRandomNumberGenerator.h"


class ForceXchanger {

	//.. for each worm particle stores exchange binding
	int * dev_xlist;

	//.. data on host when necessary
	int * xlist;

	//.. store any cuda runtime errors
	std::vector<cudaError_t> errorState;

	//.. random number generator
	GRNG &rng;
	//curandGenerator_t * curandGen;
	//float * curandNumsx;

	//.. block-thread structure
	int Threads_Per_Block;
	int Blocks_Per_Kernel;

	//.. references to things to exchange forces between
	Worms &worms;
	Fluid &fluid;

public:
	
	ForceXchanger(Worms &w, Fluid &f, GRNG &RNG);
	~ForceXchanger();

	void Init();
	void UpdateXList();
	void XchangeForces();
	void XListToHost();

	int* GetXList(){ return this->xlist; }

	void DisplayErrors();

private:

	void AllocateGPUMemory();
	void FreeGPUMemory();
	void AllocateHostMemory();
	void FreeHostMemory();
	void CheckSuccess(cudaError_t err);
	void CheckSuccess(curandStatus_t stat);
};

////////////////////////////// PUBLIC METHODS ///////////////////////////////////

ForceXchanger::ForceXchanger(Worms &w, Fluid &f, GRNG &RNG) 
	: worms(w), fluid(f), rng(RNG), dev_xlist(NULL)
	/*curandGen(new curandGenerator_t), curandNumsx(NULL),*/{
	this->Threads_Per_Block = f.Threads_Per_Block;
	this->Blocks_Per_Kernel = f.Blocks_Per_Kernel;
}

ForceXchanger::~ForceXchanger(){
	this->FreeGPUMemory();
	this->FreeHostMemory();
	//curandDestroyGenerator(*this->curandGen);
	//delete this->curandGen;
}

void ForceXchanger::Init(){
	//curandCreateGenerator(this->curandGen, CURAND_RNG_PSEUDO_DEFAULT);
	this->AllocateGPUMemory();
	this->AllocateHostMemory();
	CheckSuccess(cudaMemset(this->dev_xlist, -1, this->fluid.nfluid_int_alloc));
}

void ForceXchanger::UpdateXList(){
	//CheckSuccess(curandGenerateUniform(*this->curandGen, this->curandNumsx, 2 * nfluid_float_alloc));
	UpdateXchangeListKernel <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >> >(this->fluid.dev_X, this->fluid.dev_Y, this->worms.dev_X, this->worms.dev_Y, this->dev_xlist, this->rng.Get(2*this->fluid.parameters->_NFLUID)); 
	ErrorHandler(cudaGetLastError());
}

void ForceXchanger::XchangeForces(){
	XchangeForceKernel <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >>>(this->fluid.dev_Fx, this->fluid.dev_Fy, this->worms.dev_Fx, this->worms.dev_Fy, this->fluid.dev_X, this->fluid.dev_Y, this->worms.dev_X, this->worms.dev_Y, this->dev_xlist);
	ErrorHandler(cudaGetLastError());
}

void ForceXchanger::XListToHost(){
	CheckSuccess(cudaMemcpy(this->xlist, this->dev_xlist, this->fluid.nfluid_int_alloc, cudaMemcpyDeviceToHost));
}

void ForceXchanger::DisplayErrors(){
	if (this->errorState.size() > 0){
		printf("\nFORCE XCHANGER ERRORS:\n");
		for (int i = 0; i < this->errorState.size(); i++)
			printf("%i -- %s\n", i, cudaGetErrorString(this->errorState[i]));
		this->errorState.empty();
	}
}

/////////////////////////////// PRIVATE METHODS //////////////////////////////////

void ForceXchanger::AllocateGPUMemory(){
	CheckSuccess(cudaMalloc((void**)&(this->dev_xlist), this->fluid.nfluid_int_alloc));
	//CheckSuccess(cudaMalloc((void**)&(this->curandNumsx), 2 * nfluid_float_alloc));
}

void ForceXchanger::FreeGPUMemory(){
	CheckSuccess(cudaFree(this->dev_xlist));
	//CheckSuccess(cudaFree(this->curandNumsx));
}

void ForceXchanger::AllocateHostMemory(){
	this->xlist = new int[this->fluid.parameters->_NFLUID];
}

void ForceXchanger::FreeHostMemory(){
	delete[] this->xlist;
}

void ForceXchanger::CheckSuccess(cudaError_t err){
	if (err != cudaSuccess) this->errorState.push_back(err);
}

void ForceXchanger::CheckSuccess(curandStatus_t stat){
	if (stat == CURAND_STATUS_SUCCESS) return;
	else if (stat == CURAND_STATUS_ALLOCATION_FAILED)
		printf("\nAllocation of random numbers failed.\n");
	else if (stat == CURAND_STATUS_ARCH_MISMATCH)
		printf("\nArchitecture mismatch for random numbers.\n");
	else if (stat == CURAND_STATUS_DOUBLE_PRECISION_REQUIRED)
		printf("\nDouble precision requied for random numbers.\n");
	else if (stat == CURAND_STATUS_INITIALIZATION_FAILED)
		printf("\nInitialization of random numbers failed.\n");
	else if (stat == CURAND_STATUS_INTERNAL_ERROR)
		printf("\nInternal error generating random numbers.\n");
	else if (stat == CURAND_STATUS_LAUNCH_FAILURE)
		printf("\nLaunch failure generating random numbers.\n");
	else if (stat == CURAND_STATUS_LENGTH_NOT_MULTIPLE)
		printf("\nIncorrect length requested for random numbers.\n");
	else if (stat == CURAND_STATUS_NOT_INITIALIZED)
		printf("\nContainer not initialized before generation.\n");
	else if (stat == CURAND_STATUS_OUT_OF_RANGE)
		printf("\nOut of range error generating random numbers.\n");
	else if (stat == CURAND_STATUS_PREEXISTING_FAILURE)
		printf("\nPreexisting error caused failure generating random numbers.\n");
	else if (stat == CURAND_STATUS_TYPE_ERROR)
		printf("\nType error generating random numbers.\n");
	else if (stat == CURAND_STATUS_VERSION_MISMATCH)
		printf("\nVersion mismatch encountered generating random numbers.\n");
	printf("'ForceXchanger' instance.\n");
}

#endif