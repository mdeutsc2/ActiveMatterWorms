
#ifndef __FLUID_H__
#define __FLUID_H__

#include "curand_kernel.h"
#include "NWmain.h"
#include "NWFluidKernels.h"
#include "NWParams.h"
#include "NWRandomNumberGenerator.h"
#include "NWFluidParameters.h"
#include "NWSimulationParameters.h"
#include "NWBoundaryFunctions.h"

class ForceXchanger;
// ------------------------------------------------------------------------------------------
//	This class defines a set of Lennard-Jones particles
// ------------------------------------------------------------------------------------------
//.. data structure for LJ fluid
class Fluid {

	//.. on cuda device
	float * dev_X;
	float * dev_Y;
	float * dev_Vx;
	float * dev_Vy;
	float * dev_Fx;
	float * dev_Fy;
	float * dev_Fx_old;
	float * dev_Fy_old;
	
	//.. neighbors list within fluid
	int * dev_nlist;

	//.. container of errors
	std::vector<cudaError_t> errorState;

	//.. randon number generator state
	GRNG * rng;

	//.. parameters for fluid
	FluidParameters * parameters;

	//.. environmental parameters
	SimulationParameters * envirn;

	//.. block-thread structure
	int Threads_Per_Block;
	int Blocks_Per_Kernel;

	//.. allocation sizes
	size_t nfluid_int_alloc;
	size_t nfluid_float_alloc;

public:
	//.. on host (so print and setup on host is quicker)
	float * X;
	float * Y;
	float * Vx;
	float * Vy;

	//.. contstruct with thread-block structure
	Fluid();
	~Fluid();

	//.. User access to running simulation
	void Init(GRNG *, FluidParameters *, SimulationParameters *, int);
	void InternalForces();
	void LJForces();
	void Update();
	void DataDeviceToHost();
	void DataHostToDevice();
	void ResetNeighborsList();
	void ResetNeighborsList(int itime);
	void ClearAll();

	//.. for Debugging
	void DisplayNList();
	void DisplayErrors();

	//.. so the exchanger can do what it needs
	friend class ForceXchanger;

private:
	//.. methods for internal use
	void AllocateHostMemory();
	void FreeHostMemory();
	void AllocateGPUMemory();
	void FreeGPUMemory();
	void DistributeFluidOnHost();
	void AdjustDistribute();
	void ZeroHost();
	void ZeroGPU();
	void CheckSuccess(cudaError_t err);
	void CheckSuccess(curandStatus_t stat);
	void FigureBlockThreadStructure(int threadsPerBlock);
};
// -------------------------------------------------------------------------------------------
//	PUBLIC METHODS 
// -------------------------------------------------------------------------------------------
//Fluid::Fluid(int BPK, int TPB, GRNG &RNG, FluidParameters &fluidParameters) :
Fluid::Fluid() :
X(NULL), Y(NULL), Vx(NULL), Vy(NULL), 
//rng(RNG), parameters(fluidParameters),
dev_X(NULL), dev_Y(NULL), dev_Vx(NULL), dev_Vy(NULL),
dev_Fx(NULL), dev_Fy(NULL), dev_Fx_old(NULL), dev_Fy_old(NULL){

	//.. do necessities
	srand(time(NULL));
	/*this->Threads_Per_Block = TPB;
	this->Blocks_Per_Kernel = BPK;
	this->nfluid_int_alloc = fluidParameters._NFLUID * sizeof(int);
	this->nfluid_float_alloc = fluidParameters._NFLUID * sizeof(float);*/
}

Fluid::~Fluid(){
	this->ClearAll();
}

void Fluid::ClearAll(){
	this->FreeGPUMemory();
	this->FreeHostMemory();
	this->X = this->Y = NULL;
	this->Vx = this->Vy = NULL;
	this->dev_X = this->dev_Y = NULL;
	this->dev_Vx = this->dev_Vy = NULL;
	this->dev_Fx = this->dev_Fy = NULL;
	this->dev_Fx_old = this->dev_Fy_old = NULL;
	this->rng = NULL;
	this->parameters = NULL;
	this->envirn = NULL;
}

//.. initializes system on host and device
void Fluid::Init(GRNG * gaussianRandomNumberGenerator,
				FluidParameters * fluidParameters,
				SimulationParameters * simParameters, 
				int threadsPerBlock = 512){
	this->rng = gaussianRandomNumberGenerator;
	this->parameters = fluidParameters;
	this->envirn = simParameters;
	this->FigureBlockThreadStructure(threadsPerBlock);
	this->AllocateHostMemory();
	this->ZeroHost();
	this->AllocateGPUMemory();
	this->ZeroGPU();
	this->DistributeFluidOnHost();
	this->DataHostToDevice();
}

void Fluid::InternalForces(){
	FluidNoiseKernel <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >>>(this->dev_Fx, this->dev_Fy, this->dev_Vx, this->dev_Vy, this->dev_X, this->dev_Y, this->rng->Get(2*parameters->_NFLUID));
	ErrorHandler(cudaGetLastError());
}

void Fluid::LJForces(){
	LJFluidNListKernel <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >> >(this->dev_Fx, this->dev_Fy, this->dev_X, this->dev_Y, this->dev_nlist);
	ErrorHandler(cudaGetLastError());
}

void Fluid::Update(){
	UpdateFluidKernel <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >>>(this->dev_Fx, this->dev_Fy, this->dev_Fx_old, this->dev_Fy_old, this->dev_Vx, this->dev_Vy, this->dev_X, this->dev_Y);
	ErrorHandler(cudaGetLastError());
}

void Fluid::DataDeviceToHost(){
	CheckSuccess(cudaMemcpy(this->X, this->dev_X, this->nfluid_float_alloc, cudaMemcpyDeviceToHost));
	CheckSuccess(cudaMemcpy(this->Y, this->dev_Y, this->nfluid_float_alloc, cudaMemcpyDeviceToHost));
}

void Fluid::DataHostToDevice(){
	CheckSuccess(cudaMemcpy(this->dev_X, this->X, this->nfluid_float_alloc, cudaMemcpyHostToDevice));
	CheckSuccess(cudaMemcpy(this->dev_Y, this->Y, this->nfluid_float_alloc, cudaMemcpyHostToDevice));
	CheckSuccess(cudaMemcpy(this->dev_Vx, this->Vx, this->nfluid_float_alloc, cudaMemcpyHostToDevice));
	CheckSuccess(cudaMemcpy(this->dev_Vy, this->Vy, this->nfluid_float_alloc, cudaMemcpyHostToDevice));
}

//.. sets neighbors list between fluid-fluid and worm-fluid
void Fluid::ResetNeighborsList(){
	CheckSuccess(cudaMemset((void**)this->dev_nlist, -1, this->parameters->_NMAX * this->nfluid_int_alloc));
	SetFluidNList_N2Kernel <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >>> (this->dev_X, this->dev_Y, this->dev_nlist);
	ErrorHandler(cudaGetLastError());
}

void Fluid::ResetNeighborsList(int itime){
	if (itime % this->parameters->_LISTSETGAP == 0)
		this->ResetNeighborsList();
}

void Fluid::DisplayNList(){
	int * nlist = new int[this->parameters->_NMAX * this->parameters->_NFLUID];
	cudaMemcpy(nlist, this->dev_nlist, this->parameters->_NMAX * this->nfluid_int_alloc,cudaMemcpyDeviceToHost);

	for (int p = 0; p < this->parameters->_NFLUID; p++){
		std::cout << std::endl << p << ": ";
		for (int n = 0; n < this->parameters->_NMAX; n++){
			int id = p * this->parameters->_NMAX + n;
			std::cout << nlist[id] << " ";
		}
	}
	delete[] nlist;
}

void Fluid::DisplayErrors(){
	if (this->errorState.size() > 0){
		printf("\nFLUID ERRORS:\n");
		for (int i = 0; i < this->errorState.size(); i++)
			printf("%i -- %s\n",i,cudaGetErrorString(this->errorState[i]));
		this->errorState.empty();
	}
}
// ------------------------------------------------------------------------------------------
//	PRIVATE METHODS
// ------------------------------------------------------------------------------------------
void Fluid::AllocateHostMemory(){
	this->X = new float[this->parameters->_NFLUID];
	this->Y = new float[this->parameters->_NFLUID];
	this->Vx = new float[this->parameters->_NFLUID];
	this->Vy = new float[this->parameters->_NFLUID];
}

void Fluid::FreeHostMemory(){
	delete[] this->X, this->Y, this->Vx, this->Vy;
}

void Fluid::AllocateGPUMemory(){
	CheckSuccess(cudaMalloc((void**)&(this->dev_X), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Y), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Vx), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Vy), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Fx), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Fy), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Fx_old), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_Fy_old), this->nfluid_float_alloc));
	CheckSuccess(cudaMalloc((void**)&(this->dev_nlist), this->parameters->_NMAX * this->nfluid_int_alloc));
}

void Fluid::FreeGPUMemory(){
	CheckSuccess(cudaFree(this->dev_X));
	CheckSuccess(cudaFree(this->dev_Y));
	CheckSuccess(cudaFree(this->dev_Vx));
	CheckSuccess(cudaFree(this->dev_Vy));
	CheckSuccess(cudaFree(this->dev_Fx));
	CheckSuccess(cudaFree(this->dev_Fy));
	CheckSuccess(cudaFree(this->dev_Fx_old));
	CheckSuccess(cudaFree(this->dev_Fy_old));
	CheckSuccess(cudaFree(this->dev_nlist));
}

void Fluid::DistributeFluidOnHost(){
	for (int j = 0; j < this->parameters->_YDIM; j++){
		for (int i = 0; i < this->parameters->_XDIM; i++){
			float x = float(i) / float(this->parameters->_XDIM);
			float y = float(j) / float(this->parameters->_YDIM);

			int id = j * this->parameters->_XDIM + i;
			this->X[id] = 0.001 + x * this->envirn->_XBOX;
			this->Y[id] = 0.001 + y * this->envirn->_YBOX;
		}
	}
}

void Fluid::ZeroHost(){
	for (int i = 0; i < this->parameters->_NFLUID; i++)
		this->X[i] = this->Y[i] = this->Vx[i] = this->Vy[i] = 0.0f;
}

void Fluid::ZeroGPU(){
	CheckSuccess(cudaMemset((void**)this->dev_X, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Y, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Vx, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Vy, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Fx, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Fy, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Fx_old, 0, this->nfluid_float_alloc));
	CheckSuccess(cudaMemset((void**)this->dev_Fy_old, 0, this->nfluid_float_alloc));
}

void Fluid::CheckSuccess(cudaError_t err){
	if (err != cudaSuccess) this->errorState.push_back(err);
}

void Fluid::CheckSuccess(curandStatus_t stat){
	if (stat == CURAND_STATUS_SUCCESS) return;
	else if (stat == CURAND_STATUS_ALLOCATION_FAILED)
		printf("\nAllocation of random numbers failed ");
	else if (stat == CURAND_STATUS_ARCH_MISMATCH)
		printf("\nArchitecture mismatch for random numbers ");
	else if (stat == CURAND_STATUS_DOUBLE_PRECISION_REQUIRED)
		printf("\nDouble precision requied for random numbers ");
	else if (stat == CURAND_STATUS_INITIALIZATION_FAILED)
		printf("\nInitialization of random numbers failed ");
	else if (stat == CURAND_STATUS_INTERNAL_ERROR)
		printf("\nInternal error generating random numbers ");
	else if (stat == CURAND_STATUS_LAUNCH_FAILURE)
		printf("\nLaunch failure generating random numbers ");
	else if (stat == CURAND_STATUS_LENGTH_NOT_MULTIPLE)
		printf("\nIncorrect length requested for random numbers ");
	else if (stat == CURAND_STATUS_NOT_INITIALIZED)
		printf("\nContainer not initialized before generation ");
	else if (stat == CURAND_STATUS_OUT_OF_RANGE)
		printf("\nOut of range error generating random numbers ");
	else if (stat == CURAND_STATUS_PREEXISTING_FAILURE)
		printf("\nPreexisting error caused failure generating random numbers ");
	else if (stat == CURAND_STATUS_TYPE_ERROR)
		printf("\nType error generating random numbers ");
	else if (stat == CURAND_STATUS_VERSION_MISMATCH)
		printf("\nVersion mismatch encountered generating random numbers ");
	printf("'Fluid' instance.\n");
}

void Fluid::FigureBlockThreadStructure(int tpb){
	if (tpb >= 0) this->Threads_Per_Block = tpb;
	else this->Threads_Per_Block = 256; // good default value
	this->Blocks_Per_Kernel = this->parameters->_NFLUID / this->Threads_Per_Block + 1; // add one to guarentee enough
	this->nfluid_float_alloc = this->parameters->_NFLUID * sizeof(float);
	this->nfluid_int_alloc = this->parameters->_NFLUID * sizeof(int);
}

#endif