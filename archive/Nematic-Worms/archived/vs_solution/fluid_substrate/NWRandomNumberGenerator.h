
#ifndef __RANDOM_NUMBER_GENERATOR_H__
#define __RANDOM_NUMBER_GENERATOR_H__
// ----------------------------------------------------------------------------
//	INIT and GENERATE KERNELS (when using curandState methods)
// ----------------------------------------------------------------------------
__global__ void initKernel(curandState *dev_state, int max){
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < max){
		curand_init(1337, id, 0, &dev_state[id]);
	}
}
// ----------------------------------------------------------------------------
__global__ void generateKernel_2d(curandState *dev_state, float *dev_rand, int number){
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int idy = threadIdx.y + blockDim.y * blockIdx.y;
	int width = blockDim.x * gridDim.x;
	int tid = idy*width + idx;
	if (tid < number){
		dev_rand[tid] = curand_normal(&dev_state[tid]);
	}
}
// ----------------------------------------------------------------------------
__global__ void generateKernel_1d(curandState *dev_state, float *dev_rand, int number){
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < number){
		dev_rand[idx] = curand_normal(&dev_state[idx]);
	}
}
// ----------------------------------------------------------------------------
#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
#include "curand_kernel.h"
#include "NWmain.h"
/* ----------------------------------------------------------------------------
*	brief: Gaussian Random Number Generator
*
*	Allows calls to one static random number generator.  Has one fixed size
*	ptr on the device that fills when calls and returns that ptr.  Used in
*	other classes when random numbers are needed.  This way there are not lots
*	of different number generators and device ptrs full of random numbers.
------------------------------------------------------------------------------*/
class GRNG {
	//.. device ptr for numbers
	float * dev_n;
	//curandState *dev_s;

	//.. actual number generator
	curandGenerator_t * generator;

	//.. status of the generator
	std::vector<curandStatus_t> status;
	std::vector<cudaError_t> errorState;

	//.. size of allocation on device
	const size_t alloc_size;
	int max_size;
	int Threads_Per_Block;
	int Blocks_Per_Kernel;

	//.. Gaussian distribution specs
	const float mean;
	const float stddev;

public:
	GRNG(int maxCallSize, float distributionMean, float standardDeviation);
	~GRNG();

	float* Get(unsigned count, bool uniform);
	float* Generate(unsigned count);
	float* GenerateAll();
	void DisplayErrors();

private:
	void AllocateGPUMemory(int N);
	void FreeGPUMemory();
	void CheckStatus(curandStatus_t stat);
	void CheckSuccess(cudaError_t err);
};
// ---------------------------------------------------------------------------
//	PUBLIC METHODS
// ---------------------------------------------------------------------------
GRNG::GRNG(int maxCallSize, float distributionMean, float standardDeviation)
	: generator(new curandGenerator_t), dev_n(NULL),
	alloc_size(sizeof(float)*maxCallSize), mean(distributionMean),
	stddev(standardDeviation), max_size(maxCallSize){
	
	CheckStatus(curandCreateGenerator(this->generator, CURAND_RNG_PSEUDO_DEFAULT));
	
	this->Threads_Per_Block = 256;
	this->Blocks_Per_Kernel = (this->max_size / this->Threads_Per_Block) + 1;
	this->AllocateGPUMemory(this->max_size);
	//initKernel <<<this->Blocks_Per_Kernel, this->Threads_Per_Block >>> (this->dev_s, this->max_size);
	printf("\nRandom Number Generator (GRNG) Created");
}
//-------------------------------------------------------------------------------------------
GRNG::~GRNG(){
	this->FreeGPUMemory();
	CheckStatus(curandDestroyGenerator(*this->generator));
}
//-------------------------------------------------------------------------------------------
float* GRNG::Get(unsigned count, bool uniform = false){
	size_t call_size = count * sizeof(float);
	if (call_size > this->alloc_size){
		printf("\n\tWarning:\ttoo large of call to RNG!\n");
		count = this->alloc_size / sizeof(float);
	}
	if (uniform)
		CheckStatus(curandGenerateUniform(*this->generator, this->dev_n, count));
	else
		CheckStatus(curandGenerateNormal(*this->generator, this->dev_n, count, this->mean, this->stddev));
	CheckSuccess(cudaDeviceSynchronize());
	return this->dev_n;
}
//-------------------------------------------------------------------------------------------
float* GRNG::Generate(unsigned count){
	
	//.. protect against out of range
	if (count > this->max_size){
		printf("\n\tWarning:\ttoo large of call to RNG!\n");
		count = this->max_size;
	}

	//.. make grid and block structure
	int c_root = (int)ceil(sqrt(float(count)));
	int t_root = 32;
	dim3 blockStruct(t_root, t_root);
	dim3 gridStruct(c_root / t_root + 1, c_root / t_root + 1);

	//.. generate random numbers
	/*generateKernel_2d <<< blockStruct, gridStruct >>>
	(
		this->dev_s,
		this->dev_n,
		count
	);*/

	//.. check for errors
	CheckSuccess(cudaDeviceSynchronize());
	CheckSuccess(cudaGetLastError());

	//.. return device pointer
	return this->dev_n;
}
//-------------------------------------------------------------------------------------------
float* GRNG::GenerateAll(){

	//.. fill entire array with random numbers
	/*generateKernel_1d <<< this->Blocks_Per_Kernel, this->Threads_Per_Block >>> 
	(
		this->dev_s, 
		this->dev_n, 
		this->max_size
	);*/

	/*int c_root = (int)ceil(sqrt(float(this->max_size)));
	int t_root = 32;
	dim3 blockStruct(t_root, t_root);
	dim3 gridStruct(c_root / t_root + 1, c_root / t_root + 1);

	generateKernel_2d <<< blockStruct, gridStruct >>>
	(
		this->dev_s,
		this->dev_n,
		this->max_size
	);*/

	//.. check for errors
	CheckSuccess(cudaDeviceSynchronize());
	CheckSuccess(cudaGetLastError());

	//.. return device pointer
	return this->dev_n;
}
//-------------------------------------------------------------------------------------------
void GRNG::DisplayErrors(){
	if (this->status.size() > 0 || this->errorState.size() > 0){
		printf("\nGRNG ERRORS:\n");
		for (int i = 0; i < this->status.size(); i++){
			curandStatus_t stat = this->status[i];
			if (stat == CURAND_STATUS_ALLOCATION_FAILED)
				printf("Allocation of random numbers failed\n");
			else if (stat == CURAND_STATUS_ARCH_MISMATCH)
				printf("Architecture mismatch for random numbers\n");
			else if (stat == CURAND_STATUS_DOUBLE_PRECISION_REQUIRED)
				printf("Double precision requied for random numbers\n");
			else if (stat == CURAND_STATUS_INITIALIZATION_FAILED)
				printf("Initialization of random numbers failed\n");
			else if (stat == CURAND_STATUS_INTERNAL_ERROR)
				printf("Internal error generating random numbers\n");
			else if (stat == CURAND_STATUS_LAUNCH_FAILURE)
				printf("Launch failure generating random numbers\n");
			else if (stat == CURAND_STATUS_LENGTH_NOT_MULTIPLE)
				printf("Incorrect length requested for random numbers\n");
			else if (stat == CURAND_STATUS_NOT_INITIALIZED)
				printf("Container not initialized before generation\n");
			else if (stat == CURAND_STATUS_OUT_OF_RANGE)
				printf("Out of range error generating random numbers\n");
			else if (stat == CURAND_STATUS_PREEXISTING_FAILURE)
				printf("Preexisting error caused failure generating random numbers\n");
			else if (stat == CURAND_STATUS_TYPE_ERROR)
				printf("Type error generating random numbers\n");
			else if (stat == CURAND_STATUS_VERSION_MISMATCH)
				printf("Version mismatch encountered generating random numbers\n");
		}
		for (int i = 0; i < this->errorState.size(); i++)
			printf("%i -- %s\n", i, cudaGetErrorString(this->errorState[i]));

		this->errorState.empty();
		this->status.empty();
	}
}
// -------------------------------------------------------------------------
//	PRIVATE METHODS
// -------------------------------------------------------------------------
void GRNG::AllocateGPUMemory(int N){
	CheckSuccess(cudaMalloc((void**)&this->dev_n, this->alloc_size));
	//CheckSuccess(cudaMalloc((void**)&this->dev_s, N*sizeof(curandState)));
}
//-------------------------------------------------------------------------------------------
void GRNG::FreeGPUMemory(){
	CheckSuccess(cudaFree(this->dev_n));
	//CheckSuccess(cudaFree(this->dev_s));
}
//-------------------------------------------------------------------------------------------
void GRNG::CheckStatus(curandStatus_t stat){
	if (stat != CURAND_STATUS_SUCCESS) this->status.push_back(stat);
}
//-------------------------------------------------------------------------------------------
void GRNG::CheckSuccess(cudaError_t err){
	if (err != cudaSuccess) this->errorState.push_back(err);
}
//-------------------------------------------------------------------------------------------
#endif