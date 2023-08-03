// CUDA runtime
#include <cuda_runtime.h>

// CUDA utilities and system includes
#include <helper_functions.h>
#include <helper_cuda.h>    // includes cuda.h and cuda_runtime_api.h

// Includes
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <string>

// Our application import
#include "particleSystem.h"

#define MAX_EPSILON_ERROR 5.00f
#define THRESHOLD         0.30f

// Simulation parameters
std::string outFile = "nw_more_activity_shorter.xyzv";
uint numFilaments	= 2000;
uint filamentSize 	= 25;
float timestep 		= 0.002f;
int printRate 		= 2500;
int iterations 		= 2500000;
float gravity 		= 0.0f;
float kBend			= 200.f;
float kBond			= 57.146436f;
float hardness		= kBond;
float activity		= 0.05f;
float reverse   	= 0.0005f;
float kbT			= 2.0f;
float drag			= 0.01f;
BoundaryType boundX = BoundaryType::PERIODIC;
BoundaryType boundY = BoundaryType::WALL_NO_SLIP;
BoundaryType boundZ = BoundaryType::WALL;
uint3 gridSize 		= make_uint3(512, 128, 1);
//uint3 gridSize 		= make_uint3(1024, 128, 1);

bool srdSolvent		= true;

extern "C" void cudaInit(int argc, char **argv);


void printInfo(StopWatchInterface *timer, int numParticles, int iter)
{
	// Print some speed metrics
	float fAvgSeconds = ((float)1.0e-3 * (float)sdkGetTimerValue(&timer)/(float)iterations);
	printf("[%i / %i] Throughput = %.4f KParticles/s, Time = %.5f s, Size = %u particles, NumDevsUsed = %u, Workgroup = %u\n",
           iter, iterations, (1.0e-3 * numParticles)/fAvgSeconds, fAvgSeconds, numParticles, 1, 0);
}


void readParameters(int argc, char *argv[])
{
	if (checkCmdLineFlag(argc, (const char **) argv, "iterations"))
        iterations = getCmdLineArgumentInt(argc, (const char **) argv, "iterations");

	if (checkCmdLineFlag(argc, (const char **) argv, "numFilaments"))
        numFilaments = (uint)getCmdLineArgumentInt(argc, (const char **) argv, "numFilaments");

	if (checkCmdLineFlag(argc, (const char **) argv, "filamentSize"))
        filamentSize = (uint)getCmdLineArgumentInt(argc, (const char **) argv, "filamentSize");

	if (checkCmdLineFlag(argc, (const char **) argv, "dt"))
        timestep = getCmdLineArgumentFloat(argc, (const char **) argv, "dt");

	if (checkCmdLineFlag(argc, (const char **) argv, "printRate"))
        printRate = getCmdLineArgumentFloat(argc, (const char **) argv, "printRate");

	if (checkCmdLineFlag(argc, (const char **) argv, "gamma"))
        drag = getCmdLineArgumentFloat(argc, (const char **) argv, "gamma");

	if (checkCmdLineFlag(argc, (const char **) argv, "kbT"))
        kbT = getCmdLineArgumentFloat(argc, (const char **) argv, "kbT");

	if (checkCmdLineFlag(argc, (const char **) argv, "kBend"))
        kBend = getCmdLineArgumentFloat(argc, (const char **) argv, "kBend");

	if (checkCmdLineFlag(argc, (const char **) argv, "kBond"))
        kBond = getCmdLineArgumentFloat(argc, (const char **) argv, "kBond");

	if (checkCmdLineFlag(argc, (const char **) argv, "activity"))
        activity = getCmdLineArgumentFloat(argc, (const char **) argv, "activity");

	if (checkCmdLineFlag(argc, (const char **) argv, "reverse"))
        reverse = getCmdLineArgumentFloat(argc, (const char **) argv, "reverse");

	if (checkCmdLineFlag(argc, (const char **) argv, "hardness"))
        hardness = getCmdLineArgumentFloat(argc, (const char **) argv, "hardness");

	if (checkCmdLineFlag(argc, (const char **) argv, "gridX"))
        gridSize.x = (uint)getCmdLineArgumentInt(argc, (const char **) argv, "gridX");

	if (checkCmdLineFlag(argc, (const char **) argv, "gridY"))
        gridSize.y = (uint)getCmdLineArgumentInt(argc, (const char **) argv, "gridY");

	if (checkCmdLineFlag(argc, (const char **) argv, "gridZ"))
        gridSize.z = (uint)getCmdLineArgumentInt(argc, (const char **) argv, "gridZ");

	if (checkCmdLineFlag(argc, (const char **) argv, "boundX"))
        boundX = static_cast<BoundaryType>(getCmdLineArgumentInt(argc, (const char **) argv, "boundX"));

	if (checkCmdLineFlag(argc, (const char **) argv, "boundY"))
        boundY = static_cast<BoundaryType>(getCmdLineArgumentInt(argc, (const char **) argv, "boundY"));

	if (checkCmdLineFlag(argc, (const char **) argv, "boundZ"))
        boundZ = static_cast<BoundaryType>(getCmdLineArgumentInt(argc, (const char **) argv, "boundZ"));

	if (checkCmdLineFlag(argc, (const char **) argv, "srd"))
		srdSolvent = true;

	if (checkCmdLineFlag(argc, (const char **) argv, "o"))
	{
		char* c = &*outFile.begin();
        getCmdLineArgumentString(argc, (const char **) argv, "o", &c);
		outFile = std::string(c);
	}
}


int main(int argc, char *argv[])
{
	// Initialize CUDA
	cudaInit(argc, argv);

	// Read cmdline args
	readParameters(argc, argv);

	// Instantiate the system of particles
	auto psystem = new ParticleSystem(numFilaments, filamentSize, gridSize, srdSolvent);
	psystem->setBoundaries(boundX, boundY, boundZ);
	psystem->setDamping(drag);
	psystem->setTemperature(kbT);
	psystem->setBondBendingConstant(kBend);
	psystem->setBondSpringConstant(kBond);
	psystem->setActivity(activity);
	psystem->setReverseProbability(reverse);
	psystem->setCollideSpring(hardness);

	psystem->reset();
	psystem->writeOutputs_xyzv(outFile, 0, timestep);

	// Pre time loop actions
	StopWatchInterface *timer = NULL;
	sdkCreateTimer(&timer);
	cudaDeviceSynchronize();
    sdkStartTimer(&timer);

	// Run iterations
	for (int i = 1; i <= iterations; ++i)
    {
        psystem->update(timestep);

		if (i % printRate == 0)
		{
			psystem->writeOutputs_xyzv(outFile, i, timestep);
			// psystem->dumpSolventGrid();
			printInfo(timer, psystem->getNumParticles(), i);
		}
    }

	// Post time loop actions
	cudaDeviceSynchronize();
    sdkStopTimer(&timer);

	delete psystem;
	return 0;
}
