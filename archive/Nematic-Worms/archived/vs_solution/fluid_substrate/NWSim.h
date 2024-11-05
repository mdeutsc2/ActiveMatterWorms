/*******************************
Active Flexible Worms Simulation
********************************/
#ifndef _NWSIM_H
#define _NWSIM_H

#include "NWFluid.h"
#include "NWWorms.h"
#include "NWParams.h"
#include "NWForceXchanger.h"
#include "NWRandomNumberGenerator.h"
//#include "XYZPrinter.h"
#include "NWHostFunctions.h"

using namespace std;

/*
*	Does not work becuase the amount of memory required to do all the random numbers
*	causes errors.  Need to find a way to allocate random numbers segment by segment 
*	and give them to any objects that need them.  Possibly make an object to be referenced
*	on creation of Worms, Fluid, and ForceXchanger that hands off a device ptr when random
*	numbers are needed.  Also, the size of xlist only needs to be as big as _NPARTICLES,
*	do not need to store mostly -1's for all fluid particles.
*/
int NematicWormsSimulation(int argc, char *argv[])
{
	//.. get the shit right
	if (ProcessCommandLine(argc, argv) != 0)
		return 10;

	//.. parameters for worms
	WormsParameters * wParams = new WormsParameters;
	wParams->_XDIM = 5;
	wParams->_YDIM = 50;
	wParams->_NP = 10;
	wParams->_NWORMS = 250;
	wParams->_NPARTICLES = 2500;
	wParams->_LISTSETGAP = 100;
	wParams->_NMAX = 24;
	wParams->_EPSILON = 1.0f;
	wParams->_SIGMA = 1.0f;
	wParams->_DRIVE = 1.0f;
	wParams->_K1 = 57.146f;
	wParams->_K2 = wParams->_K1 * 10.0f;
	wParams->_K3 = 2.0f * wParams->_K2 / 3.0f;
	wParams->_L1 = 0.8f;
	wParams->_L2 = 1.6f;
	wParams->_L3 = 2.4f;
	wParams->_KBT = 1.0f;
	wParams->_GAMMA = 1.0f;
	wParams->_DAMP = 1.0f;
	wParams->_BUFFER = 1.0f;
	ParametersToDevice(*wParams);

	//.. parameters for simulaton
	SimulationParameters * sParams = new SimulationParameters;
	sParams->_DT = 0.001f;
	sParams->_FRAMERATE = 1000;
	sParams->_FRAMESPERFILE = 150;
	sParams->_NSTEPS = 1000000;
	sParams->_XBOX = 200.0f;
	sParams->_YBOX = 200.0f;
	ParametersToDevice(*sParams);

	//.. fluid substrait parameters
	FluidParameters * fParams = new FluidParameters;
	fParams->_XDIM = 100;
	fParams->_YDIM = 100;
	fParams->_NFLUID = 10000;
	fParams->_LISTSETGAP = 100;
	fParams->_NMAX = 24;
	fParams->_EPSILON = 1.0f;
	fParams->_SIGMA = 1.0f;
	fParams->_KBT = 1.0f;
	fParams->_GAMMA = 5.0f;
	fParams->_BUFFER = 1.0f;
	ParametersToDevice(*fParams);

	//.. save parameters for simulation
	SaveSimulationConfiguration();

	//.. random number generator
	GRNG randnumgen(2 * fParams->_NFLUID, 0.0f, 1.0f);

	//.. the physical stuff
	Worms worms;
	Fluid fluid;

	//.. means of force exchange between sets
	ForceXchanger xchanger(worms, fluid, randnumgen);

	worms.Init(&randnumgen, wParams, sParams);
	fluid.Init(&randnumgen, fParams, sParams);
	xchanger.Init();

	worms.ResetNeighborsList();
	fluid.ResetNeighborsList();
	printf("\n\t-------Init--------");
	DisplayErrors(worms, fluid, xchanger, randnumgen);

	PrintXYZ("test1", 
		worms.X, 
		worms.Y, 
		fluid.X, 
		fluid.Y, 
		xchanger.GetXList(),
		*wParams,
		*sParams,
		*fParams);

	/////////////////////////////////////////////////////////////////
	//////////////////// SIMULATION TIME LOOP ///////////////////////
	/////////////////////////////////////////////////////////////////
	
	//int frame = 0;
	std::clock_t start_t = std::clock();
	std::clock_t intra_t = 0.0;
	std::clock_t inter_t = 0.0;
	
	for (int itime = 0; itime < sParams->_NSTEPS; itime++){

		//.. Error Handling
		//ErrorHandler(cudaGetLastError());
		printf("\n\t-------Begin--------");
		ErrorHandler(cudaGetLastError());
		DisplayErrors(worms, fluid, xchanger, randnumgen);

		//.. start clock
		std::clock_t intra_start_t = std::clock();

		//.. calculate internal forces
		worms.InternalForces();
		fluid.InternalForces();

		printf("\n\t-------Internal--------");
		ErrorHandler(cudaGetLastError());
		DisplayErrors(worms, fluid, xchanger, randnumgen);

		//.. set neighbor lists
		if (itime % wParams->_LISTSETGAP == 0)
		{
			//.. reset neighbor list for N2
			worms.ResetNeighborsList();
			fluid.ResetNeighborsList();
			xchanger.UpdateXList();

			printf("\n\t-------SetNlist--------");
			ErrorHandler(cudaGetLastError());
			DisplayErrors(worms, fluid, xchanger, randnumgen);
		}

		//.. clock midpoint
		std::clock_t midpoint = std::clock();
		intra_t += (midpoint - intra_start_t);
		std::clock_t inter_start_t = midpoint;

		//.. calculate lennard jones forces
		worms.LJForces();
		fluid.LJForces();

		printf("\n\t--------LJ---------");
		ErrorHandler(cudaGetLastError());
		DisplayErrors(worms, fluid, xchanger, randnumgen);

		//.. calculate xchange forces
		xchanger.XchangeForces();

		printf("\n\t-------Xchange--------");
		ErrorHandler(cudaGetLastError());
		DisplayErrors(worms, fluid, xchanger, randnumgen);

		inter_t += (std::clock() - inter_start_t);
		
		//cudaDeviceSynchronize();
		worms.Update();
		fluid.Update();

		printf("\n\t-------Update--------");
		ErrorHandler(cudaGetLastError());
		DisplayErrors(worms, fluid, xchanger, randnumgen);

		// Frame calculations
		if (itime%sParams->_FRAMERATE == 0)
		{
			worms.DataDeviceToHost();
			fluid.DataDeviceToHost();
			xchanger.XListToHost();

			printf("\n\t-------Print--------");
			ErrorHandler(cudaGetLastError());
			DisplayErrors(worms, fluid, xchanger, randnumgen);
			
			PrintXYZ("test1",
				worms.X,
				worms.Y,
				fluid.X,
				fluid.Y,
				xchanger.GetXList(),
				*wParams,
				*sParams,
				*fParams);

			//PrintXYZ(worms.X, worms.Y, fluid.X, fluid.Y, xchanger.GetXList());
			//ClockWorks(frame, start_t, intra_t, inter_t);
		}
	}
	return EXIT_SUCCESS;
}

#endif