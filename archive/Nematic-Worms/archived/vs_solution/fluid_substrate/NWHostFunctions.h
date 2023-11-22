// Host CPU functions
// 5.12.15
// Mike Varga

#ifndef _HOST_FUNCTIONS_H
#define _HOST_FUNCTIONS_H
#include "NWmain.h"
#include "NWWorms.h"
#include "NWFluid.h"
#include "NWParams.h"
#include <vector>
#include <string>
#include <cstring>

using namespace std;

//  *****************************************************
// ***************** SIMULATION CONTROL ******************
//  *****************************************************

//.. script needs to mkdir datafolder, otherwise this will fail
__host__ int ProcessCommandLine(int argc, char *argv[])
{
	//.. Not writen yet
	printf("Processing cmdline arguments...\t");

	return 0; // OpenDataFiles();
}

//.. open data files such as *.xyz

//__host__ int OpenDataFiles(void)
//{
//	int numerr = 0;
//
//	//.. generate file names
//	stringstream ssxyz;
//	ssxyz << key << "_" << xyzfile << "_1" << xyz;
//	string sxyz = ssxyz.str();
//	printf("\nOpening *.xyz file:\t '%s'\n", sxyz.c_str());
//	routput.open(sxyz.c_str());
//
//	if (!routput.is_open()) numerr++;
//	if (numerr != 0)
//		printf("Completed with %d errors.\n", numerr);
//	else
//		printf("Complete.\n");
//
//	return numerr;
//}

//__host__ void ErrorHandler(cudaError_t Status)
//{
//	if (cudaSuccess != Status)
//	{
//		fprintf(stderr, "\nfail msg:  '%s'\n", cudaGetErrorString(Status));
//		//fprintf(stderr, "\npress any key to clear memory...");
//		//cin.get();
//		//ShutDownGPUDevice();
//		//CleanUpHost();
//		abort();
//	}
//}

__host__ void DisplayErrors(Worms &w, Fluid &f, ForceXchanger &x, GRNG &r)
{
	char a = 'a';
	f.DisplayErrors();
	w.DisplayErrors();
	x.DisplayErrors();
	r.DisplayErrors();
	std::cin.get(a);
}

//  *****************************************************************
// ******************* SIMULATION UTILITIES **************************
//  *****************************************************************

//__host__ void MovementBC(float &pos, float L)
//{
//	if (pos > L) pos -= L;
//	if (pos < 0) pos += L;
//}

//__host__ void PBC(float &dr, float L)
//{
//	if (dr > L / 2.0) dr -= L;
//	if (dr < -L / 2.0) dr += L;
//}

//.. memory must be pre-copied from device to host
__host__ void PrintXYZ(string fileName, float *wx, float *wy, float *flx, float *fly, int *xlist, 
	WormsParameters &wparams, SimulationParameters &sparams, FluidParameters &fparams)
{
	static int framesinfile = 0;
	static int filenum = 1;
	static const int ntypes = 4;
	static const char ptype[ntypes] = { 'A', 'B', 'C', 'D'};
	static string key = fileName;

	//.. close then open new file when full
	if (framesinfile >= sparams._FRAMESPERFILE)
	{
		//.. adjust counting and current file
		framesinfile = 0;
		filenum++;
		routput.close();

		//.. make new file name and open
		stringstream ssxyz;
		ssxyz << key << "_" << xyzfile << "_" << filenum << xyz;
		string sxyz = ssxyz.str();
		routput.open(sxyz.c_str());

		//.. announce change
		if (routput.is_open())
			printf("\nChanging output file to %s\n\n", sxyz.c_str());
		else
			printf("\nError openning new file\n\n");
	}

	// Write to file
	framesinfile++;

	//.. top of the frame data
	routput << wparams._NPARTICLES + fparams._NFLUID + 4 << endl;
	routput << wparams._NP << " " << wparams._K2 << " " 
		<< sparams._XBOX << " " << sparams._YBOX << endl;

	//.. print worm particles
	//bool * bound = new bool[_NFLUID];
	for (int i = 0; i < wparams._NPARTICLES; i++)
	{
		//.. exit if blown up
		if (isnan(wx[i]))
		{
			fprintf(stderr, "\nWORM BLEW UP!!!\n\nQUITTING IMMEDIATELY\n");
			cudaDeviceReset();
			abort();
		}

		//.. choose type
		int t = (i / wparams._NP) % ntypes;

		//.. print type and positions
		routput << ptype[t] << " " << wx[i] << " " << wy[i] << " 0 " << endl;
	}

	//.. print fluid particles (shifted in z direction by 1.0)
	for (int i = 0; i < fparams._NFLUID; i++)
	{
		if (isnan(flx[i]))
		{
			fprintf(stderr, "\nFLUID BLEW UP!!!\n\nQUITTING IMMEDIATELY\n");
			cudaDeviceReset();
			abort();
		}

		//.. type based on binding status
		if (xlist[i] >= 0) routput << "E ";
		else routput << "G ";

		//.. print positions of fluid particle
		routput << flx[i] << " " << fly[i] << " -1.0 " << endl;
	}

	//Mark box corners
	routput << "G " << "0 " << "0 " << "0" << endl;
	routput << "E " << sparams._XBOX << " 0 " << "0" << endl;
	routput << "E " << "0 " << sparams._YBOX << " 0" << endl;
	routput << "E " << sparams._XBOX << " " << sparams._YBOX << " 0" << endl;
}

//.. save to *_cfg.dat file parameters used
__host__ void SaveSimulationConfiguration(void)
{
	/*//.. file name and open
	std::ofstream save;
	stringstream ssname;
	ssname << key << "_cfg" << dat;
	string name = ssname.str();
	save.open(name.c_str());

	if (!save.is_open()) return;

	printf("Parameter Configuration Saved to...\n '%s'\n", name.c_str());

	//.. this is the necessary format
	save << "version="		<< _VERSION << endl;
	save << "key="			<< key << endl;
	save << "np="			<< _NP << endl;
	save << "nworms="		<< _NWORMS << endl;
	save << "nsteps="		<< _NSTEPS << endl;
	save << "nparticles="	<< _NPARTICLES << endl;
	save << "dt="			<< _DT << endl;
	save << "mu="			<< _MU << endl;
	save << "xbox="			<< _XBOX << endl;
	save << "ybox="			<< _YBOX << endl;
	//save << "blocksperkernel=" << _BLOCKS_PER_KERNEL << endl;
	//save << "threadsperblock=" << _THREADS_PER_BLOCK << endl;
	save << "epsilon="		<< _EPSILON << endl;
	save << "sigma="		<< _SIGMA << endl;
	save << "drive="		<< _DRIVE << endl;
	//save << "criticaldot="	<< _EXT_CRIT_DOT << endl;
	save << "k1="			<< _K1 << endl;
	save << "k2="			<< _K2 << endl;
	save << "k3="			<< _K3 << endl;
	save << "l1="			<< _L1 << endl;
	save << "l2="			<< _L2 << endl;
	save << "l3="			<< _L3 << endl;
	save << "rcut="			<< _RCUT << endl;
	save << "framerate="	<< _FRAMERATE << endl;
	save << "nlistmax="		<< _NMAX << endl;
	save << "nlistsetgap="	<< _LISTSETGAP << endl;

#ifdef _EXTENSILE
	save << "extensile=1" << endl;
#else
	save << "exensile=0" << endl;
#endif

#ifdef _STICKY
	save << "sticky=1" << endl;
#else
	save << "sticky=0" << endl;
#endif

#ifdef _NOISE
	save << "noise=" << _KBT << endl;
#else
	save << "noise=0" << endl;
#endif

#ifdef _DRAG
	save << "drag=" << _GAMMA << endl;
#else
	save << "drag=0" << endl;
#endif

#ifdef _DAMP
	save << "damp=" << _DAMP << endl;
#else
	save << "damp=0" << endl;
#endif

#ifdef __DEBUG__
	save << "...DEBUGING..." << endl;
#endif

	save.close();*/
}
/*
__host__ void ClockWorks(int &frame, std::clock_t &start_t, std::clock_t &intra_t, std::clock_t &inter_t)
{
	//.. Estimate time remaining
	double frame_t = (std::clock() - start_t) / (double)(CLOCKS_PER_SEC*60.0);
	double frames_to_go = double(_NSTEPS / _FRAMERATE) - (double)(frame - 1);
	double min_to_go = frame_t * frames_to_go;
	double _hours = min_to_go / 60.0;
	int hours = int(_hours);
	int minutes = int((_hours - (double)hours) * 60);
	cout << "--------------------------------------------" << endl;
	cout << "\n ***\tFrame " << frame++ << " saved." << endl;
	cout << "\n ***\tEstimated time until completion...\n\t"
		<< hours << "hrs " << minutes << "min" << endl << endl;
	start_t = std::clock();

	// Write statistics
	cout << " ***\tAverage clocks spent calculating...\n"
		<< "\tIntra-forces: " << intra_t / (double)_FRAMERATE
		<< "\n\tInter-forces: " << inter_t / (double)_FRAMERATE << endl;

	intra_t = 0.0;
	inter_t = 0.0;
}*/


#endif 