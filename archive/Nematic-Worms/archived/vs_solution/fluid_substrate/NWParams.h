// Global Parameters for Simulation
// 5.12.15
// Mike Varga

#ifndef _PARAMS_H
#define _PARAMS_H

#include <math.h>
#include <vector>
#include <fstream>

//.. mathamatical constants
#define PI		3.141592654f
#define PIo2	1.570796327f
#define _216	1.12246205f

const char* _VERSION = "v1.1";
//const std::string key = "BS4";

//.. initial placemenet dimensions (worms is per head)
//#define _XDIM		4
//#define _YDIM		32
//#define _FL_XDIM	64
//#define _FL_YDIM	64

//.. system defining parameters for worms
//#define _NP			10
//#define _NWORMS		70 //(_XDIM * _YDIM)
//#define _NSTEPS		20000000 //53333500 
//#define _NPARTICLES 700 //2560  //(_NWORMS*_NP)
//#define _NFLUID		4096 //(_FL_XDIM * _FL_YDIM)
//#define _DT			0.00075f 
//#define _MU			0.15f
//#define _XBOX		800.0f
//#define _YBOX		800.0f
//#define _LISTSETGAP 80
//#define _FRAMERATE	5000

//.. LJ parameters
//#define _EPSILON	0.5f 
//#define _SIGMA		1.0f
//#define _FL_EPS		1.0f
//#define _FL_SIG		0.5f

//.. landscape parameters
#define _N		8
//#define _KK		(2 * PI * _N / _XBOX)
#define _LD		1.0f

//.. driving mechanism parameters
//#define _DRIVE		 0.5f  
#define _KBIND		 200.0f
#define _MAKEBINDING 1.8f // compared to stddev of gaussian
#define _LOSSBINDING 1.2f
#define _WALKBINDING 0.000f

//.. inter worm spring constants
//#define _K1		57.146f
//#define _K2		(_K1 * 20.0f) 
//#define _K3		(2.0f * _K2/3.0f)
//#define _L1		0.7f
//#define _L2		(2.0f * _L1)
//#define _L3		(3.0f * _L1)

//.. random variations on all particles
#define _NOISE
//#define _KBT		0.5f
//#define _FL_KBT		0.5f

//.. Langevin therostat
#define _DRAG
//#define _GAMMA		2.0f
//#define _FL_GAM		2.0f

//.. dampeners for inter worm springs
#define _DAMPING1
#define _DAMPING2
//#define _DAMP	3.0f
#define _DAMPING

//.. turns on polarity similar attractive interactions
//#define _STICKY
#define _XCUT2		0.25f

//.. interaction pre-calculations (worms)
//#define _SIGMA6		1.0f
//#define _2SIGMA6	(2.0f * _SIGMA6)
//#define _LJ_AMP		(24.0f * _EPSILON * _SIGMA6)
//#define _RMIN		(_216 * _SIGMA)  // only for sigma = 1
//#define _R2MIN		(_RMIN * _RMIN)
//#define _RCUT		(2.50000f * _SIGMA) //_RMIN //
//#define _R2CUT		(_RCUT * _RCUT)
//#define _LENGTH		(_RCUT + (float)(_NP - 1.0f)*_L1)

//.. interaction pre-calculations (fluid)
//#define _FL_SIG6	0.015625f
//#define _FL_2SIG6	(2.0f * _FL_SIG6)
//#define _FL_LJ_AMP	(24.0f * _FL_EPS * _FL_SIG6)
//#define _FL_RMIN	(_216 * _FL_SIG) //1.12246205f // only for sigma = 1
//#define _FL_R2MIN	(_FL_RMIN * _FL_RMIN)
//#define _FL_RCUT	_FL_RMIN
//#define _FL_R2CUT	(_FL_RCUT * _FL_RCUT)


//.. simulation interfacing
//#define _FRAMESPERFILE 120

std::ifstream load;
std::ofstream routput;
const std::string xyzfile = "nwsim";
const std::string xyz = ".xyz";
const std::string xyzv = ".xyzv";
const std::string xyzc = ".xyzc";
const std::string csv = ".csv";
const std::string dat = ".dat";
const std::string datafolder = "./app_data/";

//.. max # allowed in neighbors lists
//#define _NMAX		28
//#define _BUFFER		1.25f

//.. data allocation numbers
//static const size_t nparticles_float_alloc = _NPARTICLES*sizeof(float);
//static const size_t nparticles_int_alloc   = _NPARTICLES*sizeof(int);
//static const size_t nfluid_float_alloc	   = _NFLUID*sizeof(float);
//static const size_t nfluid_int_alloc	   = _NFLUID*sizeof(int);

#endif