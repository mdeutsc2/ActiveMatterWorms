
#ifndef __SIMULATION_PARAMETERS_H__
#define __SIMULATION_PARAMETERS_H__
// 2D
#include "cuda.h"
#include "cuda_runtime.h"
#include <string>
#include <cstdlib>
#include "NWmain.h"
/* ------------------------------------------------------------------------
*	Data structure containing the parameters needed to run a simulation.  
*	Intended to exist inside a simulaltion object, NOT taken by reference
*	at initialization such as WormsParameters and FluidParameters.
----------------------------------------------------------------------------*/
typedef struct {

	//.. # of time steps
	int _NSTEPS;

	//.. rate of printing frames, and max frames per output file
	int _FRAMERATE, _FRAMESPERFILE;

	//.. long time integration constant
	float _DT;

	//.. number of inner time steps for fast potential calculations
	int _NSTEPS_INNER;

	//..system physical size
	float _XBOX, _YBOX, _ZBOX;
	float _BOX[_D_];

	//.. flags for sim prodecures
	bool _LMEM, // linear memory on gpu 
		_PBC[_D_], // periodic boundary conditions
		_SBC[_D_], // soft wall boundary conditions
		_HBC[_D_], // hard wall boundary conditions
		_SPHERE; // circular boundaries 2d, spherical in 3d

	//.. strength of wall potentials
	float _Kw;


} SimulationParameters;
/* ------------------------------------------------------------------------
*	This is the actual symbol that the parameters will be saved to
*	on the GPU device.  Accompanying this global parameter is the
*	function to allocate the parameters to this symbol.
---------------------------------------------------------------------------*/
__constant__ SimulationParameters dev_simParams;
/* -------------------------------------------------------------------------
*	Function to allocate SimulationParameters on GPU device.  Returns 
*	the errors directly.  Basically acts as a typedef for cudaMemcpy(..)
----------------------------------------------------------------------------*/
cudaError_t ParametersToDevice(SimulationParameters &params){
	return cudaMemcpyToSymbol(dev_simParams, &params, sizeof(SimulationParameters));
}
/*--------------------------------------------------------------------------
*	Default values for all parameter values in SimulationParameters.
--------------------------------------------------------------------------*/
namespace DEFAULT {
	namespace SIM {
		static const int NSTEPS = 100000;
		static const int NSTEPS_INNER = 10;
		static const int FRAMERATE = 1000;
		static const int FRAMESPERFILE = 100;
		static const float DT = 0.005f;
		static const float XBOX = 100.0f;
		static const float YBOX = 100.0f;
		static const float ZBOX = 100.0f;
		static const std::string FILENAME = "output.xyz";
		static const bool LMEM = false;
		static const bool PBC[_D_] = { true, true
#if _D_ == 3 
			, true // for 3rd dim
#endif
		};
		static const bool SBC[_D_] = { false, false
#if _D_ == 3
			, true
#endif
		};
		static const float KWALL = 5.0f;
		static const bool HBC[_D_] = { false, false
#if _D_ == 3
			, true
#endif
		};
		static const bool SPHERE = false;
		static const float Kw = 5.0f;
	}
}
//--------------------------------------------------------------------------
void set_unset_D_bool_vectors(const std::string& val, bool set[_D_], bool unset1[_D_], bool unset2[_D_]){
	if (val == "x"){
		set[0] = true;
		unset1[0] = unset2[0] = false;
	}
	else if (val == "y"){
		set[1] = true;
		unset1[1] = unset2[1] = false;
	}
#if _D_ == 3
	else if (val == "z"){
		set[2] = true;
		unset1[2] = unset2[2] = false;
	}
#endif
	else if (val == "xy" || val == "yx"){
		set[0] = set[1] = true;
		unset1[0] = unset2[0] = false;
		unset1[1] = unset2[1] = false;
	}
	else if (val == "xz" || val == "zx") {
		set[0] = true;
		unset1[0] = unset2[0] = false;
#if _D_ == 3
		set[2] = true;
		unset1[2] = unset2[2] = false;
#endif
	}
	else if (val == "yz" || val == "zy") {
		set[1] = true;
		unset1[1] = unset2[1] = false;
#if _D_ == 3
		set[2] = true;
		unset1[2] = unset2[2] = false;
#endif
	}
	else if (val == "xyz") {
		set[0] = set[1] = true;
		unset1[0] = unset2[0] = false;
		unset1[1] = unset2[1] = false;
#if _D_ == 3
		set[2] = true;
		unset1[2] = unset2[2] = false;
#endif
	}
	else printf("\n[ ERROR ] : Unrecognized option argument ('%s')", val.c_str());
}
//--------------------------------------------------------------------------
void GrabParameters(SimulationParameters * parameters, int argc, char *argv[], std::string &outfile){

	//.. cycle through arguments
	for (int i = 1; i < argc; i++){
		std::string arg = argv[i];
		std::string val;
		if (arg == "-dt"){
			if (i + 1 < argc){
				parameters->_DT = std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-xbox"){
			if (i + 1 < argc){
				parameters->_XBOX = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-ybox"){
			if (i + 1 < argc){
				parameters->_YBOX = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-zbox"){
			if (i + 1 < argc){
				parameters->_ZBOX = (int)std::strtof(argv[1 + i++], NULL);
				if (_D_ != 3) 
					printf("\n[ ERROR ] : Can not assign size of 3rd dimension in 2D simulation");
			}
		}
		else if (arg == "-nsteps"){
			if (i + 1 < argc){
				parameters->_NSTEPS = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-nsteps-inner"){
			if (i + 1 < argc){
				parameters->_NSTEPS_INNER = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-framerate"){
			if (i + 1 < argc){
				parameters->_FRAMERATE = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-framesperfile"){
			if (i + 1 < argc){
				parameters->_FRAMESPERFILE = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-o" || arg == "--output"){
			if (i + 1 < argc){
				outfile = std::string(argv[1 + i++]);
			}
		}
		else if (arg == "-lmem"){
			parameters->_LMEM = true;
		}
		else if (arg == "-softwalls"){
			if (i + 1 < argc){
				val = argv[1 + i++];
				set_unset_D_bool_vectors(val, parameters->_SBC, parameters->_HBC, parameters->_PBC);
			}
			else {
				for_D_{
					parameters->_SBC[d] = true;
					parameters->_PBC[d] = false;
					parameters->_HBC[d] = false;
				}
			}
		}
		else if (arg == "-kw"){
			if (i + 1 < argc){
				parameters->_Kw = (int)std::strtof(argv[1 + i++], NULL);
			}
		}
		else if (arg == "-hardwalls"){
			if (i + 1 < argc){
				val = argv[1 + i++];
				set_unset_D_bool_vectors(val, parameters->_HBC, parameters->_SBC, parameters->_PBC);
			}
			else {
				for_D_{
					parameters->_HBC[d] = true;
					parameters->_PBC[d] = false;
					parameters->_SBC[d] = false;
				}
			}
		}
		else if (arg == "-periodic"){
			if (i + 1 < argc){
				val = argv[1 + i++];
				set_unset_D_bool_vectors(val, parameters->_PBC, parameters->_HBC, parameters->_SBC);
			}
			else {
				for_D_{
					parameters->_PBC[d] = true;
					parameters->_SBC[d] = false;
					parameters->_HBC[d] = false;
				}
			}
		}
		else if (arg == "-spherical"){
			parameters->_SPHERE = true;
		}
	}
}
//--------------------------------------------------------------------------
void Init(SimulationParameters * parameters, int argc, char *argv[], std::string &outfile){

	//.. init with default parameters
	parameters->_DT = DEFAULT::SIM::DT;
	parameters->_XBOX = DEFAULT::SIM::XBOX;
	parameters->_YBOX = DEFAULT::SIM::YBOX;
	parameters->_ZBOX = DEFAULT::SIM::ZBOX;
	parameters->_NSTEPS = DEFAULT::SIM::NSTEPS;
	parameters->_NSTEPS_INNER = DEFAULT::SIM::NSTEPS_INNER;
	parameters->_FRAMERATE = DEFAULT::SIM::FRAMERATE;
	parameters->_FRAMESPERFILE = DEFAULT::SIM::FRAMESPERFILE;
	outfile = DEFAULT::SIM::FILENAME;
	parameters->_LMEM = DEFAULT::SIM::LMEM;
	for_D_ parameters->_PBC[d] = DEFAULT::SIM::PBC[d];
	for_D_ parameters->_HBC[d] = DEFAULT::SIM::HBC[d];
	for_D_ parameters->_SBC[d] = DEFAULT::SIM::SBC[d];
	parameters->_SPHERE = DEFAULT::SIM::SPHERE;
	parameters->_Kw = DEFAULT::SIM::Kw;

	//.. get assign cmdline parameters
	GrabParameters(parameters, argc, argv, outfile);

	//.. calculate parameters
	const float box[3] = { parameters->_XBOX, parameters->_YBOX, parameters->_ZBOX };
	for_D_ parameters->_BOX[d] = box[d];

	//.. put on GPU and check for error
	cudaError_t err;
	err = ParametersToDevice(*parameters);
	std::cout << "\nSimulation parameters cudaMemcpyToSymbol returned:     \t" << cudaGetErrorString(err);
}
//--------------------------------------------------------------------------
#endif