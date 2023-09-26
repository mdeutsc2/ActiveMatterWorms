
#ifndef __FLUID_PARAMETERS_H__
#define __FLUID_PARAMETERS_H__

#include "cuda.h"
#include "cuda_runtime.h"
/* ------------------------------------------------------------------------
*	Data structure containing the parameters of a Fluid object.  Intended
*	to exist at same scope level and Fluid object because Fluid takes a
*	reference at initialization.
--------------------------------------------------------------------------*/
typedef struct {
	
	//.. default setup config
	int	_XDIM, _YDIM;

	//.. # of particles
	int	_NFLUID;

	//.. time steps between resetting neighbors list, max # of neighbors 
	int	_LISTSETGAP, _NMAX;

	//.. LJ energy and length scales
	float	_EPSILON, _SIGMA;

	//.. Langevin thermostat specs
	float	_KBT, _GAMMA;

	//.. LJ pre-calculations
	float	_SIGMA6, _2SIGMA6;
	float	_LJ_AMP;
	float	_RMIN, _R2MIN;
	float	_RCUT, _R2CUT;

	//.. buffer length setting neighbors lists
	float _BUFFER;

	//.. scales interaction w/ landscape
	float _LANDSCALE;

} FluidParameters;
/* ------------------------------------------------------------------------
*	This is the actual symbol that the parameters will be saved to
*	on the GPU device.  Accompanying this global parameter is the
*	function to allocate the parameters to this symbol.
--------------------------------------------------------------------------*/
__constant__ FluidParameters dev_flParams;
/* -----------------------------------------------------------------------
*	Function to allocate FluidParameters on GPU device.  Returns the
*	errors directly.  Basically acts as a typedef for cudaMemcpy(..)
--------------------------------------------------------------------------*/
cudaError_t ParametersToDevice(FluidParameters &params, bool attractivePotentialCutoff = false){

	//.. pre-calculated variables
	params._SIGMA6 = powf(params._SIGMA, 6.0f);
	params._2SIGMA6 = 2.0f * params._SIGMA6;
	params._RMIN = powf(2.0f, 1.0f / 6.0f) * params._SIGMA;
	params._R2MIN = params._RMIN * params._RMIN;
	if (attractivePotentialCutoff) params._RCUT = 2.5f * params._SIGMA;
	else params._RCUT = params._RMIN;
	params._R2MIN = params._RMIN * params._RMIN;
	params._LJ_AMP = 24.0f * params._EPSILON * params._SIGMA6;

	return cudaMemcpyToSymbol(dev_flParams, &params, sizeof(FluidParameters));
}
/*------------------------------------------------------------------------
*	Default values for all parameter values in WormsParameters.
--------------------------------------------------------------------------*/
namespace DEFAULT {
	namespace FLUID {
		static const int	XDIM		= 80;
		static const int	YDIM		= 80;
		static const int	NPARTICLES	= XDIM * YDIM;
		static const int	LISTSETGAP	= 50;
		static const int	NMAX		= 128;
		static const float	EPSILON		= 1.0f;
		static const float	SIGMA		= 1.0f;
		static const float	KBT			= 0.05f;
		static const float	GAMMA		= 2.0f;
		static const float	DAMP		= 3.0f;
		static const float	BUFFER		= 0.25f;
		static const float	LANDSCALE	= 1.0f;
	}
}
//--------------------------------------------------------------------------
#endif