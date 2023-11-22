
#ifndef __HALLOWEEN_SIMULATION_H__
#define __HALLOWEEN_SIMULATION_H__

#include "NWmain.h"
#include "NWWorms.h"
#include "NWRandomNumberGenerator.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
#include "NWKernels.h"
// --------------------------------------------------------------------------------------------
//	Class That Runs An Active Worms Simulations Spelling "Happy Halloween"
// --------------------------------------------------------------------------------------------
class HalloweenSimulation {
	
	//.. particles and RNG
	Worms				 * worms;
	WormsParameters		 * params;
	GRNG				 * rng;
	SimulationParameters * simparams;

	//.. file stream
	std::ofstream fxyz;

public:
	HalloweenSimulation();
	~HalloweenSimulation();
	void Run();

private:
	void PlaceWorms();
	void XYZPrint(int itime);
	void DisplayErrors();

	//.. for setting up the letters
	float _S; 		// gap between letters
	float _LW;		// width of letter
	float _LMH;		// verticle height that horizontal worms can start
	float _LMV;		// verticle height that verticel worms can start
	float _LT;		// verticle height that top worms can start
	float _2ND;		// y0 height of the "HAPPY" line
};
// ------------------------------------------------------------------------------------------
//	PUBLIC METHODS
// ------------------------------------------------------------------------------------------
HalloweenSimulation::HalloweenSimulation(){

	//.. reset device before starting
	cudaDeviceReset();

	//.. setup parameters (should be done with user input)
	this->params = new WormsParameters();
	this->params->_XDIM		  = 11;
	this->params->_YDIM		  = 10;
	this->params->_NP		  = 10;
	this->params->_NWORMS	  = params->_XDIM * params->_YDIM;
	this->params->_NPARTICLES = params->_NP * params->_NWORMS;
	this->params->_LISTSETGAP = 200;
	this->params->_NMAX		  = 28;
	this->params->_EPSILON	  = 1.0f;
	this->params->_SIGMA	  = 1.0f;
	this->params->_DRIVE	  = 0.5f;
	this->params->_K1		  = 57.146f;
	this->params->_K2		  = 10.0f * params->_K1;
	this->params->_K3		  = 2.0f * params->_K2 / 3.0f;
	this->params->_L1		  = 0.8f;
	this->params->_L2		  = 2.0f * params->_L1;
	this->params->_L3		  = 3.0f * params->_L1;
	this->params->_KBT		  = 0.5f;
	this->params->_GAMMA	  = 2.0f;
	this->params->_DAMP		  = 3.0f;
	this->params->_BUFFER	  = 1.5f;
	this->params->_SIGMA6	  = powf(params->_SIGMA, 6.0f);
	this->params->_2SIGMA6	  = params->_SIGMA6 * 2.0f;
	this->params->_LJ_AMP	  = 24.0f * params->_EPSILON * params->_SIGMA6;
	this->params->_RMIN		  = _216 * params->_SIGMA;
	this->params->_R2MIN	  = params->_RMIN * params->_RMIN;
	this->params->_RCUT		  = params->_RMIN; // 2.5f * params->_SIGMA;
	this->params->_R2CUT	  = params->_RCUT * params->_RCUT;

	//.. parameters to device
	cudaError_t err;
	err = ParametersToDevice(*params);
	std::cout << "Worms parameters cudaMemcpyToSymbol returned: " << cudaGetErrorString(err) << std::endl;

	this->simparams = new SimulationParameters();
	this->simparams->_DT			= 0.001f;
	this->simparams->_FRAMERATE		= 4000;
	this->simparams->_FRAMESPERFILE = 200;
	this->simparams->_NSTEPS		= 1000000;
	this->simparams->_XBOX			= 110.0f;
	this->simparams->_YBOX			= 120.0f;

	//.. parameters to device
	err = ParametersToDevice(*simparams);
	std::cout << "Simulation parameters cudaMemcpyToSymbol returned: " << cudaGetErrorString(err) << std::endl;

	CheckParametersOnDevice<<< 1, 1 >>>();
	cudaDeviceSynchronize();

	//.. setup random number generator
	this->rng = new GRNG(2 * params->_NPARTICLES, 0.0f, 1.0f);

	//.. define the worms
	this->worms = new Worms();

	//.. open output xyz file for the simulation
	this->fxyz.open("output//memory_test1.xyz");

	//.. assign letter sizing and characteristics
	this->_S = 3.0f;
	this->_LW = ((params->_NP - 1) * params->_L1 + params->_RMIN);
	this->_LMH = _LW;
	this->_LMV = _LMH + params->_RMIN;
	this->_LT = 2.0f * _LMV - params->_RMIN;
	this->_2ND = 2.0f * _S + _LT;
}

HalloweenSimulation::~HalloweenSimulation(){
	delete this->worms;
	delete this->rng;
	delete this->params;
	delete this->simparams;
	this->fxyz.close();
}

/*
*	This initiates and runs the simulation.  All must be set first
*	in order for this to run properly.
*/
void HalloweenSimulation::Run(){
	//this->PlaceWorms();
	this->worms->Init(this->rng, this->params, this->simparams);
	this->DisplayErrors();
	for (int itime = 0; itime < this->simparams->_NSTEPS; itime++){
		this->worms->InternalForces();
		this->worms->ResetNeighborsList(itime);
		this->worms->LJForces();
		this->worms->AutoDriveForces();
		this->worms->Update();
		this->XYZPrint(itime);
	}
}
// ------------------------------------------------------------------------------------------
//	PRIVATE METHODS
// ------------------------------------------------------------------------------------------
void HalloweenSimulation::PlaceWorms(){

	float *x = new float[params->_NWORMS];
	float *y = new float[params->_NWORMS];
	float *t = new float[params->_NWORMS];
	float x0, y0;

	////////////////// SPELLING HALLOWEEN ///////////////////////

	//.. H of Halloween
	x[0] = _S;			y[0] = _S;			t[0] = PIo2;
	x[1] = _S;			y[1] = _S + _LMH;	t[1] = 0.0f;
	x[2] = _S;			y[2] = _S + _LMV;	t[2] = PIo2;
	x[3] = _S + _LW;	y[3] = _S;			t[3] = PIo2;
	x[4] = _S + _LW;	y[4] = _S + _LMV;	t[4] = PIo2;

	//.. A fo Halloween
	x0 = 2 * _S + _LW;
	y0 = _S;
	x[5] = x0;						y[5] = y0;			t[5] = PIo2;
	x[6] = x0;						y[6] = y0 + _LMH;	t[6] = 0.0f;
	x[7] = x0;						y[7] = y0 + _LMV;	t[7] = PIo2;
	x[8] = x0 + params->_RMIN/2;	y[8] = y0 + _LT;	t[8] = 0.0f;
	x[9] = x0 + _LW;				y[9] = y0;			t[9] = PIo2;
	x[10] = x0 + _LW;				y[10] = y0 + _LMV;	t[10] = PIo2;

	//.. 1st L of Halloween
	x0 = 3 * _S + 2 * _LW;
	y0 = _S;
	x[11] = x0;			y[11] = y0 - params->_RMIN;	t[11] = 0.0f;
	x[12] = x0;			y[12] = y0 + _LMH;				t[12] = -PIo2;
	x[13] = x0;			y[13] = y0 + _LMV;				t[13] = PIo2;

	//.. 2nd L of Halloween
	x0 = 4 * _S + 3 * _LW;
	y0 = _S;
	x[14] = x0;			y[14] = y0 - params->_RMIN;	t[14] = 0.0f;
	x[15] = x0;			y[15] = y0 + _LMH;				t[15] = -PIo2;
	x[16] = x0;			y[16] = y0 + _LMV;				t[16] = PIo2;

	//.. O of Halloween
	x0 = 5 * _S + 4 * _LW;
	y0 = _S;
	x[17] = x0;							y[17] = y0 - params->_RMIN;	t[17] = 0.0f;
	x[18] = x0;							y[18] = y0;						t[18] = PIo2;
	x[19] = x0;							y[19] = y0 + _LMV;				t[19] = PIo2;
	x[20] = x0 + params->_RMIN/ 2;	y[20] = y0 + _LT;				t[20] = 0.0f;
	x[21] = x0 + _LW;					y[21] = y0;						t[21] = PIo2;
	x[22] = x0 + _LW;					y[22] = y0 + _LMV;				t[22] = PIo2;

	//.. W of Halloween
	x0 = 6 * _S + 5 * _LW;
	y0 = _S;
	x[23] = x0;				y[23] = y0 - params->_RMIN;	t[23] = 0.0f;
	x[24] = x0 + _LW;		y[24] = y0 - params->_RMIN;	t[24] = 0.0f;
	x[25] = x0;				y[25] = y0;						t[25] = PIo2;
	x[26] = x0;				y[26] = y0 + _LMV;				t[26] = PIo2;
	x[27] = x0 + _LW;		y[27] = y0;						t[27] = PIo2;
	x[28] = x0 + 2 * _LW;	y[28] = y0;						t[28] = PIo2;
	x[29] = x0 + 2 * _LW;	y[29] = y0 + _LMV;				t[29] = PIo2;

	//.. 1st E of Halloween
	x0 = 7 * _S + 7 * _LW;
	y0 = _S;
	x[30] = x0;			y[30] = y0 - params->_RMIN;	t[30] = 0.0f;
	x[31] = x0;			y[31] = y0;						t[31] = PIo2;
	x[32] = x0;			y[32] = y0 + _LMH;				t[32] = 0.0f;
	x[33] = x0;			y[33] = y0 + _LMV;				t[33] = PIo2;
	x[34] = x0;			y[34] = y0 + _LT;				t[34] = 0.0f;

	//.. 2nd E of Halloween
	x0 = 8 * _S + 8 * _LW;
	y0 = _S;
	x[35] = x0;			y[35] = y0 - params->_RMIN;	t[35] = 0.0f;
	x[36] = x0;			y[36] = y0;						t[36] = PIo2;
	x[37] = x0;			y[37] = y0 + _LMH;				t[37] = 0.0f;
	x[38] = x0;			y[38] = y0 + _LMV;				t[38] = PIo2;
	x[39] = x0;			y[39] = y0 + _LT;				t[39] = 0.0f;

	//.. N of Halloween
	x0 = 9 * _S + 9 * _LW;
	y0 = _S;
	x[40] = x0;				y[40] = y0;				t[40] = PIo2;
	x[41] = x0;				y[41] = y0 + _LMH;		t[41] = PIo2;
	x[42] = x0;				y[42] = y0 + _LT;		t[42] = atan2f(-_LW, _LW / 2.0f);
	x[43] = x0 + _LW / 2.0;	y[43] = y0 + _LMH;		t[43] = atan2f(-_LW, _LW / 2.0f);
	x[44] = x0 + _LW;		y[44] = y0;				t[44] = PIo2;
	x[45] = x0 + _LW;		y[45] = y0 + _LMH;		t[45] = PIo2;

	//////////////////////// SPELLING HAPPY /////////////////////////////

	//.. H of Happy
	x0 = 3 * _S + 2 * _LW;
	y0 = _2ND;
	x[46] = x0;			y[46] = y0;			t[46] = PIo2;
	x[47] = x0;			y[47] = y0 + _LMH;	t[47] = 0.0f;
	x[48] = x0;			y[48] = y0 + _LMV;	t[48] = PIo2;
	x[49] = x0 + _LW;	y[49] = y0;			t[49] = PIo2;
	x[50] = x0 + _LW;	y[50] = y0 + _LMV;	t[50] = PIo2;

	//.. A fo Happy
	x0 = 4 * _S + 3 * _LW;
	y0 = _2ND;
	x[51] = x0;							y[51] = y0;			t[51] = PIo2;
	x[52] = x0;							y[52] = y0 + _LMH;	t[52] = 0.0f;
	x[53] = x0;							y[53] = y0 + _LMV;	t[53] = PIo2;
	x[54] = x0 + params->_RMIN / 2;	y[54] = y0 + _LT;	t[54] = 0.0f;
	x[55] = x0 + _LW;					y[55] = y0;			t[55] = PIo2;
	x[56] = x0 + _LW;					y[56] = y0 + _LMV;	t[56] = PIo2;

	//.. 1st P of Happy
	x0 = 5 * _S + 4 * _LW;
	y0 = _2ND;
	x[57] = x0;			y[57] = y0;			t[57] = PIo2;
	x[58] = x0;			y[58] = y0 + _LMH;	t[58] = 0.0f;
	x[59] = x0;			y[59] = y0 + _LT;	t[59] = 0.0f;
	x[60] = x0 + _LW;	y[60] = y0 + _LMV;	t[60] = PIo2;
	x[61] = x0;			y[61] = y0 + _LMV;	t[61] = PIo2;

	//.. 2nd P of Happy
	x0 = 6 * _S + 5 * _LW;
	y0 = _2ND;
	x[62] = x0;			y[62] = y0;			t[62] = PIo2;
	x[63] = x0;			y[63] = y0 + _LMH;	t[63] = 0.0f;
	x[64] = x0;			y[64] = y0 + _LT;	t[64] = 0.0f;
	x[65] = x0 + _LW;	y[65] = y0 + _LMV;	t[65] = PIo2;
	x[66] = x0;			y[66] = y0 + _LMV;	t[66] = PIo2;

	//.. Y of Happy
	x0 = 7 * _S + 6 * _LW;
	y0 = _2ND;
	x[67] = x0 + _LW / 2;						y[67] = y0;			t[67] = PIo2;
	x[68] = x0 + _LW / 2 - params->_RMIN / 2;	y[68] = y0 + _LMH;	t[68] = 2.35619f;
	x[69] = x0 + _LW / 2 + params->_RMIN / 2;	y[69] = y0 + _LMH;	t[69] = 0.785398f;

	worms->CustomInit(x, y, t);
	delete[] x, y, t;
}

void HalloweenSimulation::XYZPrint(int itime){
	//.. only print when it itime == 0
	if (itime % simparams->_FRAMERATE != 0) return;

	//.. count and grab errors
	static int frame = 1;
	printf("\nPrinting frame %i", frame++);
	this->worms->DislayTheta();
	this->DisplayErrors();
	this->worms->DataDeviceToHost();
	this->DisplayErrors();

	//.. print to ntypes
	static const int ntypes = 4;
	static const char ptypes[ntypes] = { 'A', 'B', 'C', 'D' };
	this->fxyz << params->_NPARTICLES  + 4 << std::endl;
	this->fxyz << "Happy Halloween!!" << std::endl;
	for (int i = 0; i < params->_NPARTICLES; i++){
		int t = (i / params->_NP) % ntypes;
		this->fxyz << ptypes[t] << " " << worms->X[i] << " " << worms->Y[i] << " 0 " << std::endl;
	}
	this->fxyz << "E " << 0 << " " << 0 << " 0 " << std::endl;
	this->fxyz << "E " << simparams->_XBOX << " " << 0 << " 0 " << std::endl;
	this->fxyz << "E " << 0 << " " << simparams->_YBOX << " 0 " << std::endl;
	this->fxyz << "E " << simparams->_XBOX << " " << simparams->_YBOX << " 0 " << std::endl;
}

void HalloweenSimulation::DisplayErrors(){
	this->rng->DisplayErrors();
	this->worms->DisplayErrors();
}

#endif