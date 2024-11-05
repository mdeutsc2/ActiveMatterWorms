// Functions only for execution on GPU device
// 5.12.15
// Mike Varga
// 2D
#ifndef _DEVICE_FUNCTIONS_H
#define _DEVICE_FUNCTIONS_H
#include "NWmain.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
//-----------------------------------------------------------------------------------
__device__ void BC_dr(float &dR, float L, const int dim)
{
	if (dev_simParams._PBC[dim]){
		if (dR > L / 2.0f) dR -= L;
		else if (dR < -L / 2.0f) dR += L;
	}
}
//-----------------------------------------------------------------------------------
__device__ void BC_dr(float _dR[_D_], float _box[_D_]){
	for_D_ BC_dr(_dR[d], _box[d], d);
}
//-----------------------------------------------------------------------------------
__device__ void BC_r(float &F, float &R, float L, const int dim) // actually implements the BCs
{
	if (dev_simParams._PBC[dim]){ // move to otherside of box
		if (R > L) R -= L;
		else if (R < 0) R += L;
	}
	else if (dev_simParams._SBC[dim]){ // apply force from soft wall
		if (R > L || R < 0.0f) F -= dev_simParams._Kw * (R - L);
	}
}
//-----------------------------------------------------------------------------------
__device__ void BC_r(float _F[_D_], float _R[_D_], float _box[_D_]){
	for_D_ BC_r(_F[d], _R[d], _box[d], d);
}
//-----------------------------------------------------------------------------------
__device__ bool InList(int arg, int* list, int listSize)
{
	for (int i = 0; i < listSize; i++)
	{
		if (arg == list[i]) return true;
	}
	return false;
}
// ----------------------------------------------------------------------------------
__device__ float dot(const float v1[_D_], const float v2[_D_]){
	float result = 0.0f;
	for_D_ result += (v1[d] * v2[d]);
	return result;
}
// --------------------------------------------------------------------------------
__device__ float mag(const float v[_D_]){
	return sqrt(dot(v,v));
}
//-----------------------------------------------------------------------------------
__device__ void diff_dir(const float this_v[_D_], const float minus_v[_D_], float _unit_dir[_D_]){
	float _dif[_D_];
	for_D_ _dif[d] = this_v[d] - minus_v[d];
	BC_dr(_dif, dev_simParams._BOX);
	float _mag = mag(_dif);
	for_D_ _unit_dir[d] = _dif[d] / _mag;
}
//-----------------------------------------------------------------------------------
__device__ float CalculateRR(const float _rid[_D_], const float _rnab[_D_], float _dr[_D_]){
	float _r_[_D_];
	for_D_ _r_[d] = _rnab[d] - _rid[d];
	BC_dr(_r_, dev_simParams._BOX);
	for_D_ _dr[d] = _r_[d];
	return dot(_r_, _r_);
}
//-----------------------------------------------------------------------------------
__device__ float CalculateLJ(const float _rr){
	return dev_Params._LJ_AMP*((dev_Params._2SIGMA6 / (_rr*_rr*_rr*_rr*_rr*_rr*_rr)) - (1.00000f / (_rr*_rr*_rr*_rr)));
}
//-----------------------------------------------------------------------------------
__device__ float CalculateLJ(const float _rr, const float sig, const float eps){
	const float _2sig6 = pow(sig, 6.0f);
	const float _lj_amp = 12.0f * eps * _2sig6;
	return _lj_amp*((_2sig6 / (_rr*_rr*_rr*_rr*_rr*_rr*_rr)) - (1.00000f / (_rr*_rr*_rr*_rr)));
}
//---------------------------------------------------------------------------------
__device__ void Rotate2D(float _v[2], const float theta){
	const float costh = cosf(theta), sinth = sinf(theta);
	const float R[2][2] = {
		{ costh, -sinth },
		{ sinth, costh }
	};
	float _nv[2] = { 0, 0 }; 
	for (int i = 0; i < 2; i++) // for x and y
		for (int j = 0; j < 2; j++) // sum over j
		_nv[i] += R[i][j] * _v[j];
	
	_v[0] = _nv[0]; _v[1] = _nv[1]; // assign
}
#endif