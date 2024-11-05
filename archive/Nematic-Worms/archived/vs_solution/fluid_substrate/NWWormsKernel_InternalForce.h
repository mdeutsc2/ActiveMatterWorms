#ifndef __WORMS_KERNEL__INTERNAL_FORCE_H__
#define __WORMS_KERNEL__INTERNAL_FORCE_H__
// 2D
// ------------------------------------------------------------------------------------------
#include "NWDeviceFunctions.h"
#include "NWParams.h"
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"

// ------------------------------------------------------------------------------------------
__global__ void InterForceKernel(float *f,
								 int fshift,
								 float *v,
								 int vshift,
								 float *r,
								 int rshift,
								 float *randNum,
								 float noiseScaler)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id < dev_Params._NPARTICLES){

		//.. particel number in worm
		int p = id % dev_Params._NP;
		int w = id / dev_Params._NP;

		//.. local memory
		float rid[_D_], fid[_D_], vid[_D_];

		//.. init local as needed
		for_D_{
			rid[d] = r[id + d*rshift];
			fid[d] = 0.0f;
			vid[d] = v[id + d*vshift];
		}
			
		//.. 1st neighbor spring forces ahead
		if (p < (dev_Params._NP - 1))
		{
			int pp1 = id + 1;
			float rnab[_D_], dr[_D_];
			float _r, _f;
			for_D_ rnab[d] = r[pp1 + d*rshift];
			_r = sqrt(CalculateRR(rid, rnab, dr));
			_f = -(dev_Params._K1 * (_r - dev_Params._L1)) / _r;
			for_D_ fid[d] -= _f * dr[d];
		}

		//.. 1st neighbor spring forces behind
		if (p > 0)
		{
			int pm1 = id - 1;
			float rnab[_D_], dr[_D_];
			float _r, _f;
			for_D_ rnab[d] = r[pm1 + d*rshift];
			_r = sqrt(CalculateRR(rid, rnab, dr));
			_f = -(dev_Params._K1 * (_r - dev_Params._L1)) / _r;
			for_D_ fid[d] -= _f * dr[d];
		}

		//.. 2nd neighbor spring forces ahead
//		if (p < (dev_Params._NP - 2))
//		{
//			int pp2 = id + 2;
//			float rnab[3];
//			float dr[3];
//			float _r, _f;
//			rnab[0] = r[pp2 + 0 * rshift];
//			rnab[1] = r[pp2 + 1 * rshift];
//			rnab[2] = r[pp2 + 2 * rshift];
//			_r = sqrtf(CalculateRR_3d(rid, rnab, dr));
//			_f = -(dev_Params._K2 * (_r - dev_Params._L2)) / _r;
//			for (int d = 0; d < 3; d++)
//				fid[d] -= _f * dr[d];
//			//float dx = x[pp2] - xid;
//			//float dy = y[pp2] - yid;
//			//DevicePBC(dx, dev_simParams._XBOX);
//			//DevicePBC(dy, dev_simParams._YBOX);
//			//float rr = dx*dx + dy*dy;
//			//float _r = sqrtf(rr);
//			//float ff = -dev_Params._K2 * (_r - dev_Params._L2) / _r;
//			//float ffx = ff * dx;
//			//float ffy = ff * dy;
//			//fxid -= ffx;
//			//fyid -= ffy;
////#ifdef _DAMPING1
////			float dvx = vx[pp2] - vxid;
////			float dvy = vy[pp2] - vyid;
////			fxid += dev_Params._DAMP * dvx;
////			fyid += dev_Params._DAMP * dvy;
////#endif
//		}
//
//		//.. 2nd neighbor spring forces behind
//		if (p > 1)
//		{
//			int pm2 = id - 2;
//			float rnab[3];
//			float dr[3];
//			float _r, _f;
//			rnab[0] = r[pm2 + 0 * rshift];
//			rnab[1] = r[pm2 + 1 * rshift];
//			rnab[2] = r[pm2 + 2 * rshift];
//			_r = sqrtf(CalculateRR_3d(rid, rnab, dr));
//			_f = -(dev_Params._K2 * (_r - dev_Params._L2)) / _r;
//			for (int d = 0; d < 3; d++)
//				fid[d] -= _f * dr[d];
//			//float dx = x[pm2] - xid;
//			//float dy = y[pm2] - yid;
//			//DevicePBC(dx, dev_simParams._XBOX);
//			//DevicePBC(dy, dev_simParams._YBOX);
//			//float rr = dx*dx + dy*dy;
//			//float _r = sqrtf(rr);
//			//float ff = -dev_Params._K2 * (_r - dev_Params._L2) / _r;
//			//float ffx = ff * dx;
//			//float ffy = ff * dy;
//			//fxid -= ffx;
//			//fyid -= ffy;
////#ifdef _DAMPING1
//			//float dvx = vx[pm2] - vxid;
//			//float dvy = vy[pm2] - vyid;
//			//fxid += dev_Params._DAMP * dvx;
//			//fyid += dev_Params._DAMP * dvy;
////#endif
//		}

		//.. 3nd neighbor spring forces ahead
//		if (p < (dev_Params._NP - 3))
//		{
//			int pp3 = id + 3;
//			float rnab[3];
//			float dr[3];
//			float _r, _f;
//			rnab[0] = r[pp3 + 0 * rshift];
//			rnab[1] = r[pp3 + 1 * rshift];
//			rnab[2] = r[pp3 + 2 * rshift];
//			_r = sqrtf(CalculateRR_3d(rid, rnab, dr));
//			_f = -(dev_Params._K3 * (_r - dev_Params._L3)) / _r;
//			for (int d = 0; d < 3; d++)
//				fid[d] -= _f * dr[d];
//			//float dx = x[pp3] - xid;
//			//float dy = y[pp3] - yid;
//			//DevicePBC(dx, dev_simParams._XBOX);
//			//DevicePBC(dy, dev_simParams._YBOX);
//			//float rr = dx*dx + dy*dy;
//			//float _r = sqrtf(rr);
//			//float ff = -dev_Params._K3 * (_r - dev_Params._L3) / _r;
//			//float ffx = ff * dx;
//			//float ffy = ff * dy;
//			//fxid -= ffx;
//			//fyid -= ffy;
//
////#ifdef _DAMPING2
////			float dvx = vx[pp3] - vxid;
////			float dvy = vy[pp3] - vyid;
////			fxid += dev_Params._DAMP * dvx;
////			fyid += dev_Params._DAMP * dvy;
////#endif
//		}
//
//		//.. 3nd neighbor spring forces behind
//		if (p > 2)
//		{
//			int pm3 = id - 3;
//			float rnab[3];
//			float dr[3];
//			float _r, _f;
//			rnab[0] = r[pm3 + 0 * rshift];
//			rnab[1] = r[pm3 + 1 * rshift];
//			rnab[2] = r[pm3 + 2 * rshift];
//			_r = sqrtf(CalculateRR_3d(rid, rnab, dr));
//			_f = -(dev_Params._K3 * (_r - dev_Params._L3)) / _r;
//			for (int d = 0; d < 3; d++)
//				fid[d] -= _f * dr[d];
//			//float dx = x[pm3] - xid;
//			//float dy = y[pm3] - yid;
//			//DevicePBC(dx, dev_simParams._XBOX);
//			//DevicePBC(dy, dev_simParams._YBOX);
//			//float rr = dx*dx + dy*dy;
//			//float _r = sqrtf(rr);
//			//float ff = -dev_Params._K3 * (_r - dev_Params._L3) / _r;
//			//float ffx = ff * dx;
//			//float ffy = ff * dy;
//			//fxid -= ffx;
//			//fyid -= ffy;
//
////#ifdef _DAMPING2
////			float dvx = vx[pm3] - vxid;
////			float dvy = vy[pm3] - vyid;
////			fxid += dev_Params._DAMP * dvx;
////			fyid += dev_Params._DAMP * dvy;
////#endif
//		}

		//.. LJ between inter-worm particles
		for (int id2 = w*dev_Params._NP; id2 < (w + 1)*dev_Params._NP; id2++)
		{
			int p2 = id2 % dev_Params._NP;
			int sep = abs(p2 - p);
			if (sep <= 2) continue; //.. ignore 1st and 2nd neighbors
			float rnab[3], dr[3];
			float rr, _f;
			for_D_ rnab[d] = r[id2 + d*rshift];
			rr = CalculateRR(rid, rnab, dr);
			if (rr > dev_Params._R2MIN) continue; //.. repulsive only
			_f = CalculateLJ(rr);
			for_D_ fid[d] -= _f * dr[d];
		}

		//.. viscous drag
		for_D_ fid[d] -= dev_Params._GAMMA * vid[d];

		//.. thermal fluctuations
		for_D_ fid[d] += noiseScaler * randNum[id + d*dev_Params._NPARTICLES];

		//.. assign temp fxid and fyid to memory
		for_D_ f[id + d*fshift] += fid[d];
	}
}

#endif