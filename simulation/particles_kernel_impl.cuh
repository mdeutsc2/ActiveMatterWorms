/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/*
 * CUDA particle system kernel code.
 */

#ifndef _PARTICLES_KERNEL_H_
#define _PARTICLES_KERNEL_H_

#include <stdio.h>
#include <math.h>
#include <cooperative_groups.h>

namespace cg = cooperative_groups;
#include "helper_math.h"
#include "math_constants.h"
#include "particles_kernel.cuh"

// simulation parameters in constant memory
__constant__ SimParams params;
__constant__ SolventParams srd;



__device__ void periodicBoundary(float &r, const float &L0, const float &L)
{
    if (r < L0) r += L;
    else if (r > L + L0) r -= L;
}

__device__ void dampedWallBoundary(float &r, float &v, const float &coef, const float &L0, const float &L, const float &R)
{
    if (r < L0 + R)
    {
        r = L0 + R;
        v *= coef;
    }
    else if (r > (L + L0) - R)
    {
        r = (L + L0) - R;
        v *= coef;
    }
}

__device__ void noslipWallBoundary(float3 &pos, float &r, float3 &vel, const float3 &dr, const float3 N, const float &L0, const float &L, const float &R)
{
    if (r < L0 + R)
    {
        pos -= dr;
        r = L0 + R;
        vel = make_float3(-vel.x * N.x, -vel.y * N.y, -vel.z * N.z);
    }
    else if (r > (L + L0) - R)
    {
        pos -= dr;
        r = (L + L0) - R;
        vel = make_float3(-vel.x * N.x, -vel.y * N.y, -vel.z * N.z);
    }
}


struct filament_integrator
{
    float dt;

    __host__ __device__
    filament_integrator(float delta_time) : dt(delta_time) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        volatile float4 posData = thrust::get<0>(t);
        volatile float4 velData = thrust::get<1>(t);
        volatile float4 forceData = thrust::get<2>(t);
        volatile float4 forceOldData = thrust::get<3>(t);
        float3 pos   = make_float3(posData);
        float3 vel   = make_float3(velData);
        float3 f     = make_float3(forceData);
        float3 f_old = make_float3(forceOldData);

        // f += params.gravity;

        // Velocity Verlet update
        float3 dv = 0.5f * (f + f_old) * dt;
        vel += dv;
        float3 dr = vel * dt + 0.5f * f * dt * dt;
        pos += dr;

        if (params.boundaryX == BoundaryType::PERIODIC)
        {
            periodicBoundary(pos.x, params.origin.x, params.boxSize.x);
        }
        else if(params.boundaryX == BoundaryType::WALL)
        {
            dampedWallBoundary(pos.x, vel.x, params.boundaryDamping, params.origin.x, params.boxSize.x, params.particleRadius);
        }
        else if(params.boundaryX == BoundaryType::WALL_NO_SLIP)
        {
            const float3 N = make_float3(1, 0, 0);
            noslipWallBoundary(pos, pos.x, vel, dr, N, params.origin.x, params.boxSize.x, params.particleRadius);
        }

        if (params.boundaryY == BoundaryType::PERIODIC)
        {
            periodicBoundary(pos.y, params.origin.y, params.boxSize.y);
        }
        else if(params.boundaryY == BoundaryType::WALL)
        {
            dampedWallBoundary(pos.y, vel.y, params.boundaryDamping, params.origin.y, params.boxSize.y, params.particleRadius);
        }
        else if(params.boundaryY == BoundaryType::WALL_NO_SLIP)
        {
            const float3 N = make_float3(0, 1, 0);
            noslipWallBoundary(pos, pos.y, vel, dr, N, params.origin.y, params.boxSize.y, params.particleRadius);
        }

        if (params.boundaryZ == BoundaryType::PERIODIC)
        {
            periodicBoundary(pos.z, params.origin.z, params.boxSize.z);
        }
        else if(params.boundaryZ == BoundaryType::WALL)
        {
            dampedWallBoundary(pos.z, vel.z, params.boundaryDamping, params.origin.z, params.boxSize.z, params.particleRadius);
        }

        // store new position, velocity, and forces
        thrust::get<0>(t) = make_float4(pos, posData.w);
        thrust::get<1>(t) = make_float4(vel, velData.w);
        thrust::get<2>(t) = make_float4(0.f);
        thrust::get<3>(t) = make_float4(f, forceOldData.w);
    }
};


struct solvent_integrator
{
    float dt;

    __host__ __device__
    solvent_integrator(float delta_time) : dt(delta_time) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        volatile float4 posData = thrust::get<0>(t);
        volatile float4 velData = thrust::get<1>(t);
        float3 pos = make_float3(posData);
        float3 vel = make_float3(velData);

        float3 dr = dt * vel;
        pos += dr;

        if (params.boundaryX == BoundaryType::PERIODIC)
        {
            periodicBoundary(pos.x, params.origin.x, params.boxSize.x);
        }
        else if(params.boundaryX == BoundaryType::WALL)
        {
            dampedWallBoundary(pos.x, vel.x, params.boundaryDamping, params.origin.x, params.boxSize.x, 0.0f);
        }
        else if(params.boundaryX == BoundaryType::WALL_NO_SLIP)
        {
            const float3 N = make_float3(1, 0, 0);
            noslipWallBoundary(pos, pos.x, vel, dr, N, params.origin.x, params.boxSize.x, 0.0f);
        }

        if (params.boundaryY == BoundaryType::PERIODIC)
        {
            periodicBoundary(pos.y, params.origin.y, params.boxSize.y);
        }
        else if(params.boundaryY == BoundaryType::WALL)
        {
            dampedWallBoundary(pos.y, vel.y, params.boundaryDamping, params.origin.y, params.boxSize.y, 0.0f);
        }
        else if(params.boundaryY == BoundaryType::WALL_NO_SLIP)
        {
            const float3 N = make_float3(0, 1, 0);
            noslipWallBoundary(pos, pos.y, vel, dr, N, params.origin.y, params.boxSize.y, 0.0f);
        }

        if (params.boundaryZ == BoundaryType::PERIODIC)
        {
            periodicBoundary(pos.z, params.origin.z, params.boxSize.z);
        }
        else if(params.boundaryZ == BoundaryType::WALL)
        {
            dampedWallBoundary(pos.z, vel.z, params.boundaryDamping, params.origin.z, params.boxSize.z, 0.0f);
        }

        // store new position, velocity, and forces
        thrust::get<0>(t) = make_float4(pos, posData.w);
        thrust::get<1>(t) = make_float4(vel, velData.w);
    }
};


__device__ float3 positionRelativeToBox(float3 p)
{
    float3 pos = make_float3(p.x - params.origin.x, p.y - params.origin.y, p.z - params.origin.z);

    if (params.boundaryX == BoundaryType::PERIODIC)
    {
        if (pos.x >= params.boxSize.x) pos.x -= params.boxSize.x;
        else if (pos.x < 0.0f) pos.x += params.boxSize.x;
    }
    if (params.boundaryY == BoundaryType::PERIODIC)
    {
        if (pos.y >= params.boxSize.y) pos.y -= params.boxSize.y;
        else if (pos.y < 0.0f) pos.y += params.boxSize.y;
    }
    if (params.boundaryZ == BoundaryType::PERIODIC)
    {
        if (pos.z >= params.boxSize.z) pos.z -= params.boxSize.z;
        else if (pos.z < 0.0f) pos.z += params.boxSize.z;
    }
    return pos;
}

// calculate position in uniform grid
__device__ int3 calcGridPos(float3 p)
{
    int3 gridPos;
    gridPos.x = floor((p.x - params.origin.x) / params.cellSize.x);
    gridPos.y = floor((p.y - params.origin.y) / params.cellSize.y);
    gridPos.z = floor((p.z - params.origin.z) / params.cellSize.z);
    return gridPos;
}

// calculate position in uniform grid (SRD)
__device__ int3 calcGridPosSRD(float3 p)
{
    int3 gridPos;
    float3 pos = positionRelativeToBox(p + srd.randOffset);
    gridPos.x = floor(pos.x / srd.cellSize.x);
    gridPos.y = floor(pos.y / srd.cellSize.y);
    gridPos.z = floor(pos.z / srd.cellSize.z);
    return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__ uint calcGridHash(int3 gridPos)
{
    gridPos.x = gridPos.x & (params.gridSize.x-1);  // wrap grid, assumes size is power of 2
    gridPos.y = gridPos.y & (params.gridSize.y-1);
    gridPos.z = gridPos.z & (params.gridSize.z-1);
    return __umul24(__umul24(gridPos.z, params.gridSize.y), params.gridSize.x) + __umul24(gridPos.y, params.gridSize.x) + gridPos.x;
}

// calculate address in grid from position (clamping to edges)
__device__ uint calcGridHashSRD(int3 gridPos)
{
    if (params.boundaryX == BoundaryType::PERIODIC)
    {
        gridPos.x = gridPos.x & (srd.gridSize.x-1);  // wrap grid, assumes size is power of 2
    }
    if (params.boundaryY == BoundaryType::PERIODIC)
    {
        gridPos.y = gridPos.y & (srd.gridSize.y-1);
    }
    if (params.boundaryZ == BoundaryType::PERIODIC)
    {
        gridPos.z = gridPos.z & (srd.gridSize.z-1);
    }
    // TODO: change __umul24->__umul32 for CC>2.0
    return __umul24(__umul24(gridPos.z, srd.gridSize.y), srd.gridSize.x) + __umul24(gridPos.y, srd.gridSize.x) + gridPos.x;
}

// calculate address in grid from position (clamping to edges)
__device__ float lengthPeriodic(float3& dr)
{
    float3 L = params.boxSize;
    float3 Lo2 = L / 2.0f;

    if (params.boundaryX == BoundaryType::PERIODIC)
    {
        if (dr.x > Lo2.x) dr.x -= L.x;
        else if (dr.x < -Lo2.x) dr.x += L.x;
    }
    if (params.boundaryY == BoundaryType::PERIODIC)
    {
        if (dr.y > Lo2.y) dr.y -= L.y;
        else if (dr.y < -Lo2.y) dr.y += L.y;
    }
    if (params.boundaryZ == BoundaryType::PERIODIC)
    {
        if (dr.z > Lo2.z) dr.z -= L.z;
        else if (dr.z < -Lo2.z) dr.z += L.z;
    }

    return length(dr);
}

__device__
float3 getTangent(float3 posA, float3 posB)
{
    float3 relPos = posB - posA;
    float dist = lengthPeriodic(relPos);
    return relPos / dist;
}

// calculate grid hash value for each particle
__global__
void calcHashD(uint   *gridParticleHash,  // output
               uint   *gridParticleIndex, // output
               float4 *pos,               // input: positions
               uint    numParticles,
               bool    srd = false)
{
    uint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    volatile float4 p = pos[index];
    int3 gridPos;
    uint hash;

    // get address in grid
    if (srd)
    {
        gridPos = calcGridPosSRD(make_float3(p.x, p.y, p.z));
        hash = calcGridHashSRD(gridPos);
    }
    else
    {
        gridPos = calcGridPos(make_float3(p.x, p.y, p.z));
        hash = calcGridHash(gridPos);
    }

    // store grid hash and particle index
    gridParticleHash[index] = hash;
    gridParticleIndex[index] = index;
}

// rearrange particle data into sorted order, and find the start of each cell
// in the sorted hash array
__global__
void reorderDataAndFindCellStartD(uint   *cellStart,        // output: cell start index
                                  uint   *cellEnd,          // output: cell end index
                                  float4 *sortedPos,        // output: sorted positions
                                  float4 *sortedVel,        // output: sorted velocities
                                  float4 *sortedTangent,    // output: sorted tangents
                                  uint   *gridParticleHash, // input: sorted grid hashes
                                  uint   *gridParticleIndex,// input: sorted particle indices
                                  float4 *pos,              // input: position array
                                  float4 *vel,              // input: velocity array
                                  float4 *tangent,          // input: tangents array
                                  uint    numParticles)
{
    // Handle to thread block group
    cg::thread_block cta = cg::this_thread_block();
    extern __shared__ uint sharedHash[];    // blockSize + 1 elements
    uint index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

    uint hash;

    // handle case when no. of particles not multiple of block size
    if (index < numParticles)
    {
        hash = gridParticleHash[index];

        // Load hash data into shared memory so that we can look
        // at neighboring particle's hash value without loading
        // two hash values per thread
        sharedHash[threadIdx.x+1] = hash;

        if (index > 0 && threadIdx.x == 0)
        {
            // first thread in block must load neighbor particle hash
            sharedHash[0] = gridParticleHash[index-1];
        }
    }

    cg::sync(cta);

    if (index < numParticles)
    {
        // If this particle has a different cell index to the previous
        // particle then it must be the first particle in the cell,
        // so store the index of this particle in the cell.
        // As it isn't the first particle, it must also be the cell end of
        // the previous particle's cell

        if (index == 0 || hash != sharedHash[threadIdx.x])
        {
            cellStart[hash] = index;

            if (index > 0)
                cellEnd[sharedHash[threadIdx.x]] = index;
        }

        if (index == numParticles - 1)
        {
            cellEnd[hash] = index + 1;
        }

        // Now use the sorted index to reorder data array
        uint sortedIndex = gridParticleIndex[index];
        sortedPos[index] = pos[sortedIndex];

        // optionally sort velocities
        if (vel != NULL && sortedVel != NULL)
        {
            sortedVel[index] = vel[sortedIndex];
        }

        // optionally sort filament tangents per particle
        if (tangent != NULL && sortedTangent != NULL)
        {
            sortedTangent[index] = tangent[sortedIndex];
        }
    }


}

// collide two spheres using DEM method
__device__
float3 collideParticles(float3 posA, float3 posB,
                      float3 velA, float3 velB,
                      float radiusA, float radiusB,
                      float attraction)
{
    // calculate relative position
    float3 relPos = posB - posA;

    float dist = lengthPeriodic(relPos);
    float collideDist = radiusA + radiusB;

    float3 force = make_float3(0.0f);

    if (dist < collideDist)
    {
        float3 norm = relPos / dist;

        // relative velocity
        float3 relVel = velB - velA;

        // relative tangential velocity
        float3 tanVel = relVel - (dot(relVel, norm) * norm);

        // spring force
        force = -params.spring*(collideDist - dist) * norm;
        // dashpot (damping) force
        force += params.damping*relVel;
        // tangential shear force
        force += params.shear*tanVel;
        // attraction
        // force += attraction*relPos;
    }

    return force;
}


// bond to particles with Hookean spring
__device__
float3 bondHookean(float3 posA, float3 posB, float k, float L)
{
    float3 relPos = posB - posA;
    float dist = lengthPeriodic(relPos);
    float3 norm = relPos / dist;
    return -k * (L - dist) * norm;
}


// compute extensile drive forces
__device__
float3 extensileForce(float3 posA, float3 posB, float3 tangA, float3 tangB, float A)
{
    float3 relPos = posB - posA;
    float dist = lengthPeriodic(relPos);
    if (dot(tangA, tangB) < -0.5f)
    {
        return -A * (tangA - tangB) / (2.0f * dist);
    }
    return make_float3(0.0f);
}



// collide a particle against all other particles in a given cell
__device__
float3 collideCell(int3    gridPos,
                   uint    index,
                   float3  r,
                   float3  v,
                   float3  t,
                   float4 *pos,
                   float4 *vel,
                   float4 *tangent,
                   uint   *cellStart,
                   uint   *cellEnd,
                   uint   *gridParticleIndex)
{
    uint gridHash = calcGridHash(gridPos);

    // get start of bucket for this cell
    uint  startIndex     = cellStart[gridHash];
    uint  filamentSize   = params.filamentSize;
    uint  filamentIndex  = gridParticleIndex[index] / filamentSize;
    uint  chainIndex     = gridParticleIndex[index] % filamentSize;
    uint  filamentIndex2, chainIndex2;
    bool  sameFilament;

    float3 force = make_float3(0.0f);

    if (startIndex != 0xffffffff)          // cell is not empty
    {
        // iterate over particles in this cell
        uint endIndex = cellEnd[gridHash];

        for (uint j=startIndex; j<endIndex; j++)
        {
            if (j != index)                // check not colliding with self
            {
                float3 r2      = make_float3(pos[j]);
                float3 v2      = make_float3(vel[j]);
                float3 t2      = make_float3(tangent[j]);
                filamentIndex2 = gridParticleIndex[j] / filamentSize;
                chainIndex2    = gridParticleIndex[j] % filamentSize;
                sameFilament   = filamentIndex == filamentIndex2;

                // filament bonding
                if (filamentIndex == filamentIndex2 && ((chainIndex + 1 == chainIndex2) || (chainIndex - 1 == chainIndex2)))
                {
                    force += bondHookean(r, r2, params.bondSpringK, params.bondSpringL);
                }
                else // collide particles
                {
                    force += collideParticles(r, r2, v, v2, params.particleRadius, params.particleRadius, params.attraction);
                }
                // if (!sameFilament || (chainIndex > chainIndex2 + 1) || (chainIndex < chainIndex2 - 1))
                // {
                //     force += collideParticles(r, r2, v, v2, particleRadius, particleRadius, 0.0f);
                // }

                if (!sameFilament)
                {
                    force += extensileForce(r, r2, t, t2, params.activity);
                }
            }
        }
    }

    return force;
}


// Launch one per solvent cell
__global__
void cellCenterOfMomentumSRD(float4 *cellCOM,
                             float4 *sortedVel,
                             uint   *cellStart,
                             uint   *cellEnd,
                             uint   *gridParticleIndex,
                             uint    numCells)
{
    uint cellHash = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (cellHash>= numCells) return;

    uint startIndex = cellStart[cellHash]; // index in sorted arrays where cell data starts
    float3 com      = make_float3(0.0f);

    if (startIndex != 0xffffffff)
    {
        uint endIndex = cellEnd[cellHash]; // index in sorted arrays where cell data end
        float M = 0;
        for (uint i = startIndex; i < endIndex; i++) // Iterate from start to end
        {
            float3 v = make_float3(sortedVel[i]);
            float m  = sortedVel[i].w; // mass is stored in w-component of float4
            com += m * v; // aggregate momentum
            M += m;
        }

        if (M > 0) // take mean
        {
            com *= (1.0f / M);
        }
    }

    cellCOM[cellHash] = make_float4(com, 0.0f);
}


__device__
float kineticEnergy(float m, float3 v)
{
    return 0.5f * m * dot(v, v);
}


// Launch one per solvent cell
__global__
void isokineticThermostatSRD(float4 *vel,
                             float4 *sortedVel,
                             uint   *cellStart,
                             uint   *cellEnd,
                             uint   *gridParticleIndex,
                             uint    numCells)
{
    uint cellHash = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (cellHash>= numCells) return;

    uint startIndex = cellStart[cellHash]; // index in sorted arrays where cell data starts
    if (startIndex != 0xffffffff)
    {
        uint endIndex = cellEnd[cellHash]; // index in sorted arrays where cell data end
        float E_kinetic = 0;
        for (uint i = startIndex; i < endIndex; i++)
        {
            E_kinetic += kineticEnergy(sortedVel[i].w, make_float3(sortedVel[i]));
        }

        const float lambda = sqrtf(srd.kbT / E_kinetic);

        uint index;
        for (uint i = startIndex; i < endIndex; i++)
        {
            index = gridParticleIndex[i];
            vel[index] = make_float4(lambda * make_float3(vel[index]), sortedVel[i].w);
            sortedVel[i] = vel[index];
        }
    }
}


__device__ float3 rotate2D(float3 v, const float alpha)
{
    const float sin_alpha = sin(alpha), cos_alpha = cos(alpha);
    return make_float3(
        v.x * cos_alpha - v.y * sin_alpha,
        v.x * sin_alpha + v.y * cos_alpha,
        v.z
    );
}


// Launch one for ALL particles (filaments + solvent)
// 2D-XY only!
__global__
void collideSolventKernel(float4 *pos,float4 *vel, float4 *vforces, uint *cellStart, uint *cellEnd, uint numParticles) {
    //TODO: add cellStart, cellEnd, gridParticleIndex
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;
    if (index < params.numParticles) return; // only update solvent particles

    // read particle data from sorted arrays
    float3 r1 = make_float3(pos[index]);
    float3 v = make_float3(vel[index]);
    int3 gridPos = calcGridPos(r1);
    float3 force = make_float3(0.0f);

    //pre-calculating some coefficients
    float eps = 1.3806e-16;
    float sigma = 3.4;
    float constA = (48*eps)/(sigma*sigma);

    for (int z=-1;z<=1;z++) {
        for (int y=-1;y<=1;y++) {
            for (int x=-1;x<=1;x++) {
                int3 neighborPos = gridPos + make_int3(x,y,z);
                // handling collisions within cell and neighboring cells
                uint gridHash = calcGridHash(neighborPos);
                // get start of bucked for this cell
                uint startIndex = cellStart[gridHash];
                if (startIndex != 0xffffffff){
                    // iterate over particles in this cell
                    uint endIndex = cellEnd[gridHash];
                    for (uint j = startIndex; j<=endIndex; j++){
                        if (j!=index) {
                            // not colliding with self
                            float3 r2 = make_float3(pos[j]);
                            float3 v2 = make_float3(vel[j]);
                            // colliding particles
                            float3 rel_pos = r2 - r1;
                            float dist = lengthPeriodic(rel_pos);
                            float collide_dist = 2*params.particleRadius;
                            if (dist < collide_dist) {
                                float sorij = sigma/dist;
                                force += constA*dist*(pow(sorij,14.0) - 0.5*pow(sorij,8.0));
                            }


                        }
                    }
                }
                force += 0.0f; // collideSolvent
            }
        }
    }
    // write velocity back
    
}
__global__
void collideSolventKernel_old(float4 *pos,
                          float4 *vel,
                          float4 *cellCOM,
                          float  *uniform,
                          uint    numParticles)
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    // only update solvent particles
    if (index < params.numParticles) return;

    float3 r = make_float3(pos[index]);
    float3 v = make_float3(vel[index]);

    int3 gridPos = calcGridPosSRD(r);
    uint cellIndex = calcGridHashSRD(gridPos);
    float3 CoM = make_float3(cellCOM[cellIndex]);

    // The SRD collision equation.
    int sign = 2 * int(uniform[cellIndex] > 0.5f) - 1;
    v = CoM + rotate2D(v - CoM, sign * srd.alpha);

    // if (index == 0)
    // {
    //     int sign = 2 * int(uniform[cellIndex] > 0.5f) - 1;
    //     printf("[%i] in cell %i\n", index, cellIndex);
    //     printf("  v_before = {%f %f}\n", v.x, v.y);
    //     printf("  CoM      = {%f %f}\n", CoM.x, CoM.y);
    //     printf("  sign     = %i\n", sign);

    //     v = CoM + rotate2D(v - CoM, sign * srd.alpha);
    //     printf("  v_after  = {%f %f}\n", v.x, v.y);
    // }

    vel[index] = make_float4(v, vel[index].w);
}


__global__
void collideFilamentsKernel(float4 *forces,          // update: unsorted forces array
                            float4 *sortedPos,            // input: sorted positions
                            float4 *sortedVel,            // input: sorted velocities
                            float4 *sortedTangent,        // input: sorted tangents
                            uint   *gridParticleIndex,    // input: sorted particle indices
                            uint   *cellStart,
                            uint   *cellEnd,
                            uint    numParticles)
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    // read particle data from sorted arrays
    float3 pos     = make_float3(sortedPos[index]);
    float3 vel     = make_float3(sortedVel[index]);
    float3 tangent = make_float3(sortedTangent[index]);

    // get address in grid
    int3 gridPos = calcGridPos(pos);

    // examine neighbouring cells
    float3 force = make_float3(0.0f);

    for (int z=-1; z<=1; z++)
    {
        for (int y=-1; y<=1; y++)
        {
            for (int x=-1; x<=1; x++)
            {
                int3 neighbourPos = gridPos + make_int3(x, y, z);
                force += collideCell(neighbourPos,
                                     index,
                                     pos,
                                     vel,
                                     tangent,
                                     sortedPos,
                                     sortedVel,
                                     sortedTangent,
                                     cellStart,
                                     cellEnd,
                                     gridParticleIndex);
            }
        }
    }

    // write new velocity back to original unsorted location
    uint originalIndex = gridParticleIndex[index];
    forces[originalIndex] += make_float4(force, 0.0f);
}


//
// Implementation of the GROMACS harmonic angle potential
// http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#harmonic-angle-potential
//
__global__
void filamentKernel(float4 *forces,     // update: particle forces
                    float4 *tangent,    // output: filament tangent at particle
                    float4 *pos,        // input:  particle positions
                    uint numFilaments)
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numFilaments) return;

    const uint size    = params.filamentSize;
    const float k_bend = params.bondBendingK;

    uint i = 0;
    uint head_i = index * size;

    float3 r_i  = make_float3(0.0f), r_j  = make_float3(0.0f), r_k  = make_float3(0.0f);
    float3 f_i  = make_float3(0.0f), f_j  = make_float3(0.0f), f_k  = make_float3(0.0f);
    float3 r_ij = make_float3(0.0f), r_jk = make_float3(0.0f), r_ik = make_float3(0.0f);
    float A = 0.0f, d_ij = 0.0f, d_jk = 0.0f, d_ij_jk = 0.0f;
    float3 B = make_float3(0.0f), C = make_float3(0.0f);

    // Bond bending forces
    for (int p = 0; p < size - 2; p++)
    {
        i    = head_i + p;
        r_i  = make_float3(pos[i]);
        r_j  = make_float3(pos[i + 1]);
        r_k  = make_float3(pos[i + 2]);

        r_ij = r_j - r_i;
        r_jk = r_k - r_j;
        d_ij = lengthPeriodic(r_ij);
        d_jk = lengthPeriodic(r_jk);
        d_ij_jk = dot(r_ij, r_jk);

        // Derivation 1
        // A = (k / (d_ij*d_ij*d_jk*d_jk));
        // f_i = -A * ( (d_ij_jk * d_ij_jk / (d_ij*d_ij)) * r_ij - r_jk );
        // f_k = -A * ( (d_ij_jk * d_ij_jk / (d_jk*d_jk)) * r_jk - r_ij );
        // f_j = -f_i - f_k;

        // Derivation 2
        A = k_bend / (d_ij * d_jk);
        B = r_ij * d_ij_jk / dot(r_ij, r_ij);
        C = r_jk * d_ij_jk / dot(r_jk, r_jk);

        f_i = -A * (r_jk - B);
        f_k = -A * (C - r_ij);
        f_j = -f_i - f_k;

        forces[i]     += make_float4(f_i, 0.0f);
        forces[i + 1] += make_float4(f_j, 0.0f);
        forces[i + 2] += make_float4(f_k, 0.0f);

        // Compute tangent at j
        tangent[i + 1] = make_float4(getTangent(r_i, r_k), 0.0f);

        if (p == 0) // Head tangent
        {
            tangent[i] = make_float4(getTangent(r_i, r_j), 0.0f);
        }
        else if (p == size - 3) // last iteration do tail as well
        {
            tangent[i + 2] = make_float4(getTangent(r_j, r_k), 0.0f);
        }
    }
}

__global__
void reverseFilamentsKernel(float4 *force,      // update: particle forces
                            float4 *forceOld,   // update: particle old forces
                            float4 *vel,        // update: particle velocities
                            float4 *pos,        // update: particle positions
                            float *uniform,     // input: random numbers drawn from uniform dist
                            uint numFilaments   // input: number of filaments
                            )
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numFilaments) return;
    if (uniform[index] > params.reverseProbability) return;

    const uint size = params.filamentSize;
    uint head = index * size;
    uint i,j;

    float4 tmp_r, tmp_v, tmp_f, tmp_fOld;

    for (uint p = 0; p < size / 2; p++)
    {
        i = head            + p;
        j = head + size - 1 - p;
        tmp_r       = pos[i];
        tmp_v       = vel[i];
        tmp_f       = force[i];
        tmp_fOld    = forceOld[i];
        pos[i]      = pos[j];
        vel[i]      = vel[j];
        force[i]    = force[j];
        forceOld[i] = forceOld[j];
        pos[j]      = tmp_r;
        vel[j]      = tmp_v;
        force[j]    = tmp_f;
        forceOld[j] = tmp_fOld;
    }
}

__global__
void langevinKernel(float4 *force,     // update: particle forces
                    float4 *vel,       // input: particle velocities
                    float4 *gaussian,  // input: random numbers drawn from normal dist
                    uint numParticles  // input: number of particles
                 )
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    float3 f = make_float3(0.0f);
    f -= params.gamma * make_float3(vel[index]);
    f += make_float3(gaussian[index]);
    force[index] += make_float4(f, 0.0f);
}

#endif
