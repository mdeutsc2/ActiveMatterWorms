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

// This file contains C wrappers around the some of the CUDA API and the
// kernel functions so that they can be called from "particleSystem.cpp"

#if defined(__APPLE__) || defined(MACOSX)
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#include <cstdlib>
#include <cstdio>
#include <string.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>

#include <helper_cuda.h>
#include <helper_functions.h>
#include "thrust/device_ptr.h"
#include "thrust/for_each.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/sort.h"

#include "particles_kernel_impl.cuh"

extern "C"
{
    curandGenerator_t rng;

    void cudaInit(int argc, char **argv)
    {
        int devID;

        // use command-line specified CUDA device, otherwise use device with highest Gflops/s
        devID = findCudaDevice(argc, (const char **)argv);

        if (devID < 0)
        {
            printf("No CUDA Capable devices found, exiting...\n");
            exit(EXIT_SUCCESS);
        }

        curandCreateGenerator(&rng, CURAND_RNG_PSEUDO_DEFAULT);
    }

    void randomizeUniform(float *device, uint n)
    {
        curandGenerateUniform(rng, device, n);
    }

    void randomizeGaussian(float *device, uint n, float mean, float stddev)
    {
        curandGenerateNormal(rng, device, n, mean, stddev);
    }

    void allocateArray(void **devPtr, size_t size)
    {
        checkCudaErrors(cudaMalloc(devPtr, size));
    }

    void freeArray(void *devPtr)
    {
        checkCudaErrors(cudaFree(devPtr));
    }

    void threadSync()
    {
        checkCudaErrors(cudaDeviceSynchronize());
    }

    void copyArrayToDevice(void *device, const void *host, int offset, int size)
    {
        checkCudaErrors(cudaMemcpy((char *) device + offset, host, size, cudaMemcpyHostToDevice));
    }

    void copyArrayFromDevice(void *host, const void *device,
                             struct cudaGraphicsResource **cuda_vbo_resource, int size)
    {
        checkCudaErrors(cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost));
    }

    void setParameters(SimParams *hostParams)
    {
        // copy parameters to constant memory
        checkCudaErrors(cudaMemcpyToSymbol(params, hostParams, sizeof(SimParams)));
    }

    void setSolventParameters(SolventParams *solventParams)
    {
        // copy parameters to constant memory
        checkCudaErrors(cudaMemcpyToSymbol(srd, solventParams, sizeof(SolventParams)));
    }

    //Round a / b to nearest higher integer value
    uint iDivUp(uint a, uint b)
    {
        return (a % b != 0) ? (a / b + 1) : (a / b);
    }

    // compute grid and thread block size for a given number of elements
    void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
    {
        numThreads = min(blockSize, n);
        numBlocks = iDivUp(n, numThreads);
    }

    void integrateFilaments(float *pos,
                            float *vel,
                            float *f,
                            float *f_old,
                            float deltaTime,
                            uint numParticles)
    {
        thrust::device_ptr<float4> d_pos4((float4 *)pos);
        thrust::device_ptr<float4> d_vel4((float4 *)vel);
        thrust::device_ptr<float4> d_f4((float4 *)f);
        thrust::device_ptr<float4> d_f_old4((float4 *)f_old);

        auto start = thrust::make_tuple(d_pos4, d_vel4, d_f4, d_f_old4);
        auto end = thrust::make_tuple(d_pos4+numParticles, d_vel4+numParticles, d_f4+numParticles, d_f_old4+numParticles);

        thrust::for_each(thrust::make_zip_iterator(start),thrust::make_zip_iterator(end), filament_integrator(deltaTime));
    }


    void integrateSolvent(float *pos,
                          float *vel,
                          float deltaTime,
                          uint startIndex,
                          uint numSolvent)
    {
        thrust::device_ptr<float4> d_pos4((float4 *)pos);
        thrust::device_ptr<float4> d_vel4((float4 *)vel);

        auto start = thrust::make_tuple(d_pos4 + startIndex, d_vel4 + startIndex);
        auto end = thrust::make_tuple(d_pos4 + startIndex + numSolvent, d_vel4 + startIndex + numSolvent);

        thrust::for_each(thrust::make_zip_iterator(start),thrust::make_zip_iterator(end), solvent_integrator(deltaTime));
    }

    void calcHash(uint  *gridParticleHash,
                  uint  *gridParticleIndex,
                  float *pos,
                  int    numParticles,
                  bool   srd)
    {
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // execute the kernel
        calcHashD<<< numBlocks, numThreads >>>(gridParticleHash,
                                               gridParticleIndex,
                                               (float4 *) pos,
                                               numParticles,
                                               srd);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed");
    }

    void reorderDataAndFindCellStart(uint  *cellStart,
                                     uint  *cellEnd,
                                     float *sortedPos,
                                     float *sortedVel,
                                     float *sortedTangent,
                                     uint  *gridParticleHash,
                                     uint  *gridParticleIndex,
                                     float *pos,
                                     float *vel,
                                     float *tangent,
                                     uint   numParticles,
                                     uint   numCells)
    {
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // set all cells to empty
        checkCudaErrors(cudaMemset(cellStart, 0xffffffff, numCells*sizeof(uint)));

        uint smemSize = sizeof(uint)*(numThreads+1);
        reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(
            cellStart,
            cellEnd,
            (float4 *) sortedPos,
            (float4 *) sortedVel,
            (float4 *) sortedTangent,
            gridParticleHash,
            gridParticleIndex,
            (float4 *) pos,
            (float4 *) vel,
            (float4 *) tangent,
            numParticles);

        getLastCudaError("Kernel execution failed: reorderDataAndFindCellStartD");
    }

    void collideFilaments(float *force,
                          float *sortedPos,
                          float *sortedVel,
                          float *sortedTangent,
                          uint  *gridParticleIndex,
                          uint  *cellStart,
                          uint  *cellEnd,
                          uint   numParticles,
                          uint   numCells)
    {

        // thread per particle
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // execute the kernel
        collideFilamentsKernel<<< numBlocks, numThreads >>>((float4 *)force,
                                                            (float4 *)sortedPos,
                                                            (float4 *)sortedVel,
                                                            (float4 *)sortedTangent,
                                                            gridParticleIndex,
                                                            cellStart,
                                                            cellEnd,
                                                            numParticles);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed: collide");

    }

    void collideSolvent(float *pos,
                        float *vel,
                        float *vforces,
                        uint *cellStart,
                        uint *cellEnd,
                        uint   numSolventCells,
                        uint   numParticles)
    {
        //randomizeUniform(random, numSolventCells);

        // thread per particle
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // execute the kernel
        collideSolventKernel<<< numBlocks, numThreads >>>((float4 *)pos,
                                                          (float4 *)vel,
                                                          (float4 *)vforces,
                                                          cellStart,
                                                          cellEnd,
                                                          numParticles);
    }

    void filamentForces(float *force,
                        float *tangent,
                        float *pos,
                        float numFilaments)
    {
        // thread per filament
        uint numThreads, numBlocks;
        computeGridSize(numFilaments, 256, numBlocks, numThreads);

        filamentKernel<<< numBlocks, numThreads >>>((float4 *)force,
                                                    (float4 *)tangent,
                                                    (float4 *)pos,
                                                    numFilaments);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed: filamentForces");
    }


    void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles)
    {
        thrust::sort_by_key(thrust::device_ptr<uint>(dGridParticleHash),
                            thrust::device_ptr<uint>(dGridParticleHash + numParticles),
                            thrust::device_ptr<uint>(dGridParticleIndex));
    }

    void reverseFilaments(float *force,
                          float *forceOld,
                          float *vel,
                          float *pos,
                          float *random,
                          uint numFilaments)
    {
        randomizeUniform(random, numFilaments);

        // thread per filament
        uint numThreads, numBlocks;
        computeGridSize(numFilaments, 256, numBlocks, numThreads);

        reverseFilamentsKernel<<< numBlocks, numThreads >>>((float4 *)force,
                                                            (float4 *)forceOld,
                                                            (float4 *)vel,
                                                            (float4 *)pos,
                                                            random,
                                                            numFilaments);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed: reverseFilaments");
    }

    void langevinThermostat(float *force,
                            float *vel,
                            float *random,
                            float gamma,
                            float kbT,
                            uint numParticles)
    {
        const float stddev = sqrtf(2.0f * gamma * kbT);
        randomizeGaussian(random, 4 * numParticles, 0.0f, stddev);

        // thread per particle
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        langevinKernel<<< numBlocks, numThreads >>>((float4 *)force,
                                                    (float4 *)vel,
                                                    (float4 *)random,
                                                    numParticles);

        getLastCudaError("Kernel execution failed: langevinThermostat");
    }


    void solventCellCentorOfMomentum(float *solventCellCOM,
                                     float *sortedVel,
                                     uint  *cellStart,
                                     uint  *cellEnd,
                                     uint  *gridParticleIndex,
                                     uint   numSolventCells)
    {
        uint numThreads, numBlocks;
        computeGridSize(numSolventCells, 256, numBlocks, numThreads);

        cellCenterOfMomentumSRD<<< numBlocks, numThreads >>>((float4 *)solventCellCOM,
                                                             (float4 *)sortedVel,
                                                             cellStart,
                                                             cellEnd,
                                                             gridParticleIndex,
                                                             numSolventCells);

        getLastCudaError("Kernel execution failed: solventCellCentorOfMomentum");
    }

    void solventIsokineticThermostat(float *vel,
                                     float *sortedVel,
                                     uint  *cellStart,
                                     uint  *cellEnd,
                                     uint  *gridParticleIndex,
                                     uint   numSolventCells)
    {
        uint numThreads, numBlocks;
        computeGridSize(numSolventCells, 256, numBlocks, numThreads);

        isokineticThermostatSRD<<< numBlocks, numThreads >>>((float4 *)vel,
                                                             (float4 *)sortedVel,
                                                             cellStart,
                                                             cellEnd,
                                                             gridParticleIndex,
                                                             numSolventCells);

        getLastCudaError("Kernel execution failed: solventCellCentorOfMomentum");
    }

}   // extern "C"
