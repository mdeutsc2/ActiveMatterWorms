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

#ifndef PARTICLES_KERNEL_H
#define PARTICLES_KERNEL_H

#include "vector_types.h"
typedef unsigned int uint;

#define PBC_X 1
#define PBC_Y 0
#define PBC_Z 0

// simulation parameters
struct SimParams
{
    float3 gravity;
    float particleRadius;

    uint3 gridSize;
    uint numCells;
    float3 origin;
    float3 boxSize;
    float3 cellSize;

    uint numParticles;
    uint numFilaments;
    uint filamentSize;

    float bondSpringK;
    float bondSpringL;
    float bondBendingK;

    float activity;
    float reverseProbability;
    float kbT;
    float gamma;

    float spring;
    float damping;
    float shear;
    float attraction;
    float boundaryDamping;
};

#endif
