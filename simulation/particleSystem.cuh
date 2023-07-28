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

extern "C"
{
    void cudaInit(int argc, char **argv);
    void randomizeUniform(float *device, uint n);
    void randomizeGaussian(float *device, uint n, float mean, float stddev);

    void allocateArray(void **devPtr, int size);
    void freeArray(void *devPtr);

    void threadSync();

    void copyArrayFromDevice(void *host, const void *device, struct cudaGraphicsResource **cuda_vbo_resource, int size);
    void copyArrayToDevice(void *device, const void *host, int offset, int size);

    void setParameters(SimParams *hostParams);
    void setSolventParameters(SolventParams *solventParams);

    void integrateFilaments(float *pos,
                            float *vel,
                            float *f,
                            float *f_old,
                            float deltaTime,
                            uint numParticles);

    void integrateSolvent(float *pos,
                          float *vel,
                          float *f,
                          float *f_old,
                          float deltaTime,
                          uint startIndex,
                          uint numSolvent);

    void calcHash(uint  *gridParticleHash,
                  uint  *gridParticleIndex,
                  float *pos,
                  int    numParticles,
                  bool   srd = false);

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
                                     uint   numCells);

    void collideFilaments(float *force,
                          float *sortedPos,
                          float *sortedVel,
                          float *sortedTangent,
                          uint  *gridParticleIndex,
                          uint  *cellStart,
                          uint  *cellEnd,
                          uint   numParticles,
                          uint   numCells);

    void collideSolvent(float *pos,
                        float *vel,
                        float *vforces,
                        uint *cellStart,
                        uint *cellEnd,
                        uint *gridParticleIndex,
                        uint   numParticles);

    void filamentForces(float *force,
                        float *tangent,
                        float *pos,
                        float numFilaments);

    void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles);

    void reverseFilaments(float *force,
                          float *forceOld,
                          float *vel,
                          float *pos,
                          float *uniform,
                          uint numFilaments);

    void langevinThermostat(float *force,
                            float *vel,
                            float *random,
                            float gamma,
                            float kbT,
                            uint numParticles);

    void solventCellCentorOfMomentum(float *solventCellCOM,
                                     float *sortedVel,
                                     uint  *cellStart,
                                     uint  *cellEnd,
                                     uint  *gridParticleIndex,
                                     uint   numCells);

    void solventIsokineticThermostat(float *vel,
                                     float *sortedVel,
                                     uint  *cellStart,
                                     uint  *cellEnd,
                                     uint  *gridParticleIndex,
                                     uint   numCells);
}
