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

// OpenGL Graphics includes
// #define HELPERGL_EXTERN_GL_FUNC_IMPLEMENTATION
// #include <helper_gl.h>

#include "particleSystem.h"
#include "particleSystem.cuh"
#include "particles_kernel.cuh"

#include <cuda_runtime.h>

#include <helper_functions.h>
#include <helper_math.h>
#include <helper_cuda.h>

#include <assert.h>
#include <math.h>
#include <memory.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <fstream>

#ifndef CUDART_PI_F
#define CUDART_PI_F         3.141592654f
#endif

ParticleSystem::ParticleSystem(uint numFilaments, uint filamentSize, uint3 gridSize, bool srdSolvent) :
    m_bInitialized(false),
    m_updateCount(0),
    m_numFilaments(numFilaments),
    m_filamentSize(filamentSize),
    m_numParticles(numFilaments * filamentSize),
    m_hPos(0),
    m_hVel(0),
    m_hForce(0),
    m_hvForce(0),
    m_hTangent(0),
    m_dPos(0),
    m_dVel(0),
    m_dTangent(0),
    m_dForce(0),
    m_dvForce(0),
    m_dForceOld(0),
    m_dvForceOld(0),
    m_dRandom(0),
    m_gridSize(gridSize),
    m_solventEnabled(srdSolvent),
    m_dSolventCellCOM(0),
    m_numSolvent(0),
    m_numSolventCells(0),
    m_solventUpdateCount(0),
    m_solventUpdateGap(0),
    m_timer(NULL)
{
    // Filament parameters
    m_params.numFilaments = m_numFilaments;
    m_params.filamentSize = m_filamentSize;
    m_params.numParticles = m_numParticles;
    m_params.particleRadius = 0.5f;
    m_params.particleMass = 1.0f;
    m_params.solSigma = 3.4f; // sigma for LJ r_cutoff interaction

    // System size/boundaries
    m_gridSortBits = 18;    // increase this for larger grids
    m_numGridCells = m_gridSize.x*m_gridSize.y*m_gridSize.z;
    m_params.gridSize = m_gridSize;
    m_params.numCells = m_numGridCells;
    float cellSize = m_params.particleRadius * 2.0f;  // cell size equal to particle diameter
    m_params.cellSize = make_float3(cellSize, cellSize, cellSize);
    m_params.boxSize = make_float3(gridSize.x*cellSize, gridSize.y*cellSize, gridSize.z*cellSize);
    m_params.origin = -m_params.boxSize / 2.0f;
    m_params.boundaryX = BoundaryType::PERIODIC;
    m_params.boundaryY = BoundaryType::PERIODIC;
    m_params.boundaryZ = BoundaryType::PERIODIC;

    // Particle-Particle bonding
    m_params.bondSpringK = 57.f;
    m_params.bondSpringL = m_params.particleRadius * 0.8f;
    m_params.bondBendingK = 0.001f;

    // Particle-Particle forces
    m_params.spring = 0.5f;
    m_params.damping = 0.02f;
    m_params.shear = 0.1f;
    m_params.attraction = 0.0f;
    m_params.boundaryDamping = -0.5f;

    // Active forces
    m_params.activity = 0.1f;
    m_params.reverseProbability = 0.0f;

    // Langevin thermostat
    m_params.kbT = 1.0f;
    m_params.gamma = 0.01;

    // Global interations
    m_params.gravity = make_float3(0.0f);

    // SRD Solvent particle parameters
    _initSolventParams();

    _initialize(m_numParticles, m_numSolvent);
}

ParticleSystem::~ParticleSystem()
{
    _finalize();
    m_numParticles = m_numFilaments = m_filamentSize = 0;
}

inline float frand()
{
    return rand() / (float) RAND_MAX;
}

void
ParticleSystem::_initSolventParams()
{
    if (m_solventEnabled)
    {
        // TODO: add solvent radius 
        m_solvent.enabled = m_solventEnabled;
        m_solvent.alpha = 2.26893f; // 130 deg TODO: make settings
        // m_solvent.gridSize = make_uint3(m_params.gridSize.x, m_params.gridSize.y, max(m_params.gridSize.z, 1));
        // m_solvent.cellSize = m_params.cellSize;
        m_solvent.gridSize = make_uint3(m_params.gridSize.x / 2, m_params.gridSize.y / 2, max(m_params.gridSize.z / 2, 1));
        m_solvent.cellSize = m_params.cellSize * 2;
        m_solvent.numCells = m_solvent.gridSize.x * m_solvent.gridSize.y * m_solvent.gridSize.z;
        m_solvent.randOffset = make_float3(0.f, 0.f, 0.f);
        m_solvent.particleMass = 1.0f;
        m_solvent.kbT = m_params.kbT;

        float density2D = 5.f; // TODO: make setting
        m_solventUpdateGap = 1; // TODO: make settings
        m_numSolvent =  uint(density2D * m_solvent.numCells);
        m_numSolventCells = m_solvent.numCells;
        m_solvent.numParticles = m_numSolvent;

        printf("SRD Solvent Enabled.\n");
    }
}

// create a color ramp
void colorRamp(float t, float *r)
{
    const int ncolors = 7;
    float c[ncolors][3] =
    {
        { 1.0, 0.0, 0.0, },
        { 1.0, 0.5, 0.0, },
        { 1.0, 1.0, 0.0, },
        { 0.0, 1.0, 0.0, },
        { 0.0, 1.0, 1.0, },
        { 0.0, 0.0, 1.0, },
        { 1.0, 0.0, 1.0, },
    };
    t = t * (ncolors-1);
    int i = (int) t;
    float u = t - floor(t);
    r[0] = lerp(c[i][0], c[i+1][0], u);
    r[1] = lerp(c[i][1], c[i+1][1], u);
    r[2] = lerp(c[i][2], c[i+1][2], u);
}

void
ParticleSystem::_initialize(int numParticles, int numSolvent)
{
    assert(!m_bInitialized);

    m_numParticles = numParticles;
    int N = numParticles + numSolvent;

    // allocate host storage
    m_hPos     = new float[N*4];
    m_hVel     = new float[N*4];
    m_hForce   = new float[m_numParticles*4];
    m_hvForce  = new float[numSolvent*4];
    m_hTangent = new float[m_numParticles*4];
    memset(m_hPos, 0, N*4*sizeof(float));
    memset(m_hVel, 0, N*4*sizeof(float));
    memset(m_hForce, 0, m_numParticles*4*sizeof(float));
    memset(m_hvForce, 0, numSolvent*4*sizeof(float));
    memset(m_hTangent, 0, m_numParticles*4*sizeof(float));

    m_hCellStart = new uint[m_numGridCells];
    memset(m_hCellStart, 0, m_numGridCells*sizeof(uint));

    m_hCellEnd = new uint[m_numGridCells];
    memset(m_hCellEnd, 0, m_numGridCells*sizeof(uint));

    // allocate GPU data
    unsigned int memSize = sizeof(float) * 4 * m_numParticles;
    unsigned int memSizeS = sizeof(float) * 4 * numSolvent;
    unsigned int memSizeN = sizeof(float) * 4 * N;
    allocateArray((void **)&m_dPos, memSizeN);
    allocateArray((void **)&m_dVel, memSizeN);
    allocateArray((void **)&m_dTangent, memSize);
    allocateArray((void **)&m_dForce, memSize);
    allocateArray((void **)&m_dvForce, memSizeS);
    allocateArray((void **)&m_dForceOld, memSize);
    allocateArray((void **)&m_dvForceOld, memSizeS);
    allocateArray((void **)&m_dRandom, memSize);
    allocateArray((void **)&m_dSortedPos, memSizeN);
    allocateArray((void **)&m_dSortedVel, memSizeN);
    allocateArray((void **)&m_dSortedTangent, memSize);
    allocateArray((void **)&m_dGridParticleHash, N*sizeof(uint));
    allocateArray((void **)&m_dGridParticleIndex, N*sizeof(uint));
    allocateArray((void **)&m_dCellStart, m_numGridCells*sizeof(uint));
    allocateArray((void **)&m_dCellEnd, m_numGridCells*sizeof(uint));

    if (m_solventEnabled)
    {
        allocateArray((void **)&m_dSolventCellCOM, m_solvent.numCells*4*sizeof(float));
        setSolventParameters(&m_solvent);
    }

    sdkCreateTimer(&m_timer);
    setParameters(&m_params);

    m_bInitialized = true;
}

void
ParticleSystem::_finalize()
{
    assert(m_bInitialized);

    delete [] m_hPos;
    delete [] m_hVel;
    delete [] m_hTangent;
    delete [] m_hCellStart;
    delete [] m_hCellEnd;

    freeArray(m_dPos);
    freeArray(m_dVel);
    freeArray(m_dTangent);
    freeArray(m_dSortedPos);
    freeArray(m_dSortedVel);

    freeArray(m_dGridParticleHash);
    freeArray(m_dGridParticleIndex);
    freeArray(m_dCellStart);
    freeArray(m_dCellEnd);
}

// step the simulation
void
ParticleSystem::update(float deltaTime)
{
    assert(m_bInitialized);

    // update constants
    setParameters(&m_params);

    integrateFilaments(
        m_dPos,
        m_dVel,
        m_dForce,
        m_dForceOld,
        deltaTime,
        m_numParticles);

    if (m_params.reverseProbability > 0.0f)
    {
        reverseFilaments(
            m_dForce,
            m_dForceOld,
            m_dVel,
            m_dPos,
            m_dRandom,
            m_numFilaments);
    }

    if (!m_solventEnabled)
    {
        langevinThermostat(
            m_dForce,
            m_dVel,
            m_dRandom,
            m_params.gamma,
            m_params.kbT,
            m_numParticles);
    }

    // forces of filament bonds
    filamentForces(
        m_dForce,
        m_dTangent,
        m_dPos,
        m_numFilaments);

    // calculate grid hash for filament particles only
    calcHash(
        m_dGridParticleHash,
        m_dGridParticleIndex,
        m_dPos,
        m_numParticles);

    // sort filament particles based on hash
    sortParticles(m_dGridParticleHash, m_dGridParticleIndex, m_numParticles);

    // reorder filament particles into sorted order and
    // find start and end of each cell
    reorderDataAndFindCellStart(
        m_dCellStart,
        m_dCellEnd,
        m_dSortedPos,
        m_dSortedVel,
        m_dSortedTangent,
        m_dGridParticleHash,
        m_dGridParticleIndex,
        m_dPos,
        m_dVel,
        m_dTangent,
        m_numParticles,
        m_numGridCells);

    // MD collision forces
    collideFilaments(
        m_dForce,
        m_dSortedPos,
        m_dSortedVel,
        m_dSortedTangent,
        m_dGridParticleIndex,
        m_dCellStart,
        m_dCellEnd,
        m_numParticles,
        m_numGridCells);

    bool doSolventUpdate = m_solventEnabled && (m_updateCount % m_solventUpdateGap) == 0;
    if (doSolventUpdate)
    {
        _solventUpdate(deltaTime * m_solventUpdateGap);
    }

    m_updateCount++;
}

void
ParticleSystem::_solventUpdate(float dt)
{
    // random solvent cell offsets
    m_solvent.randOffset.x = frand() * m_solvent.cellSize.x;
    m_solvent.randOffset.y = frand() * m_solvent.cellSize.y;
    // m_solvent.randOffset.z = randf() * m_solvent.cellSize.z;

    setSolventParameters(&m_solvent);

    // calculate grid hash including solvent particles
    calcHash(
        m_dGridParticleHash,
        m_dGridParticleIndex,
        m_dPos,
        m_numParticles + m_numSolvent,
        true);

    // sort particles based on hash including solvent particles
    sortParticles(m_dGridParticleHash, m_dGridParticleIndex, m_numParticles + m_numSolvent);

    // reorder filament and solvent particles together
    reorderDataAndFindCellStart(
        m_dCellStart,
        m_dCellEnd,
        m_dSortedPos,
        m_dSortedVel,
        NULL,
        m_dGridParticleHash,
        m_dGridParticleIndex,
        m_dPos,
        m_dVel,
        NULL,
        m_numParticles + m_numSolvent,
        m_numGridCells);

    solventIsokineticThermostat(
        m_dVel,
        m_dSortedVel,
        m_dCellStart,
        m_dCellEnd,
        m_dGridParticleIndex,
        m_numSolventCells
    );

    // solventCellCentorOfMomentum(
    //     m_dSolventCellCOM,
    //     m_dSortedVel,
    //     m_dCellStart,
    //     m_dCellEnd,
    //     m_dGridParticleIndex,
    //     m_numSolventCells
    // );

    collideSolvent(
        m_dPos,
        m_dVel,
        m_dvForce,
        m_dCellStart,
        m_dCellEnd,
        m_dGridParticleIndex,
        m_numParticles + m_numSolvent
    );

    integrateSolvent(
        m_dPos,
        m_dVel,
        m_dvForce,
        m_dvForceOld,
        dt,
        m_numParticles,
        m_numSolvent
    );

    m_solventUpdateCount++;
}

void
ParticleSystem::dumpGrid()
{
    // dump grid information
    copyArrayFromDevice(m_hCellStart, m_dCellStart, 0, sizeof(uint)*m_numGridCells);
    copyArrayFromDevice(m_hCellEnd, m_dCellEnd, 0, sizeof(uint)*m_numGridCells);
    uint maxCellSize = 0;

    for (uint i=0; i<m_numGridCells; i++)
    {
        if (m_hCellStart[i] != 0xffffffff)
        {
            uint cellSize = m_hCellEnd[i] - m_hCellStart[i];
            if (cellSize > maxCellSize)
            {
                maxCellSize = cellSize;
            }
        }
    }

    printf("maximum particles per cell = %d\n", maxCellSize);
}


void
ParticleSystem::dumpSolventGrid()
{
    // dump grid information

    copyArrayFromDevice(m_hCellStart, m_dCellStart, 0, sizeof(uint)*m_numSolventCells);
    copyArrayFromDevice(m_hCellEnd, m_dCellEnd, 0, sizeof(uint)*m_numSolventCells);

    float *cellCOM = new float[4*m_numSolventCells];
    copyArrayFromDevice(cellCOM, m_dSolventCellCOM, 0, 4*sizeof(float)*m_numSolventCells);

    uint maxCellSize = 0;

    for (uint i = 0; i < m_numSolventCells; i++)
    {
        if (m_hCellStart[i] != 0xffffffff)
        {
            uint cellSize = m_hCellEnd[i] - m_hCellStart[i];
            if (cellSize > maxCellSize)
            {
                printf("[%i] %i CoM = {%f %f} \n", i, cellSize, cellCOM[i*4], cellCOM[i*4 + 1]);
                maxCellSize = cellSize;
            }
        }
    }

    printf("maximum particles per solvent cell = %d\n", maxCellSize);
}

void
ParticleSystem::dumpParticles(uint start, uint count)
{
    copyArrayFromDevice(m_hVel, m_dPos, 0, sizeof(float)*4*count);
    copyArrayFromDevice(m_hVel, m_dVel, 0, sizeof(float)*4*count);

    for (uint i=start; i<start+count; i++)
    {
        printf("pos: (%.4f, %.4f, %.4f, %.4f)\n", m_hPos[i*4+0], m_hPos[i*4+1], m_hPos[i*4+2], m_hPos[i*4+3]);
        printf("vel: (%.4f, %.4f, %.4f, %.4f)\n", m_hVel[i*4+0], m_hVel[i*4+1], m_hVel[i*4+2], m_hVel[i*4+3]);
    }
}

void
ParticleSystem::writeOutputs(const std::string& fileName, int iteration, float timeDelta)
{
    const int N = m_numParticles + m_numSolvent;
    copyArrayFromDevice(m_hPos, m_dPos, 0, sizeof(float)*4*N);
    copyArrayFromDevice(m_hVel, m_dVel, 0, sizeof(float)*4*N);
    
    // Write SimParams as json object
    auto p = this->m_params;
    if (iteration == 0) {
        std::ofstream fout_params;
        fout_params.open("parameters.json",std::ios::out | std::ios::app);
        fout_params << "{"
         << "\"iteration\":"                 << iteration << ","
         << "\"dt\":"                        << timeDelta << ","
         << "\"time\":"                      << iteration * timeDelta << ","
         << "\"num_particles\":"             << p.numParticles << ","
         << "\"num_filaments\":"             << p.numFilaments << ","
         << "\"filament_size\":"             << p.filamentSize << ","
         << "\"bond_spring_coef\":"          << p.bondSpringK << ","
         << "\"bond_spring_length\":"        << p.bondSpringL << ","
         << "\"bond_bending_coef\":"         << p.bondBendingK << ","
         << "\"activity\":"                  << p.activity << ","
         << "\"kbT\":"                       << p.kbT << ","
         << "\"gamma\":"                     << p.gamma << ","
         << "\"particle_radius\":"           << p.particleRadius << ","
         << "\"particle_mass\":"             << p.particleMass << ","
         << "\"num_cells\":"                 << p.numCells << ","
         << "\"grid_size\":"          << "[" << p.gridSize.x << "," << p.gridSize.y << "," << p.gridSize.z << "],"
         << "\"box_size\":"           << "[" << p.boxSize.x << "," << p.boxSize.y << "," << p.boxSize.z << "],"
         << "\"cell_size\":"          << "[" << p.cellSize.x << "," << p.cellSize.y << "," << p.cellSize.z << "],"
         << "\"boundaries\":"         << "[" << p.boundaryX << "," << p.boundaryY << "," << p.boundaryZ << "],"
         << "\"reverse_probability\":"       << p.reverseProbability << ","
         << "\"collision_spring_coef\":"     << p.spring << ","
         << "\"collision_damping_coef\":"    << p.damping << ","
         << "\"collision_shear_coef\":"      << p.shear << ","
         << "\"collision_attraction_coef\":" << p.attraction << ","
         << "\"boundary_damping_coef\":"     << p.boundaryDamping
         << "}" << std::endl;
         fout_params.close();
    }

    std::ofstream fout;
    fout.open(fileName, std::ios::out | std::ios::app);
    fout << N << std::endl;


    fout << "# iteration "<< iteration << ", dt "<< timeDelta << std::endl;

    const char * type = "A";
    for (int i = 0; i < N; i++)
    {
        type = i < m_numParticles ? "A" : "B";
        fout << type << " " << m_hPos[i*4] << " " << m_hPos[i*4+1] << " " << m_hPos[i*4+2]
                     << std::endl;
    }

    fout.close();

    // writing out a velocity .vtk file
    std::string insertedString = "_vel";
    std::string fileType = ".xyz";
    std::string vel_fileName;
    std::stringstream sstm;
    char delimiter = '.';

    std::size_t pos = fileName.find(delimiter);
    std::string prefix = fileName.substr(0, pos);
    sstm << prefix << insertedString << iteration << fileType;
    vel_fileName = sstm.str();
    // std::cout << iteration << "\n";
    // std::cout << vel_fileName << "\n";

    std::ofstream fout_vel;
    fout_vel.open(vel_fileName, std::ios::out | std::ios::app);
    fout_vel << "# iteration "<< iteration << ", dt "<< timeDelta << std::endl;
    //const char * type = "A";
    for (int i = 0; i < N; i++)
    {
        type = i < m_numParticles ? "A" : "B";
        fout_vel << type << " " << m_hVel[i*4] << " " << m_hVel[i*4+1] << " " << m_hVel[i*4+2]
                    << std::endl;
    }

    fout_vel.close();
}


void
ParticleSystem::writeOutputs_xyzv(const std::string& fileName, int iteration, float timeDelta)
{
    const int N = m_numParticles + m_numSolvent;
    copyArrayFromDevice(m_hPos, m_dPos, 0, sizeof(float)*4*N);
    copyArrayFromDevice(m_hVel, m_dVel, 0, sizeof(float)*4*N);

    std::ofstream fout;
    fout.open(fileName, std::ios::out | std::ios::app);
    fout << N << std::endl;

    // Write SimParams as json object
    auto p = this->m_params;
    fout << "{"
         << "\"iteration\":"                 << iteration << ","
         << "\"dt\":"                        << timeDelta << ","
         << "\"time\":"                      << iteration * timeDelta << ","
         << "\"num_particles\":"             << p.numParticles << ","
         << "\"num_filaments\":"             << p.numFilaments << ","
         << "\"filament_size\":"             << p.filamentSize << ","
         << "\"bond_spring_coef\":"          << p.bondSpringK << ","
         << "\"bond_spring_length\":"        << p.bondSpringL << ","
         << "\"bond_bending_coef\":"         << p.bondBendingK << ","
         << "\"activity\":"                  << p.activity << ","
         << "\"kbT\":"                       << p.kbT << ","
         << "\"gamma\":"                     << p.gamma << ","
         << "\"particle_radius\":"           << p.particleRadius << ","
         << "\"particle_mass\":"             << p.particleMass << ","
         << "\"num_cells\":"                 << p.numCells << ","
         << "\"grid_size\":"          << "[" << p.gridSize.x << "," << p.gridSize.y << "," << p.gridSize.z << "],"
         << "\"box_size\":"           << "[" << p.boxSize.x << "," << p.boxSize.y << "," << p.boxSize.z << "],"
         << "\"cell_size\":"          << "[" << p.cellSize.x << "," << p.cellSize.y << "," << p.cellSize.z << "],"
         << "\"boundaries\":"         << "[" << p.boundaryX << "," << p.boundaryY << "," << p.boundaryZ << "],"
         << "\"reverse_probability\":"       << p.reverseProbability << ","
         << "\"collision_spring_coef\":"     << p.spring << ","
         << "\"collision_damping_coef\":"    << p.damping << ","
         << "\"collision_shear_coef\":"      << p.shear << ","
         << "\"collision_attraction_coef\":" << p.attraction << ","
         << "\"boundary_damping_coef\":"     << p.boundaryDamping
         << "}" << std::endl;

    const char * type = "A";
    for (int i = 0; i < N; i++)
    {
        type = i < m_numParticles ? "A" : "B";
        fout << type << " " << m_hPos[i*4] << " " << m_hPos[i*4+1] << " " << m_hPos[i*4+2]
                     << " " << m_hVel[i*4] << " " << m_hVel[i*4+1] << " " << m_hVel[i*4+2]
                     << std::endl;
    }

    fout.close();
}

float *
ParticleSystem::getArray(ParticleArray array)
{
    assert(m_bInitialized);

    float *hdata = 0;
    float *ddata = 0;
    struct cudaGraphicsResource *cuda_vbo_resource = 0;

    switch (array)
    {
        default:
        case POSITION:
            hdata = m_hPos;
            ddata = m_dPos;
            break;

        case VELOCITY:
            hdata = m_hVel;
            ddata = m_dVel;
            break;
    }

    copyArrayFromDevice(hdata, ddata, &cuda_vbo_resource, m_numParticles*4*sizeof(float));
    return hdata;
}

void
ParticleSystem::setArray(ParticleArray array, const float *data, int start, int count)
{
    assert(m_bInitialized);

    switch (array)
    {
        default:
        case POSITION:
            copyArrayToDevice(m_dPos, data, start*4*sizeof(float), count*4*sizeof(float));
            break;

        case VELOCITY:
            copyArrayToDevice(m_dVel, data, start*4*sizeof(float), count*4*sizeof(float));
            break;

        case TANGENT:
            copyArrayToDevice(m_dTangent, data, start*4*sizeof(float), count*4*sizeof(float));
            break;
    }
}

void
ParticleSystem::reset()
{
    float filamentLength = m_params.filamentSize * m_params.bondSpringL + 2.0f * m_params.particleRadius;
    int xdim_max = int(floor(m_params.boxSize.x / filamentLength));
    int ydim_max = int(floor(m_params.boxSize.y / (2 * m_params.particleRadius)));
    int xyplane_max = xdim_max * ydim_max;
    float x,y,z,x_head;
    uint w, p = 0, v = 0, t = 0;

    printf("Placing filament particles.\n");
    printf("filamentLength: %.4f\n", filamentLength);
    printf("filamentsXdim:  %i\n",   xdim_max);
    printf("filamentsYdim:  %i\n",   ydim_max);

    for (uint i = 0; i < m_numParticles; i++)
    {
        w = (i / m_params.filamentSize);
        x_head = (w % xdim_max) * filamentLength + m_params.particleRadius;
        x = x_head + (i % m_params.filamentSize) * m_params.bondSpringL;
        y = 2.f * m_params.particleRadius * ((w / xdim_max) % ydim_max + 1);
        z = 2.f * m_params.particleRadius * (w / xyplane_max + 1);
        m_hPos[p++] = x + m_params.origin.x;
        m_hPos[p++] = y + m_params.origin.y;
        m_hPos[p++] = z + m_params.origin.z;
        m_hPos[p++] = m_params.particleRadius; // radius
        m_hVel[v++] = 0.1f * (frand()*2.0f-1.0f);
        m_hVel[v++] = 0.1f * (frand()*2.0f-1.0f);
        m_hVel[v++] = 0.1f * (frand()*2.0f-1.0f);
        m_hVel[v++] = m_params.particleMass;
        m_hTangent[t++] = 1.0f;
        m_hTangent[t++] = m_hTangent[t++] = m_hTangent[t++] = 0.0f;
    }

    // Flip half the filaments
    for (w = 0; w < m_params.numFilaments; w++)
    {
        if (frand() < 0.5f) continue;

        uint i,j;
        float r;
        for (p = 0; p < m_params.filamentSize / 2; p++)
        {
            i = w      * m_params.filamentSize      + p; // starts at head
            j = ((w+1) * m_params.filamentSize) - 1 - p; // starts at tail
            x = m_hPos[i*4 + 0]; // tmp save head pos
            y = m_hPos[i*4 + 1];
            z = m_hPos[i*4 + 2];
            r = m_hPos[i*4 + 3];
            m_hPos[i*4 + 0] = m_hPos[j*4 + 0]; // Assign j pos to i
            m_hPos[i*4 + 1] = m_hPos[j*4 + 1];
            m_hPos[i*4 + 2] = m_hPos[j*4 + 2];
            m_hPos[i*4 + 3] = m_hPos[j*4 + 3];
            m_hPos[j*4 + 0] = x; // Assign i pos to j
            m_hPos[j*4 + 1] = y;
            m_hPos[j*4 + 2] = z;
            m_hPos[j*4 + 3] = r;
            m_hTangent[i*4 + 0] = -m_hTangent[i*4 + 0]; // Invert tangent i
            m_hTangent[i*4 + 1] = -m_hTangent[i*4 + 1];
            m_hTangent[i*4 + 2] = -m_hTangent[i*4 + 2];
            m_hTangent[i*4 + 3] = -m_hTangent[i*4 + 3];
            m_hTangent[j*4 + 0] = -m_hTangent[j*4 + 0]; // Invert tangent j
            m_hTangent[j*4 + 1] = -m_hTangent[j*4 + 1];
            m_hTangent[j*4 + 2] = -m_hTangent[j*4 + 2];
            m_hTangent[j*4 + 3] = -m_hTangent[j*4 + 3];
        }
    }

    uint N = m_numParticles;

    if (m_solventEnabled)
    {
        printf("Placing solvent particles.\n");
        printf("numSolvent: %i\n", m_numSolvent);
        N += m_numSolvent;

        uint i;
        const float v0 = 5.0f; // TODO: relate to kbT
        for (int s = 0; s < m_numSolvent; s++)
        {
            i = m_numParticles + s;
            m_hPos[i*4 + 0] = m_params.boxSize.x * frand() + m_params.origin.x;
            m_hPos[i*4 + 1] = m_params.boxSize.y * frand() + m_params.origin.y;
            m_hPos[i*4 + 2] = m_params.boxSize.z * frand() + m_params.origin.z;
            m_hPos[i*4 + 3] = 0.0; // solvent radius not defined
            m_hVel[i*4 + 0] = v0 * (frand() * 2.0f - 1.0f);
            m_hVel[i*4 + 1] = v0 * (frand() * 2.0f - 1.0f);
            m_hVel[i*4 + 2] = v0 * (frand() * 2.0f - 1.0f);
            m_hVel[i*4 + 3] = m_solvent.particleMass;
        }
    }

    setArray(POSITION, m_hPos, 0, N);
    setArray(VELOCITY, m_hVel, 0, N);
    setArray(TANGENT, m_hTangent, 0, m_numParticles);
}
