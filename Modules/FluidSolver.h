//
// Created by Иван on 11/17/2024.
//

#ifndef MYFLIP2D_FLUIDSOLVER_H
#define MYFLIP2D_FLUIDSOLVER_H

#include "Utility.h"
#include <vector>
#include <cmath>
#include <cassert>


struct Particle{
    float x;
    float y;
    float u;
    float v;
    //TODO: add affine C for APIC
};

class FluidSolver {

    float density;

    //число ячеек
    int numX;       //=floor(width / spacing) + 1
    int numY;       //=floor(height / spacing) + 1
    float h;        //=max(width/numX, height/numY)
    float h_inv;     // = 1.0/h
    int numCells;   //= numX * numY;

    // жидкость (значения в ячейках)
    std::vector<float> u; //size = numCells
    std::vector<float> v; //size = numCellss
    std::vector<float> du;//size = numCellss
    std::vector<float> dv;//size = numCellss
    std::vector<float> u_prev;//size = numCellss
    std::vector<float> v_prev;//size = numCellss
    std::vector<float> p;//size = numCellss
    std::vector<float> s;//size = numCellss
    std::vector<int> cellType;//size = numCellss
    std::vector<float> pressure; // size=numCells
    //std::vector<float> densities; //size = numCells
    //std::vector<int> cellColor; //3*NumCells

    // частицы
    int maxParticles;
    std::vector<Particle> particles;
    std::vector<float> particlePos; //2*maxParticles; //x1,y1, x2,y2, x3,y3, .....
    //std::vector<float> particleColor; //3*maxParticles; //r1,g1,b1, r2,g2,b2, ....
    std::vector<float> particleVel; //2*maxParticles; //u1,v1, u2,v2, u3,v3, ...
    std::vector<float> particleDens; //n=NumCells
    float particleRestDensity; //=0.0

    float particleRadius; //
    float pInvSpacing; //1.0 / (2.2 * particleRadius)
    int pNumX; //floor(width * pInvSpacing) + 1
    int pNumY; //floor(height * pInvSpacing) + 1
    int pNumCells; //pNumX * pNumY

    std::vector<int> numCellParticles; //pNumCells
    std::vector<int> firstCellParticle; //pNumCells + 1
    std::vector<int> cellParticleIds; //maxParticles

    int numParticles; //=0

    void setupSimulation();
    void relabel();
    void integrateParticles(float dt, float g);
    void handleParticleCollisions();
    void pushParticlesApart(const int numIters);
    void transferVelocitiesToGrid();
    void updateParticleDensity();
    void pressureSolve(const float dt);
    void applyPrecon(std::vector<double>& z, std::vector<double>& r, std::vector<double>& precon, std::vector<double>& Adiag,  std::vector<double>& Ax, std::vector<double>& Ay);
    void applyPressure(const float& dt);
    void applyA(std::vector<double>& z,std::vector<double>& s, std::vector<double>& Adiag, std::vector<double>& Ax, std::vector<double>& Ay);

        bool isFluid(int i, int j);
    void applyBodyForces(const float dt, const float g);
    void makeIncompressible(const int numIters, const float dt, const float overRelaxation = 1.9);
    void transferVelocitiesToParticles(const float flipCoef);
    void runFrameSimulation(const float dt, const float g, const float flipCoef,
                            const int numPressureIters, const int numParticleIters);

public:
    FluidSolver(const float dens, const float width, const float height, const float h_step,
                const float particleRadius, const int maxParticles);
    void setUpParticlesAndCells(const int particleAmount, std::vector<float> particlePositions);
    void runSimulation(const float dt, const float g, const float flipCoef, const int numFrames,
                       const int numPressureIters, const int numParticleIters,
                       const std::string outputFileName);
};


#endif //MYFLIP2D_FLUIDSOLVER_H
