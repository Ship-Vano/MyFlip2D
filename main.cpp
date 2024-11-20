#include "Modules/FluidSolver.h"

int main() {
    float gravity = -9.81;
    float dt = 1.0 / 60.0;
    float flipCoef = 0.9;
    int numPressureIters = 50;
    int numParticleIters = 2;
    //float overRelaxation = 1.9;
    float simHeight = 3.0;
    float simWidth = 3.0;
    float resolution = 100.0;
    float hh = simHeight / resolution;
    float rho = 1000.0;
    float relWaterHeight = 0.8;
    float relWaterWidth = 0.8;
    float partRadius = 0.3 * hh;
    float dx = 2.0 * partRadius;
    float dy = std::sqrt(3.0)/2.0 * dx;

    int numParticlesX = std::floor( (relWaterWidth * simWidth - 2.0 * hh - 2.0 * partRadius) / dx );
    int numParticlesY = std::floor( (relWaterHeight * simHeight - 2.0 * hh - 2.0 * partRadius) / dy);
    int maxSimParticles = numParticlesX * numParticlesY;

    FluidSolver solver(rho, simWidth, simHeight, hh,
                       partRadius, maxSimParticles);

    std::vector<float> particlePositions(2*maxSimParticles, 0.0);
    int p = 0;
    for (int i = 0; i < numParticlesX; ++i) {
        for (int j = 0; j < numParticlesY; ++j) {
            particlePositions[p++] = hh + partRadius + dx * static_cast<float>(i) + (j % 2 == 0 ? 0.0 : partRadius);
            particlePositions[p++] = hh + partRadius + dy * static_cast<float>(j);
        }
    }

    solver.setUpParticlesAndCells(maxSimParticles, particlePositions);

    solver.runSimulation(dt, gravity, flipCoef, 10,
                        numPressureIters, numParticleIters, "./OutputData/res.txt");
    return 0;
}