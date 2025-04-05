#include "Modules/FluidSolver.h"

int main() {
    int task_type = 0;
    float gravity = -5.81f;
    float dt = 1.0f / 60.0f;
    float flipCoef = 0.5f;
    int numPressureIters = 100;
    int numParticleIters = 2;
    //float overRelaxation = 1.9;
    float simHeight = 720.0f;
    float simWidth = 720.0f;
    float resolution = 100.0f;
    float hh = simHeight / resolution;
    float rho = 1000.0f;
    float relWaterHeight = 0.8f;
    float relWaterWidth = 0.6f;
    float partRadius = 0.3f * hh;
    float dx = 2.0f * partRadius;
    float dy = dx;
    //float dy = std::sqrt(3.0f)/2.0f * dx;

    if(task_type == 0){
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

        solver.runSimulation(dt, gravity, flipCoef, 2000,
                             numPressureIters, numParticleIters, "OutputData/res.txt");
    }
    else if(task_type == 1){
        std::vector<float> particlePositions; // Пустой вектор для позиций

        // Нижний слой (40% высоты)
        float lowerWaterHeight = 0.6f * simHeight;
        int numParticlesX_lower = std::floor((relWaterWidth * simWidth - 2.0 * hh - 2.0 * partRadius) / dx);
        int numParticlesY_lower = std::floor((lowerWaterHeight - 2.0 * hh - 2.0 * partRadius) / dy);

        // Генерация частиц нижнего слоя
        for (int i = 0; i < numParticlesX_lower; ++i) {
            for (int j = 0; j < numParticlesY_lower; ++j) {
                float x = hh + partRadius + dx * i + (j % 2 == 0 ? 0.0 : partRadius);
                float y = hh + partRadius + dy * j;
                particlePositions.push_back(x);
                particlePositions.push_back(y);
            }
        }

        // Шар сверху
        float ballRadius = 0.1 * simHeight; // Радиус 10% от высоты
        float ballCenterX = simWidth / 2.0f;
        float ballCenterY = simHeight * 0.7f; // Центр на 70% высоты

        // Границы квадрата, охватывающего шар
        float ballLeft = ballCenterX - ballRadius;
        float ballRight = ballCenterX + ballRadius;
        float ballBottom = ballCenterY - ballRadius;
        float ballTop = ballCenterY + ballRadius;

        // Количество частиц в квадрате
        int numBallX = std::floor((ballRight - ballLeft - 2.0 * partRadius) / dx);
        int numBallY = std::floor((ballTop - ballBottom - 2.0 * partRadius) / dy);

        // Генерация частиц шара
        for (int i = 0; i < numBallX; ++i) {
            for (int j = 0; j < numBallY; ++j) {
                float x = ballLeft + partRadius + dx * i + (j % 2 == 0 ? 0.0f : partRadius);
                float y = ballBottom + partRadius + dy * j;
                float dx_center = x - ballCenterX;
                float dy_center = y - ballCenterY;
                // Проверка принадлежности к шару
                if (dx_center * dx_center + dy_center * dy_center <= ballRadius * ballRadius) {
                    // Проверка на выход за границы аквариума
                    if (x >= hh && x <= simWidth - hh && y >= hh && y <= simHeight - hh) {
                        particlePositions.push_back(x);
                        particlePositions.push_back(y);
                    }
                }
            }
        }

        // Обновляем количество частиц
        int maxSimParticles = particlePositions.size() / 2;
        std::cout << "particles = " << maxSimParticles  << std::endl;
        std::cout << "dx = " << dx << " dy = " << dy << std::endl;

        // Создаем солвер с актуальным количеством частиц
        FluidSolver solver(rho, simWidth, simHeight, hh, partRadius, maxSimParticles);
        solver.setUpParticlesAndCells(maxSimParticles, particlePositions);

        solver.runSimulation(dt, gravity, flipCoef, 5000,
                             numPressureIters, numParticleIters, "OutputData/res.txt");

    }

    return 0;
}