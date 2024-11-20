//
// Created by Иван on 11/17/2024.
//

#include "FluidSolver.h"


// SIMULATE PARTICLES

//void (){
//
//}

void FluidSolver::integrateParticles(float dt, float g) {
    for(int i = 0; i < numParticles; ++i){
        particleVel[2 * i + 1] += dt * g;
        particlePos[2 * i] += particleVel[2 * i] * dt;
        particlePos[2 * i + 1] += particleVel[2 * i + 1]*dt;
    }
}

void FluidSolver::handleParticleCollisions(){
    float minX = h + particleRadius;
    float maxX = (numX - 1)*h - particleRadius;
    float minY = h + particleRadius;
    float maxY = (numY - 1)*h - particleRadius;
    //wall collisions:
    for(int i = 0; i < numParticles; ++i){
        float x = particlePos[2*i];
        float y = particlePos[2*i + 1];

        if(x < minX){
            x = minX;
            particleVel[2*i] = 0.0;
        }
        if(x > maxX){
            x = maxX;
            particleVel[2*i] = 0.0;
        }
        if(y < minY){
            y = minY;
            particleVel[2*i + 1] = 0.0;
        }
        if(y > maxY){
            y = maxY;
            particleVel[2*i+1] = 0.0;
        }

        particlePos[2*i] = x;
        particlePos[2*i+1] = y;
    }
}

void FluidSolver::pushParticlesApart(const int numIters){
    
    std::fill(numCellParticles.begin(), numCellParticles.end(), 0);
    // вычисляем количество частиц в каждой ячейке
    for(int i = 0; i < numParticles; ++i){
        float x = particlePos[2*i];
        float y = particlePos[2*i + 1];
        int xi = clamp(static_cast<int>(std::floor(x * h_inv)), 0, pNumX-1);
        int yi = clamp(static_cast<int>(std::floor(y * h_inv)), 0, pNumY-1);
        int cellNr = xi * pNumY + yi;
        ++numCellParticles[cellNr];
    }

    // // заполняем ячейки частицами

    //частичные суммы (utility)
    int first = 0;
    for(int i = 0; i < pNumCells; ++i){
        first += numCellParticles[i];
        firstCellParticle[i] = first;
    }
    firstCellParticle[pNumCells] = first;

    //непосредственное заполнение
    for(int i = 0; i < numParticles; ++i){
        float x = particlePos[2*i];
        float y = particlePos[2*i + 1];
        int xi = clamp(static_cast<int>(std::floor(x * h_inv)), 0, pNumX-1);
        int yi = clamp(static_cast<int>(std::floor(y * h_inv)), 0, pNumY-1);
        int cellNr = xi * pNumY + yi;
        firstCellParticle[cellNr]--;
        cellParticleIds[firstCellParticle[cellNr]] = i;
    }

    // Оттаклкиваем частицы друг от друга

    float minDist = 2.0 * particleRadius;
    float minDist2 = minDist * minDist;

    for (int iter = 0; iter < numIters; ++iter) {
        for (int i = 0; i < numParticles; ++i) {
            float px = particlePos[2 * i];
            float py = particlePos[2 * i + 1];

            int pxi = static_cast<int>(std::floor(px * pInvSpacing));
            int pyi = static_cast<int>(std::floor(py * pInvSpacing));
            int x0 = std::max(pxi - 1, 0);
            int y0 = std::max(pyi - 1, 0);
            int x1 = std::min(pxi + 1, pNumX - 1);
            int y1 = std::min(pyi + 1, pNumY - 1);

            for (int xi = x0; xi <= x1; ++xi) {
                for (int yi = y0; yi <= y1; ++yi) {
                    int cellNr = xi * pNumY + yi;
                    int ffirst = firstCellParticle[cellNr];
                    int last = firstCellParticle[cellNr + 1];
                    for (int j = ffirst; j < last; j++) {
                        int id = cellParticleIds[j];
                        if (id == i){
                            continue;
                        }
                        float qx = particlePos[2 * id];
                        float qy = particlePos[2 * id + 1];

                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;
                        if (d2 > minDist2 || d2 == 0.0){
                            continue;
                        }
                        float d = std::sqrt(d2);
                        float s = 0.5 * (minDist - d) / d;
                        dx *= s;
                        dy *= s;
                        particlePos[2 * i] -= dx;
                        particlePos[2 * i + 1] -= dy;
                        particlePos[2 * id] += dx;
                        particlePos[2 * id + 1] += dy;
                    }
                }
            }
        }
    }
    
}

void FluidSolver::transferVelocitiesToGrid() {
    int n = numY;
    float h1 = h_inv;
    float h2 = 0.5 * h;

    u_prev.swap(u);
    v_prev.swap(v);

    std::fill(du.begin(), du.end(), 0.0);
    std::fill(dv.begin(), dv.end(), 0.0);
    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);

    // ОПРЕДЕЛЯЕМ ТИП ЯЧЕЙКИ
    for(int i = 0; i < numCells; ++i){
        cellType[i] = s[i] == 0 ? SOLID_CELL : AIR_CELL;
    }

    // ЗАДАЁМ ЯЧЕЙКИ С ЖИДКОСТЬЮ
    for(int i = 0; i < numParticles; ++i){
        float x = particlePos[2 * i];
        float y = particlePos[2 * i + 1];
        int xi = clamp( static_cast<int>(std::floor(x * h1)), 0, numX - 1 );
        int yi = clamp( static_cast<int>(std::floor(y * h1)), 0, numY - 1 );
        int cellNr = xi * n + yi;
        if(cellType[cellNr] == AIR_CELL){
            cellType[cellNr] = FLUID_CELL;
        }
    }

    // цикл по двум компонентам: u и v
    for(int comp = 0; comp < 2; ++comp){
        float dx = comp == 0 ? 0.0 : h2;
        float dy = comp == 0? h2 : 0.0;

        std::vector<float> f = comp == 0 ? u : v;
        std::vector<float> f_prev = comp == 0 ? u_prev : v_prev;
        std::vector<float> d = comp == 0 ? du : dv;

        for(int i = 0; i < numParticles; ++i){
            float x = clamp(particlePos[2 * i], h, static_cast<float>(numX - 1) * h);
            float y = clamp(particlePos[2 * i + 1], h, static_cast<float>(numY - 1) * h);

            int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), numX - 2);
            float tx =  ((x - dx) - static_cast<float>(x0) * h) * h1;
            int x1 = std::min(x0 + 1, numX - 2);

            int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), numY - 2);
            float ty =  ((y - dy) - static_cast<float>(y0) * h) * h1;
            int y1 = std::min(y0 + 1, numY - 2);

            float sx = 1.0 - tx;
            float sy = 1.0 - ty;

            float d0 = sx*sy;
            float d1 = tx*sy;
            float d2 = tx*ty;
            float d3 = sx*ty;

            int nr0 = x0 * n + y0;
            int nr1 = x1 * n + y0;
            int nr2 = x1 * n + y1;
            int nr3 = x0 * n + y1;

            float pv = particleVel[2 * i + comp];
            f[nr0] += pv * d0; d[nr0] += d0;
            f[nr1] += pv * d1; d[nr1] += d1;
            f[nr2] += pv * d2; d[nr2] += d2;
            f[nr3] += pv * d3; d[nr3] += d3;
        }

        for(int i = 0; i < numX; ++i){
            if(d[i] > 0.0){
                f[i] /= d[i];
            }
        }

        // восстанавливаем "твёрдые" (solid) ячейки
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                bool isSolid = cellType[i * n + j] == SOLID_CELL;
                if (isSolid || (i > 0 && cellType[(i - 1) * n + j] == SOLID_CELL))
                    u[i * n + j] = u_prev[i * n + j];
                if (isSolid || (j > 0 && cellType[i * n + j - 1] == SOLID_CELL))
                    v[i * n + j] = v_prev[i * n + j];
            }
        }
    }

}

void FluidSolver::transferVelocitiesToParticles(const float flipCoef) {
    int n = numY;
    float h1 = h_inv;
    float h2 = 0.5 * h;

    // цикл по двум компонентам: u и v
    for(int comp = 0; comp < 2; ++comp){
        float dx = comp == 0 ? 0.0 : h2;
        float dy = comp == 0? h2 : 0.0;

        std::vector<float> f = comp == 0 ? u : v;
        std::vector<float> f_prev = comp == 0 ? u_prev : v_prev;
        std::vector<float> d = comp == 0 ? du : dv;

        for(int i = 0; i < numParticles; ++i){
            float x = clamp(particlePos[2 * i], h, static_cast<float>(numX - 1) * h);
            float y = clamp(particlePos[2 * i + 1], h, static_cast<float>(numY - 1) * h);

            int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), numX - 2);
            float tx =  ((x - dx) - static_cast<float>(x0) * h) * h1;
            int x1 = std::min(x0 + 1, numX - 2);

            int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), numY - 2);
            float ty =  ((y - dy) - static_cast<float>(y0) * h) * h1;
            int y1 = std::min(y0 + 1, numY - 2);

            float sx = 1.0 - tx;
            float sy = 1.0 - ty;

            float d0 = sx*sy;
            float d1 = tx*sy;
            float d2 = tx*ty;
            float d3 = sx*ty;

            int nr0 = x0 * n + y0;
            int nr1 = x1 * n + y0;
            int nr2 = x1 * n + y1;
            int nr3 = x0 * n + y1;

            int offset = comp == 0 ? n : 1;
			float valid0 = cellType[nr0] != AIR_CELL || cellType[nr0 - offset] != AIR_CELL ? 1.0 : 0.0;
			float valid1 = cellType[nr1] != AIR_CELL || cellType[nr1 - offset] != AIR_CELL ? 1.0 : 0.0;
			float valid2 = cellType[nr2] != AIR_CELL || cellType[nr2 - offset] != AIR_CELL ? 1.0 : 0.0;
			float valid3 = cellType[nr3] != AIR_CELL || cellType[nr3 - offset] != AIR_CELL ? 1.0 : 0.0;

			float v = particleVel[2 * i + comp];
			float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

			if (d > 0.0) {
                float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1]  \
                            + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
                float corr = (valid0 * d0 * (f[nr0] - f_prev[nr0]) + valid1 * d1 * (f[nr1] - f_prev[nr1])
                            + valid2 * d2 * (f[nr2] - f_prev[nr2]) + valid3 * d3 * (f[nr3] - f_prev[nr3])) / d;
                float flipV = v + corr;

                particleVel[2 * i + comp] = (1.0 - flipCoef) * picV + flipCoef * flipV;
			}
        }
    }
}

void FluidSolver::updateParticleDensity(){
    int n = numY;
    float h2 = 0.5 * h;
    std::fill(particleDens.begin(), particleDens.end(), 0.0);

    for(int i = 0; i < numParticles; ++i){
        float x = particlePos[2 * i];
        float y = particlePos[2 * i + 1];

        x = clamp(x, h, (numX - 1)*h);
        y = clamp(x, h, (numY - 1)*h);

        int x0 = static_cast<int>(std::floor((x-h2)*h_inv));
        float tx = ((x-h2) - x0 * h) * h_inv;
        float x1 = std::min(x0 + 1, numX - 2);

        int y0 = static_cast<int>(std::floor((y-h2)*h_inv));
        float ty = ((y-h2) - y0 * h) * h_inv;
        float y1 = std::min(y0 + 1, numY - 2);

        float sx = 1.0-tx;
        float sy = 1.0-ty;

        if(x0 < numX && y0 < numY) particleDens[x0*n + y0] = sx*sy;
        if(x1 < numX && y0 < numY) particleDens[x1*n + y0] = tx*sy;
        if(x1 < numX && y1 < numY) particleDens[x1*n + y1] = tx*ty;
        if(x0 < numX && y1 < numY) particleDens[x0*n + y1] = sx*ty;
    }

    if(particleRestDensity == 0.0){
        float sum = 0.0;
        int numFLuidCells = 0;

        for(int i = 0; i < numCells; ++i){
            if(cellType[i] == FLUID_CELL){
                sum += particleDens[i];
                numFLuidCells++;
            }
        }

        if(numFLuidCells > 0){
            particleRestDensity = sum / numFLuidCells;
        }
    }
}

void FluidSolver::makeIncompressible(const int numIters, const float dt, const float overRelaxation) {
    std::fill(p.begin(), p.end(), 0.0);
    u_prev = u;
    v_prev = v;

    int n = numY;
    float cp = density * h / dt;

    for(int iter = 0; iter < numIters; ++iter){
        for(int i = 1; i < numX - 1; ++i){
            for(int j = 1; j < numY - 1; ++j){

                // интересуют только "жидкие" ячейки
                if(cellType[i * n + j] != FLUID_CELL){
                    continue;
                }

                // индексы для удобства использования
                int center = i * n + j;   //i,j
                int left = (i - 1) * n + j; // (i-1), j
                int right = (i + 1) * n + j; // (i+1), j
                int bottom = i * n + j - 1; //i, (j-1
                int top = i * n + j + 1; //i, (j+1)

                float sx0 = s[left];
                float sx1 = s[right];
                float sy0 = s[bottom];
                float sy1 = s[top];
                float s_sum = sx0 + sx1 + sy0 + sy1;

                if(s_sum == 0.0){
                    continue;
                }

                float div = (u[right] - u[center]) + (v[top] - v[center]);

                if(particleRestDensity > 0.0){
                    float k = 1.0; //коэф жёсткости
                    float compression = particleDens[i * n + j] - particleRestDensity;
                    if(compression > 0.0){
                        div = div - k * compression;
                    }
                }

                float p_val = - div / s_sum;
                p_val *= overRelaxation; //разобраться с overRelaxation, а затем использовать
                p[center] = cp * p_val;

                u[center] -= sx0 * p_val;
                u[right]  += sx1 * p_val;
                v[center] -= sy0 * p_val;
                v[top]    += sy1 * p_val;
            }
        }
    }
}

void FluidSolver::runFrameSimulation(const float dt, const float g, const float flipCoef,
                                const int numPressureIters,
                                const int numParticleIters) {
    int numSubSteps = 1;
    float sdt = dt / static_cast<float>(numSubSteps);

    for(int step = 0; step < numSubSteps; ++step){
        integrateParticles(sdt, g);
        pushParticlesApart(numParticleIters);
        handleParticleCollisions();
        transferVelocitiesToGrid();
        updateParticleDensity();
        makeIncompressible(numPressureIters, sdt);
        transferVelocitiesToParticles(flipCoef);
    }
}

FluidSolver::FluidSolver(const float dens, const float width, const float height, const float h_step,
                         const float particleRadius, const int maxParticles):
                         density(dens), maxParticles(maxParticles), particleRadius(particleRadius){

    numX = static_cast<int>(std::floor(width / h_step)) + 1;
    numY = static_cast<int>(std::floor(width / h_step)) + 1;
    h = std::max(width / static_cast<float>(numX), height / static_cast<float>(numY));
    h_inv = static_cast<float>(1.0) / h;
    numCells = numX * numY;

    u = std::vector<float>(numCells, 0.0);
    v = std::vector<float>(numCells, 0.0);
    du = std::vector<float>(numCells, 0.0);
    dv = std::vector<float>(numCells, 0.0);
    u_prev = std::vector<float>(numCells, 0.0);
    v_prev = std::vector<float>(numCells, 0.0);
    p = std::vector<float>(numCells, 0.0);
    s = std::vector<float>(numCells, 0.0);
    cellType = std::vector<int>(numCells, 0);

    particlePos = std::vector<float>(2 * maxParticles, 0.0);
    particleVel = std::vector<float>(2 * maxParticles, 0.0);
    particleDens= std::vector<float>(numCells, 0.0);
    particleRestDensity = 0.0;

    pInvSpacing = 1.0 / (2.2 * particleRadius);
    pNumX = static_cast<int>(std::floor(width * pInvSpacing)) + 1;
    pNumY =  static_cast<int>(std::floor(height * pInvSpacing)) + 1;
    pNumCells = pNumX * pNumY; //количество ячеек с частицами

    numCellParticles = std::vector<int>(pNumCells, 0);
    firstCellParticle = std::vector<int>(pNumCells + 1, 0);
    cellParticleIds = std::vector<int>(maxParticles, 0);

    numParticles = 0;

    //cellColor =std::vector<float>(3*numCells, 0.0);
    //particleColor std::vector<float>(3 * maxParticles, 0.0);
    /*
     * for (int i = 0; i < maxParticles; ++i)
				particleColor[3 * i + 2] = 1.0;
     * */

}

void FluidSolver::setUpParticlesAndCells(const int particleAmount, std::vector<float> particlePositions) {
    numParticles = particleAmount;
    particlePos = particlePositions;

    int n = numY;

    for(int i = 0; i < numX; ++i){
        for(int j = 0; j < numY; ++j){
            float s_val = 1.0; //жидкость
            if(i == 0 || i == numX-1 || j == 0){
                s_val = 0.0; //твёрдая граница
            }
            s[i * n + j] = s_val;
        }
    }
}

void FluidSolver::runSimulation(const float dt, const float g, const float flipCoef, const int numFrames,
                                const int numPressureIters, const int numParticleIters,
                                const std::string outputFileName) {
    std::fstream file(outputFileName);
    assert(file.is_open());
    for(int frame = 0; frame < numFrames; ++frame){
        runFrameSimulation(dt, g, flipCoef, numPressureIters, numParticleIters);
        // запись состояния
        writeVectorToFile(file, particlePos);
    }
}


