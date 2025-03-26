//
// Created by Иван on 11/17/2024.
//

#ifndef MYFLIP2D_UTILITY_H
#define MYFLIP2D_UTILITY_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>

enum CELL_TYPE{
    FLUID_CELL = 0,
    AIR_CELL = 1,
    SOLID_CELL = 2
};

enum VELOCITY_TYPE{
    U_FIELD = 0,
    V_FIELD = 1
};

template<typename DT>
DT clamp(const DT x, const DT min, const DT max) {
    if(x < min){
        return min;
    }
    else if(x > max){
        return max;
    }
    return x;
}
template <typename DT>
void writeVectorToFile(std::ofstream& file, std::vector<DT> v) {
    for (int i = 0; i < v.size(); ++i)
        file << std::setprecision(16) << v[i] << " ";
    file << " " << std::endl;
}

template <typename T>
double dot(const std::vector<T>& grid1, const std::vector<T>& grid2, const int size_x, const int size_y) {
    double dotProd = 0.0;
    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            dotProd += grid1[i * size_x + j] * grid2[i * size_x + j];
        }
    }

    return dotProd;
}

template <typename T>
T max(std::vector<T> grid1, int size_x, int size_y) {
    T maxVal = std::numeric_limits<T>::lowest();
    for (int i = 0; i < size_x; i++) {
        for (int j = 0; j < size_y; j++) {
            if (grid1[i*size_x + j] > maxVal) {
                maxVal = grid1[i*size_x + j];
            }
        }
    }

    return maxVal;
}

#endif //MYFLIP2D_UTILITY_H
