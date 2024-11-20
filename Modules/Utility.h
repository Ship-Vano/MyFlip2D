//
// Created by Иван on 11/17/2024.
//

#ifndef MYFLIP2D_UTILITY_H
#define MYFLIP2D_UTILITY_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

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


#endif //MYFLIP2D_UTILITY_H
