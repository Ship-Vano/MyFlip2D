//
// Created by Иван on 11/17/2024.
//

#include "Utility.h"

std::vector<float> scale(const std::vector<float>& vec1, const float scalar){
    std::vector<float> res{scalar*vec1[0], scalar*vec1[1]};
    return res;
}

std::vector<float> add(const std::vector<float>& vec1, const std::vector<float>& vec2){
    std::vector<float> res{vec1[0]+vec2[0], vec1[1]+vec2[1]};
    return res;
}



