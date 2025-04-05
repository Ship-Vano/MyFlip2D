cmake_minimum_required(VERSION 3.27)
project(MyFlip2D)

set(CMAKE_CXX_STANDARD 17)

add_executable(MyFlip2D main.cpp
        Modules/FluidSolver.cpp
        Modules/Utility.cpp
)
