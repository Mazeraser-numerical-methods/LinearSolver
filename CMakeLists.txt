cmake_minimum_required(VERSION 4.3)
project(LinearSolver)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

add_executable(LinearSolver main.cpp
        LinearSolver.cpp
        LinearSolver.h
        demos.cpp
        demos.h)
