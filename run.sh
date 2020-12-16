#!/bin/bash

ml GCC/9.3.0
ml CMake
ml intel

cmake -S . -B build -DCMAKE_C_COMPILER=icc
cmake --build build
