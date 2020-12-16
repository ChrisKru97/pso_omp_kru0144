#!/bin/bash

ml GCC/9.3.0
ml CMake

cmake -S . -B build
cmake --build build