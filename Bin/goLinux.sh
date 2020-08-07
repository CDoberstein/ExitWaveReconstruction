#!/bin/sh
export CC="gcc"
export CXX="g++"

cmake -G"Unix Makefiles" -DC++_STANDARD=C++11 -DCMAKE_BUILD_TYPE=Release -DDYNAMIC_LINKING=1 -DPARSE_GCC_ERRORS=1 -DUSE_CIMG=1 -DUSE_PNG=1 -DUSE_OPENMP=1 -DUSE_FFTW=1 -DUSE_BOOST=1 -DUSE_OPENCL=1 -DUSE_TIFF=1 -DUSE_DOXYGEN=0 ../Src
