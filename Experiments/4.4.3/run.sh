#!/bin/sh
echo "\033[1mExperiment 1/2: no buffer zone\033[0m\n"
../../Bin/projects/ExitWaveReconstruction/Reconstruction ParameterFile1.param

echo "\n\n\033[1mExperiment 2/2: buffer zone of size W=256\033[0m\n"
../../Bin/projects/ExitWaveReconstruction/Reconstruction ParameterFile2.param
