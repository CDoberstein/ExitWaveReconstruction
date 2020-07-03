#!/bin/sh
echo "\033[1mExperiment 1/3: no scaling of individual exit wave pixels\033[0m\n"
../../Bin/projects/ExitWaveReconstruction/Reconstruction ParameterFile1.param

echo "\n\n\033[1mExperiment 2/3: initial calculation of an exit wave scale mask\033[0m\n"
../../Bin/projects/ExitWaveReconstruction/Reconstruction ParameterFile2.param

echo "\n\n\033[1mExperiment 3/3: initial calculation of an exit wave scale mask and updating the scale mask every 30 iterations\033[0m\n"
../../Bin/projects/ExitWaveReconstruction/Reconstruction ParameterFile3.param
