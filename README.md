# Joint Exit Wave Reconstruction and Multimodal Image Registration

This repository contains an implementation of the variational method for exit wave reconstruction described in 

+ [1] Christian Doberstein. Joint Exit Wave Reconstruction and Multimodal Registration of Transmission Electron Microscopy Image Series (Dissertation). *RWTH Aachen University*, 2020. [[DOI](http://doi.org/10.18154/RWTH-2020-06672)]

+ [2] Christian Doberstein and Benjamin Berkels. A least-squares functional for joint exit wave reconstruction and image registration. *Inverse Problems*, 35(5), 2019. [[DOI](https://doi.org/10.1088/1361-6420/ab0b04) | [arXiv](https://arxiv.org/abs/1812.02786)]

It is based on the [QuocMesh software library](https://archive.ins.uni-bonn.de/numod.ins.uni-bonn.de/software/quocmesh/index.html) and distributed under the terms of the [Common Development and Distribution License](LICENSE.txt).

Please feel free to contact me at <doberstein@aices.rwth-aachen.de> if you have any question or feedback regarding the implementation or the mathematical background.

## Prerequisites

The project can be built and run under Linux or macOS. In order to build the project, cmake and make as well as a C++11 capable compiler such as gcc are required. Additionally, the openmp, fftw, boost and tiff libraries are needed as well as gnuplot and a properly working OpenCL environment if some of the computations are performed on a GPU. Some of the image conversion programs in `Bin/tools/image/converter` require the cimg library and setting `-DUSE_CIMG=1` in `Bin/goLinux.sh`.

## Compiling

The repository contains three directories: `Src`, `Bin` and `Experiments`. The source code is located in the `Src` directory and the files directly related to exit wave reconstruction are contained in the `Src/projects/ExitWaveReconstruction` subdirectory. The `Bin` directory is used to build the project and the `Experiments` directory contains the parameter files for some of the example reconstructions from [1].

The source code can be compiled by calling

    cd Bin
    ./goLinux.sh
    make
    make test

in a terminal. Make sure that none of the parent directories' names contain characters that are not alphanumeric (e.g. brackets or spaces) as the build process may fail otherwise.

Upon successful compilation, the main executable called `Reconstruction` will be located in the `Bin/projects/ExitWaveReconstruction` directory. The `Bin/tools/image/converter` directory contains several executables that can be used to convert images from the native Quocmesh format `.q2bz` to a different format.

## Running

The main executable `Reconstruction` expects exactly one argument, which should be the path to a parameter file that is used for the reconstruction. The parameter file determines all program settings, the minimization algorithm, the origin of the input data, the initial guess and other settings. A template parameter file that contains all possible parameters including a detailed documentation can be found in the `Src/projects/ExitWaveReconstruction/Param` directory.

Once a parameter file has been created, the exit wave reconstruction can be started by calling

    Bin/projects/ExitWaveReconstruction/Reconstruction <path-to-parameter-file>

in a terminal. This will print information about the progress to the terminal and generate a new directory in the same directory as the parameter file, to which the program output is saved. Note that the program may generate a lot of output depending on the settings in the parameter file. In particular, the size of the images and the value of the parameter `saveDataIterations` affect the amount of required disk space the most.

The reconstruction can be stopped at any time by terminating the program (ctrl+c) and resumed later by making the same program call as above.

By default, computations are performed with double precision. In order to perform the computations with single precision, the third line of `Src/projects/ExitWaveReconstruction/Reconstruction.cpp` needs to be changed to `typedef float RealType;`. Afterwards, the source code needs to be compiled again by calling `make` in the `Bin` directory before running another reconstruction.

## Output

The output that is generated during the exit wave reconstruction is distributed across five subdirectories:

+ `InitialGuess`: contains the initial guess for the exit wave, translations and focus values.

+ `InputData`: contains the focus series and a copy of the parameter file. If the focus series was simulated from a known exit wave, this directory also contains a subdirectory `CroppedArg` that contains a cropped version of the exit wave that was used to simulate the focus series (where the buffer zone of the exit wave that was used for the simulation of the input data is removed).

+ `IntermediateResults`: contains directories with the estimates of the exit wave, the translations and the focus parameters after a given number of iterations. The names of the subdirectories indicate the iteration count. The parameter `saveDataIterations` determines how frequently intermediate results are stored.

+ `Plots`: contains plots for the residuals and updates of the objective functional's value and the exit wave, translations and focus values. The residual plots for the exit wave / translations / focus values are only available if the correct exit wave / translations / focus values are known. No plots are created for arguments that are kept constant throughout the whole minimization.

+ `TimeStepData`: contains the data points that were used for the generation of the plots as plain text files.

Some of the directories contain a script `GenerateImages.sh`, which is a shortcut for calling the appropriate GenerateImages executable from the `Bin/projects/ExitWaveReconstruction/ImageGen` directory. In order for this script to work, it is necessary to define an environment variable `QUOC_BIN_DIR` that contains the path to the `Bin` directory, which must not end with a `/`. The GenerateImages executable uses the (hidden) files in the current directory to create several images of the current estimate of the exit wave. It also simulates a focus series from the current estimate of the exit wave and generates images of the differences of the simulated images and the input images.

## Example

The template parameter file `Src/projects/ExitWaveReconstruction/Param/ParameterFile.param` can be used without alteration for a joint exit wave reconstruction and image registration. Using this parameter file, an artificial exit wave corresponding to a standard lattice of copper atoms is simulated, which in turn is used to simulate a focus series consisting of 12 images of the size 256 x 256 pixel. The reconstruction can be started by calling

    mkdir Example
    cp Src/projects/ExitWaveReconstruction/Param/* Example/
    Bin/projects/ExitWaveReconstruction/Reconstruction Example/ParameterFile.param

## Experiments

The `Experiments` directory contains the parameter files (and input data where needed) for the experiments from Chapter 4 and Sections 5.1 and 5.2 of the thesis [1]. The subdirectories are named according to the section of the thesis where the experiment is described. They additionally contain a script `run.sh` that can be used to run all of the experiments of a given section successively. Note that some of the experiments may take several hours or even days to complete. All of the experiments in the thesis were performed with double precision.
