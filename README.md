# SHG_Focused_Gaussian_Beams package

is a C++-based toolkit for simulating focused Gaussian beams second-harmonic generation (SHG) effciency using the coupled-wave equations (CWEs) that well-describe the three wave mixing (TWM) processes in a second-order nonlinear media. 
CUDA programming allows us to implement parallel computing in order to speed up calculations that typically require a considerable computational demand.

The provided software implements a solver for the CWEs including diffraction terms, linear absorption and nonlinear absorption and thermal evolution. The package is restricted to use it in the continuos-wave (CW) regime.
This code implements a scheme based on Split-Step Fourier Method (SSFM).

For running this package is necessary to have installed g++ compiler in your computer and installed the The Fastest Fourier Transform in the West (FFTW) library (https://www.fftw.org/)


## Setup and execution

To run simulations using the package clone this project typing in a terminal
```
git clone https://github.com/alfredos84/cuSHG.git
```
Once the project was cloned, the user will find a parent folder `SHG` that contains: 
- `SHG/src` which has the main file `SHG.cpp` and the bash file `SHG.sh` used to compile and execute the package by passing several simulations parameters, and
- `SHG/src/headers` which has the headers files needed to execute the package.

Note: the cloned repository is usually created at `home` folder in your Linux System.

### Bash file `SHG.sh`

The bash file is mainly used to massively perform simulations by passing the main file different parameters such as pump power, beam waist, oven temparture, etc.
Before starting the user has to allow the system execute the bash file. To do that type in the terminal
```
chmod 777 SHG.sh # enable permissions for execution
```

Finally, to execute the file execute the following command line
```
./SHG.sh         # execute the files
```

When finished a new folder named in the bash variable `FOLDERSIM` will have been created containing the output files.

In the `SHG.sh` file you will find two different the command lines for the compilation:
```
g++ SHG.cpp -DTHERMAL -I /usr/local/include -L /usr/local/include/fftwf3 -lfftw3 -lfftw3f -lm -o SHG
```
that includes the thermal calculations by using the preprocessor variable `THERMAL`, and
```
g++ SHG.cpp -I /usr/local/include -L /usr/local/include/fftwf3 -lfftw3 -lfftw3f -lm -o SHG
```
that does not include thermal calculations (this mode is faster than the first one).
Currently, there are available two crystals, namely, MgO:PPLN and MgO:sPPLT. 
The type of crystal is set in the header file `Crystal.h`, where the MgO:sPPLT nonlinear crystal is set by default. This file should be modified according to users' needs.


The flags `-lfftw3` and `-lfftw3f` tell the compiler to use the `<fftw3.h> library` that performs the Fourier transform.

Finally, the execution is done using the command line in the `SHG.sh` file is
```
./SHG <ARGUMENTS_TO_PASS>
```
where `$ARGx` and others are variables externaly passed to the main file `SHG.cpp`.
It was written in this way to make easier performing simulations massively.

### Outputs

This package returns a set of `.dat` files with the fundamental (pump) and SH (signal) electric fields, separated into real and imaginary parts.



### Contact me
For any questions or queries, do not hesitate to contact the developer by writing to alfredo.daniel.sanchez@gmail.com or alfredo.sanchez@icfo.eu.
