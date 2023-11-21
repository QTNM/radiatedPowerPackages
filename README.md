# radiatedPowerPackages
Designed to offer easy calculation of EM fields and the associated voltage signals for a moving electron.

## Requirements and Dependencies

The core code requires C++20, built with CMake (3.18+) with the following external dependencies for each package:

|   | BOOST (1.73+) | ROOT (6.14+) | libRootFftwWrapper(*)|
|:-:|:-------------:|:------------:|:--------------------:|
| Basic Functions | | x | x |
| Field Classes | | x | x |
| Signal Processing | | x | x
| Antennas | | x | |
| Electron Dynamics | x | x | |
| Scattering | x | | |

HDF5 libraries are also required for some executables purely for I/O.

(*) libRootFftwWrapper available at: https://github.com/nichol77/libRootFftwWrapper. Requires an existing FFTW install.

There are also inter-dependencies between the individual packages:

|   | Basic Functions | Field Classes | Signal Processing | Antennas | Electron Dynamics | Scattering |
|:-:|:---------------:|:-------------:|:-----------------:|:--------:|:-----------------:|:----------:|
| Basic Functions | | | | | |
| Field Classes | x | | | | |
| Signal Processing | x | x | | x | |
| Antennas | x | | | | |
| Electron Dynamics | x | | | | |
| Scattering | | | | | |

Dependencies are listed row by row, for example the Signal Processing package requires the Basic Functions, Field Classes and Antennas packages.

There are also a set of example programs (radiatedPowerPackages/executables) for which we recommend installing all packages.

A docker image containing the required external libraries is available at: https://hub.docker.com/repository/docker/sebj101/qtnm_deps.

A minimal set of instructions to install the required dependencies using Ubuntu 20.04 are:

### Basic build environment
```
$ sudo apt-get update -y
$ sudo apt-get upgrade -y
$ sudo apt-get install -y git mpich make wget python3-pip build-essential libssl-dev libfftw3-dev libgsl27
```

### Install CMake
```
$ wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz
$ tar -zxvf cmake-3.20.0.tar.gz
$ cd cmake-3.20.0
$ ./bootstrap
$ make
$ make install
```

### Download (pre-compiled) ROOT
```
$ wget https://root.cern/download/root_v6.26.04.Linux-ubuntu22-x86_64-gcc11.2.tar.gz
$ tar -xzvf root_v6.26.04.Linux-ubuntu22-x86_64-gcc11.2.tar.gz
$ echo "source /path/to/root/bin/thisroot.sh" >> ~/.bashrc
$ source ~/.bashrc
```

### BOOST
```
$ wget https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2
$ tar --bzip2 -xf boost_1_73_0.tar.bz2
$ cd "boost_1_73_0/"
$ ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options
$ ./b2 install
```

### libROOTFFTW
```
$ git clone https://github.com/nichol77/libRootFftwWrapper.git
$ cd "libRootFftwWrapper"
$ mkdir build
$ cd build
$ cmake .. 
$ make 
$ make install
```

### HDF5
```
$ sudo apt-get install -y libhdf5-serial-dev
```

## Build instructions
The package is designed to be built with CMake and the build instructions are as follows.

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build . -jN
```
where N is the number of cores used in the build

## Creating electron trajectories
All the code required to create magnetic fields and electron trajectories is contained within the ```ElectronDynamics``` folder. 
It is possible to use the contained classes to generate your own electron trajectories (this is demonstrated in the ```writeTrajectory``` executable).
However, a helper class (```ElectronTrajectoryGen``` contained within ```ElectronDynamics/TrajectoryGen.h```) exists which will create ROOT files containing the trajectories in a format that can be read by other package components.

The required arguments are:
* The output file path
* A pointer to the magnetic field map (more on this later)
* A ROOT TVector3 of the initial electron position
* A ROOT TVector3 of the initial electron velocity
* The simulation step size to use in seconds (for reference an endpoint electron in a 1T field will have a cyclotron period of about $3.7 \times 10^{-11}$ seconds so 20 steps per orbit would be $1.85 \times 10^{-12}$. 
* The total time to simulate in seconds

The optional arguments are:
* The time (in seconds) at which to start the simulation
* The energy loss (which by default is zero). For realistic energy loss this should be set to $2 r_{e} / (3 c)$

Upon successfuly creation of the ```ElectronTrajectoryGen``` object, call the ```GenerateTraj``` member function to produce the file.

### Solver
The solver used to propagate the particles is the Boris solver which is energy conserving (except for any specified energy losses). 
An explanation of the solver (named as the Boris C solver) is given in [Ref. 1][1].

### Magnetic fields
Classes representing several types of magnetic fields and trapping configurations are found in ```ElectronDynamics/QTNMFields.h```.
The most familiar ones are ```BathtubField``` and ```HarmonicField```.

## Using the signal processing

### Antennas
Different types of antenna are implemented as derived classes of the abstract base class ```IAntenna```. 
Typically the one I used for my analyses is the ```HalfWaveDipole```.
Implementations do exist for other antennas but I provide no guarantee of accuracy for them (yet).
The ```HalfWaveDipole``` requires a spatial position, the directions of two (perpendicular) cartesian axes and a central frequency. 
Additionally a time delay (in seconds) may be specified for the antenna.

### FieldPoint
If you want to view the EM fields at a given point in space the best class to use is ```FieldPoint``` located in ```FieldClasses/FieldClasses.h``` which takes an its inputs a ROOT file containing an electron trajectory and a pointer to an antenna.
Calling the function ```GenerateFields``` is required to generate the fields between two specified times.

### InducedVoltage
The ```InducedVoltage``` voltage class provides a lighter weight implementation of just the voltage induced on the specific antenna (or array of antennas). 
The required inputs are an electron trajectory file and either a pointer to an antenna or (in the case of simulating an array of antennas wired together) and vector of antennas.
Additionally, a boolean can be specified over whether or not to include the signal propagation time when generating the signal.

If one just wants to view this signal without any downmixing, downsampling or noise added then simply call the function ```GenerateVoltage```. 
A ```TGraph``` of the time series signal can then be produced using ```GetVoltageGraph```.

### Signal
Full signal processing is done using the ```Signal``` class which takes as an input the path to the electron trajectory file, a pointer to the ```IAntenna```, a ```LocalOscillator``` used to define the down-mixing and the sampling rate (in Hertz).
Additionally a vector of ```GaussianNoise``` terms can be supplied as well as a maximum acquisition time.

[1]: <https://aip.scitation.org/doi/pdf/10.1063/1.5051077>
