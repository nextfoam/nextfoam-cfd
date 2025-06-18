## About NextFOAM-cfd solver
NextFOAM-cfd solver is a free, open source computational fluid dynamics (CFD) software package released by [NEXTfoam](https://nextfoam.co.kr/foam-Introen.php) based on OpenFOAM released by OpenCFD.

## Copyright
NextFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file COPYING in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

## NextFOAM-2506 features than OpenFOAM released by OpenCFD
- Unsteady density based compressible solver UTSLAeroFOAM is included
- Solver for BARAM-25.2.3 (https://baramcfd.org) 

## NextFOAM-2501 features than OpenFOAM released by OpenCFD
- Fully compatible with OpenFOAM-v2412
- added fluent mesh(.cas/.msh) conversion utility 'fluentToFoam' to convert to multi-region meshes
- Porosity effect to energy equation for fluid is added
- Fickian/dilute approximation for mass diffusion is available
- Switch for solving species equations is added
- Density based compressible solver **TSLAeroFoam** is added
- Improvement of pressure-velocity coupling
- Improvement of velocity & density interpolation
- Improvement of under-relaxation factor dependency of navier-stokes equation
- Improvement of discretization method of pressure gradient term
- Improvement of time step dependency of transient solver
- Improvement of linearization method of production term of turbulence model
- Development of convergence judgment function of CHT solver
- Improvement of porous media model
- Improvement of MRF (Multi Reference Frame)

## Download and Installation instructions on Ubuntu Linux

NextFOAM-cfd top directory is set as `/opt/OpenFOAM` for all users. Installation directories are followings:

| application | directory |
| --- | --- |
| NextFOAM-2506 | /opt/OpenFOAM/NextFOAM-2506 |
| ThirdParty-2506 | /opt/OpenFOAM/ThirdParty-2506 |

Install required packages for building NextFOAM-2506 in the Ubuntu Linux. Run commands as root:

```
apt-get -y update
apt-get install -y build-essential flex zlib1g-dev libgmp-dev libmpfr-dev texinfo cmake
```

Download `openmpi 4.0.5` source and install on `/opt/openmpi-4.0.5` directory
```
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz
tar zxf openmpi-4.0.5.tar.gz
rm openmpi-4.0.5.tar.gz
cd openmpi-4.0.5
./configure --prefix=/opt/openmpi-4.0.5
make -j 4 all
make install
echo 'export PATH=$PATH:/opt/openmpi-4.0.5/bin' >> /etc/bash.bashrc
```
Clone `nextfoam-cfd` and move to the top directory

```
mkdir -p /opt/OpenFOAM
git clone https://github.com/nextfoam/nextfoam-cfd.git
mv nextfoam-cfd/NextFOAM-2506 /opt/OpenFOAM
mv nextfoam-cfd/ThirdParty-2506 /opt/OpenFOAM
```

Setup the environment variables in the `/opt/OpenFOAM/NextFOAM-2506/etc/bashrc`
```
vi /opt/OpenFOAM/NextFOAM-2506/etc/bashrc

export WM_PROJECT_VERSION=2506
projectDir="/opt/OpenFOAM/NextFOAM-$WM_PROJECT_VERSION"
```

Compile NextFOAM-2506.

**(Note)** If you install NextFOAM-2506 on Ubuntu 22.04, you should install `gcc-9` and `g++-9` and set `gcc-9` as the compiler. 

```
apt install gcc-9 g++-9
export WM_COMPILE_CONTROL="version=9"
echo 'export WM_COMPILE_CONTROL="version=9"' >> /etc/bash.bashrc
```

```
source /opt/OpenFOAM/NextFOAM-2506/etc/bashrc
cd /opt/OpenFOAM/NextFOAM-2506
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/fftw-3.3.10/lib
./Allwmake -j
```

Enable the environment variables at start time
```
echo 'source /opt/OpenFOAM/NextFOAM-2506/etc/bashrc' >> /etc/bash.bashrc
```

## Contact to NEXTfoam
Contact to NEXTfoam by marketing@nextfoam.co.kr
