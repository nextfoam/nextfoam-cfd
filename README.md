## About NextFOAM-CFD
NextFOAM-CFD is a free, open source computational fluid dynamics (CFD) software package released by [NEXTfoam](https://nextfoam.co.kr/foam-Introen.php) based on OpenFOAM released by OpenCFD.

## Copyright
NextFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file COPYING in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

## NextFOAM-2401 features than OpenFOAM released by OpenCFD
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

NextFOAM top directory is set as `/opt/OpenFOAM` for all users. Installation directories are followings:

| application | directory |
| --- | --- |
| NextFOAM-2401 | /opt/OpenFOAM/NextFOAM-2401 |
| ThirdParty-2401 | /opt/OpenFOAM/ThirdParty-2401 |

Install required packages for building NextFOAM-2401 in the Ubuntu Linux. Run commands as root:

```
# apt-get -y update
# apt-get -y install build-essential flex zlib1g-dev libgmp-dev libmpfr-dev cmake
```

Download `openmpi 4.0.5` source and install on `/opt/openmpi-4.0.5` directory
```
# wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz
# tar zxf openmpi-4.0.5.tar.gz
# rm openmpi-4.0.5.tar.gz
# cd openmpi-4.0.5
# ./configure --prefix=/opt/openmpi-4.0.5
# make -j 4 all
# make install
# echo 'export PATH=$PATH:/opt/openmpi-4.0.5/bin' >> /etc/bash.bashrc
```
Clone `nextfoam-cfd` and move to the top directory

```
# mkdir -p /opt/OpenFOAM
# git clone https://github.com/nextfoam/nextfoam-cfd.git
# mv nextfoam-cfd/NextFOAM-2401 /opt/OpenFOAM
# mv nextfoam-cfd/ThirdParty-2401 /opt/OpenFOAM
```

Setup the environment variables in the `/opt/OpenFOAM/NextFOAM-2401/etc/bashrc`
```
# vi /opt/OpenFOAM/NextFOAM-2401/etc/bashrc

export WM_PROJECT_VERSION=2401
projectDir="/opt/OpenFOAM/NextFOAM-$WM_PROJECT_VERSION"
export WM_PROJECT_USER_DIR="/opt/$WM_PROJECT/nextfoam-$WM_PROJECT_VERSION"
```

**(Note)** When running `decomposePar` command, the no scotch library error is shown, download `SCOTCH-6.1.0` and set this version in the `etc/config.sh/scotch`
```
# cd /opt/OpenFOAM/ThirdParty-2401
# wget https://sources.easybuild.io/s/SCOTCH/scotch_6.1.0.tar.gz
# tar zxf scotch_6.1.0.tar.gz

# vi /opt/OpenFOAM/NextFOAM-2401/etc/config.sh/scotch

SCOTCH_VERSION=scotch_6.1.0
```

Compile NextFOAM-2401. Set your number of cores at `WM_NCOMPPROCS`

**(Note)** If you install NextFOAM-2401 on Ubuntu 22.04, you should install `gcc-9` and `g++-9` and set `gcc-9` as the compiler. 

```
# apt install gcc-9 g++-9
# export WM_COMPILE_CONTROL="version=9"
# echo 'export WM_COMPILE_CONTROL="version=9"' >> /etc/bash.bashrc
```

```
# source /opt/OpenFOAM/NextFOAM-2401/etc/bashrc
# cd /opt/OpenFOAM/NextFOAM-2401
# export WM_NCOMPPROCS=4
# ./Allwmake
```

Enable the environment variables at start time
```
# echo 'source /opt/OpenFOAM/NextFOAM-2401/etc/bashrc' >> /etc/bash.bashrc
```

## Contact to NEXTfoam
http://www.nextfoam.co.kr

http://blog.nextfoam.co.kr

marketing@nextfoam.co.kr
