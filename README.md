## About NextFOAM solver
NextFOAM solver is a free, open source computational fluid dynamics (CFD) software package released by [NEXTfoam](https://nextfoam.co.kr/foam-Introen.php) based on OpenFOAM released by OpenCFD.

## Copyright
NextFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file COPYING in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

## NextFOAM-24 features than OpenFOAM released by OpenCFD
- Improvement of pressure-velocity coupling
- Improvement of velocity & density interpolation
- Improvement of under-relaxation factor dependency of navier-stokes equation
- Improvement of discretization method of pressure gradient term
- Improvement of time step dependency of transient solver
- Improvement of linearization method of production term of turbulence model
- Development of convergence judgment function of CHT solver
- Improvement of porous media model
- Improvement of MRF (Multi Reference Frame)
- Increase of convergence of radiation model

## Download and Installation instructions on Ubuntu Linux

NextFOAM top directory is set as `/opt/NextFOAM` for all users. Installation directories are followings:

| application | directory |
| --- | --- |
| NextFOAM-24 | /opt/OpenFOAM/NextFOAM-24 |
| ThirdParty-24 | /opt/OpenFOAM/ThirdParty-24 |

Install required packages for building NextFOAM-24 in the Ubuntu Linux
```
$ sudo apt-get -y update
$ sudo apt-get -y install build-essential flex zlib1g-dev libgmp-dev libmpfr-dev
```

Download `openmpi 4.0.5` source and install on `/opt/openmpi-4.0.5` directory
```
$ wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz
$ tar zxf openmpi-4.0.5.tar.gz
$ rm openmpi-4.0.5.tar.gz
$ cd openmpi-4.0.5
$ ./configure --prefix=/opt/openmpi-4.0.5
$ make -j 4 all
$ sudo make install
$ sudo "echo 'export PATH=$PATH:/opt/openmpi-4.0.5/bin' >> /etc/bash.bashrc"
```
Clone `nextfoam-solver` and move to the top directory

```
$ sudo mkdir -p /opt/OpenFOAM
$ git clone https://github.com/nextfoam/nextfoam-solver.git
$ sudo mv nextfoam-solver/NextFOAM-24 /opt/OpenFOAM
$ sudo mv nextfoam-solver/ThirdParty-24 /opt/OpenFOAM
```

Setup the environment variables in the `/opt/OpenFOAM/NextFOAM-24/etc/bashrc`
```
$ sudo vi /opt/OpenFOAM/NextFOAM-24/etc/bashrc

export WM_PROJECT_VERSION=24
projectDir="/opt/OpenFOAM/NextFOAM-$WM_PROJECT_VERSION"
export WM_PROJECT_USER_DIR="/opt/$WM_PROJECT/nextfoam-$WM_PROJECT_VERSION"
```

Download `SCOTCH-6.1.0` and set this version in the `etc/config.sh/scotch`
```
$ cd /opt/OpenFOAM/ThirdParty-24
$ sudo wget https://sources.easybuild.io/s/SCOTCH/scotch_6.1.0.tar.gz
$ tar zxf scotch_6.1.0.tar.gz

$ sudo vi /opt/OpenFOAM/NextFOAM-24/etc/config.sh/scotch

SCOTCH_VERSION=scotch_6.1.0
```

Compile NextFOAM-24. Set your number of cores at `WM_NCOMPPROCS`
```
$ sudo source /opt/OpenFOAM/NextFOAM-24/etc/bashrc
$ cd /opt/OpenFOAM/NextFOAM-24
$ export WM_NCOMPPROCS=4
```
**(Note)** If you install NextFOAM-24 on Ubuntu 22.04, you should install `gcc-9` and `g++-9` and set `gcc-9` as the compiler. 

```
$ sudo apt install gcc-9 g++-9
$ export WM_COMPILE_CONTROL="version=9"
$ sudo ./Allwmake
```

Enable the environment variables at start time
```
$ sudo "echo 'source /opt/OpenFOAM/NextFOAM-24/etc/bashrc' >> /etc/bash.bashrc"
$ sudo "echo 'export WM_COMPILE_CONTROL="version=9"' >> /etc/bash.bashrc"
```

## Contact to NEXTfoam
Contact to NEXTfoam by marketing@nextfoam.co.kr
