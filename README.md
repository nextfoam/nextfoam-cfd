## About NextFOAM
NextFOAM is a free, open source computational fluid dynamics (CFD) software package released by [NEXTfoam](https://nextfoam.co.kr/foam-Introen.php) based on OpenFOAM released by OpenCFD.

## Copyright
NextFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file COPYING in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

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

Compile NextFOAM-24. Set your number of cores at `WM_NCOMPPROCS`
```
$ sudo source /opt/OpenFOAM/NextFOAM-24/etc/bashrc
$ cd /opt/OpenFOAM/NextFOAM-24
$ export WM_NCOMPPROCS=4
$ sudo ./Allwmake
```

Enable the environment variables at start time
```
$ sudo "echo 'source /opt/OpenFOAM/NextFOAM-24/etc/bashrc' >> /etc/bash.bashrc"
```

## Contact to NEXTfoam
Contact to NEXTfoam by marketing@nextfoam.co.kr
