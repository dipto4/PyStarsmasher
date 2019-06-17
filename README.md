# Readme for PyStarsmasher
## Introduction
Traditionally, Starsmasher uses a fortran interface and input files like sph.input and sph.init to interact with the user. The purpose of this project is to let python handle all of the input and output functions for Starsmasher. The interface would allow the code to be coupled with other astrophysical codes in libraries like AMUSE. 

To run the collision script, one needs to evolve a star cluster using NBODY6. The version of NBODY6 that tracks collisions is present here: https://github.com/dipto4/NBODY6_collision_detection
## Development History
June 15, 2019 -- Second version of PyStarsmasher uploaded. Fixed MPI bug wherein a new MPI library had to be compiled with `--disable-dlopen`. The new version requires no extra steps. In addition, the user does not need to execute the script with `mpirun`. Instead, inside the script, the user is able to put the number of workers they want to run the worker code. The shared library approach has been replaced by MPI intercommunicators. 

April 15, 2019 -- First version of PyStarsmasher uploaded. Still in early stages of development. The input functions built inside fortran have been replaced with python counterparts that allow the user to interact with the code using python scripts.

## Installation and usage instructions

## New Version

### Installing dependencies
1. Anaconda is recommended as it installs all of the necessary python libraries. Download it from https://www.anaconda.com/distribution/
2.To install anaconda on gonzales, transfer the file to your gonzales home directory and type in `bash <name of anaconda installation file>`. Follow the on screen instructions. When it prompts you if you want to modify the `.bashrc` file, say yes.
3. Log out and log back in for the anaconda environment to take effect. If you would like to see if anaconda has been installed properly, type `python` into the command line and make sure the first line says `Python 2.7.13 |Anaconda 4.4.0 (64-bit)| (default, Dec 20 2016, 23:09:15)` or something similar. Type `exit()` to get out of the python command line. The scripts, currently, are only compatible with python 2.7 so it is recommended that you install the python 2.7 version of anaconda.
4. Load the openMPI module. To view what modules are available type `module list` into the terminal. 
5. Install the latest version of mpi4py from source code. Instructions and files can be found at : ()[https://mpi4py.readthedocs.io/en/stable/install.html]. Follow the instructions for installing via distutils. 

### Installing PyStarsmasher
1. Go into the `src` folder inside the repository directory and type in `make -f makefile.gonzales_ifort`.
2. PyStarsmasher should be ready now. To use it, create a script and run it. An example script has been provided below.

### Installing dependencies
1. Anaconda is recommended as it installs all of the necessary python libraries. Download it from https://www.anaconda.com/distribution/
2. To install anaconda on gonzales, transfer the file to your gonzales home directory and type in `bash <name of anaconda installation file>`. Follow the on screen instructions. When it prompts you if you want to modify the `.bashrc` file, say yes.
3. Log out and log back in for the anaconda environment to take effect. If you would like to see if anaconda has been installed properly, type `python` into the command line and make sure the first line says `Python 2.7.13 |Anaconda 4.4.0 (64-bit)| (default, Dec 20 2016, 23:09:15)` or something similar. Type `exit()` to get out of the python command line. The scripts, currently, are only compatible with python 2.7 so it is recommended that you install the python 2.7 version of anaconda.
4. There is a significant bug inside MPI that does not allow the python interface to work properly with the preinstalled version of MPI inside gonzales. One of the ways to fix it is to download OpenMPI and build it yourself. OpenMPI 1.8.1 can be downloaded from here: https://download.open-mpi.org/release/open-mpi/v1.8/openmpi-1.8.1.tar.gz
5. Transfer the file to gonzales and then extract it using `tar -xvf openmpi-1.8.1.tar.gz`. Go inside the directory and type `./configure --disable-dlopen --prefix=/home/username/openmpi_install`. The configure script will run and generate a makefile. After it has finished, type `make -j 8`. This will compile Open MPI. After it is done, type `make install`.
6. Now go into your `.bashrc` file and add the following lines `export PATH="/home/username/openmpi_install/bin:$PATH"`, `export LD_LIBRARY_PATH="/home/username/openmpi_install/lib:${LD_LIBRARY_PATH}"` at the end of the file. Then type `source ~/.bashrc` for the changes to take effect. Now, the system will use the version of OpenMPI that is compatible with the interace. Make sure you change `username` to whatever your username is. The MPI bug is pretty significant and ways are being explored that would circumvent the installation of a different version of MPI.
7. All your dependencies should be built now and you would be able to use PyStarsmasher.

### Installing PyStarsmasher
1. Go into the `src` folder inside the repository directory and type in `make -f makefile.gonzales_ifort`. This would create the shared library `libstarsmasher.so`. Copy this file into the directory above it. 
2. Now go into the tools directory and type in `make -f makefile.tools`. This would create a shared library called `tools.so`. Copy this file into the directory above it.
3. PyStarsmasher should be ready now.

### Using PyStarsmasher
Below is a sample script that highlights how PyStarsmasher works:
```
from PyStarsmasher import Starsmasher

someSimulation = Starsmasher()


someSimulation.n = 10000        # particles per simulation
someSimulation.nnopt = 23       # Nearest neighbors
someSimulation.dtout = 1        # Output frequency
someSimulation.tf = 200         # Final time (in SPH units)

someSimulation.starmass = 0.5   # units are MSun
someSimulation.starradius = 0.3 # units are RSun

someSimulation.ppn = 16
someSimulation.simulationType = '1es'   # Relaxation run

someSimulation.run(num_of_workers=16)   # setParams() and runSim() have been replaced

```



Here is a script demonstrating how the **old version** works.

```
from PyStarsmasher import Starsmasher

someSimulation = Starsmasher()

someSimulation.n = 10000        # particles per simulation
someSimulation.nnopt = 23       # Nearest neighbors
someSimulation.dtout = 1        # Output frequency
someSimulation.tf = 200         # Final time (in SPH units)

someSimulation.starmass = 0.5   # units are MSun
someSimulation.starradius = 0.3 # units are RSun

someSimulation.ppn = 16
someSimulation.simulationType = '1es'   # Relaxation run

someSimulation.setParams()      # sets these parameters

someSimulation.runsim()         # runs the simulation
 ```
 Open your favorite text editor. For example type `vi test.py` and then save it. Make sure the file is in the same directory as the other python files. 
 
 Before we run the script, we have to make some changes to the submission script. Type `vi gonzales.pbs` and make the following changes: change `python collisionScript.py` to `python test.py` and leave the other things intact. To submit the job, type `qsub gonzales.pbs`. 

### Data
Unlike traditional Starsmasher, the data generated through PyStarsmasher is put into a folder called `data` by default. If you would like to change it, use `someSimulation.dirname = <name of the folder>` and PyStarsmasher would generate the data into that folder.


## Known Bugs
1. I believe there are memory leaks within the program. Need to check...
2. The merger script is buggy. An update is coming soon.
