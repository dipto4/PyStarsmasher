# Readme for PyStarsmasher
## Introduction
Traditionally, Starsmasher uses a fortran interface and input files like sph.input and sph.init to interact with the user. The purpose of this project is to let python handle all of the input and output functions for Starsmasher. The interface would allow the code to be coupled with other astrophysical codes in libraries like AMUSE. 

## Development History
April 15, 2019 -- First version of PyStarsmasher uploaded. Still in early stages of development. The input functions built inside fortran have been replaced with python counterparts that allow the user to interact with the code using python scripts.

## Installation and usage instructions
### Installing dependencies
1. Anaconda is recommended as it installs all of the necessary python libraries. Download it from https://www.anaconda.com/distribution/
2. To install anaconda on gonzales, transfer the file to your gonzales home directory and type in `bash <name of anaconda installation file>`. Follow the on screen instructions. When it prompts you if you want to modify the `.bashrc` file, say yes.
3. Log out and log back in for the anaconda environment to take effect. If you would like to see if anaconda has been installed properly, type `python` into the command line and make sure the first line says `Python 2.7.13 |Anaconda 4.4.0 (64-bit)| (default, Dec 20 2016, 23:09:15)` or something similar. Type `exit` to get out of the python command line. The scripts, currently, are only compatible with python 2.7 so it is recommended that you install the python 2.7 version of anaconda.
