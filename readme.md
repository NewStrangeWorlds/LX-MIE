# LX-MIE #
#### Author: Daniel Kitzmann ####

This readme is currently incomplete and will be expanded in the near future.

# Overview #

LX-MIE is a Mie scattering code that calculates:
* Mie efficiencies, cross sections, and asymmetry parameters
* Scattering phase functions for a given set of scattering angles
* Representation of the phase function as a Legendre series

LX-MIE is optimised to treat very large size parameters. It has been tested to values as high as 10^9.
Details can be found in Kitzmann&Heng (2017), https://arxiv.org/abs/1710.04946 .


# Compilation #
LX-MIE can be compiled with the makefile by typing 'make all' to the
terminal. Compiler options can be found in the 'make.global_options". The makefile is currently configured to 
use g++ from the GNU Compiler Collection. Compiled object files are placed in the 'obj' folder. 
The executable (currently configured to be 'lx_mie') can be found in the root directory.


# Running LX-MIE #
Starting LX-MIE requires a parameter file as console command, e.g. 'lx_mie input/parameter.input'.
The input folder contains three examples: a simple calculation of just the Mie efficiencies/cross sections and two more complicated cases involving the calculation of the phase function and its representation as a Legendre series.
