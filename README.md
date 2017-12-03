# xscf
Scf program bassed on libcint.
This project is a small leightwieght electronic structure code that currently performs restricted open shell hartree fock 
calculations for molecules.
Future plans are:
Address open inner-shells, relevant for X-ray absorption.

to build the code you need 
-cmake
-lapack/blas library

Building instructions

mkdir build
cd build 
cmake ..
make


to run the code you need an input file that contains

specification of the molecular geometry:

geom=(
atom1 x1 y1 z1
atom2 x2 y2 z2
...
atomn xn yn zn
)

where the coordinates (xi,yi,zi) are given in atomic units (bohr radius)

specification of the basis set by

basis_set=basis_set_file

where basis_set_file contains the specification of the basis set (Gaussian format).
There are some basis sets specifications in subdirectory basis_set/

Specifications of the number of basis functions:

nmo=x

where x is the number of basis functions.
If x is not correct the program will abort.

Specifications of the occupation of the molecular orbitals.
This is a comma seperated list containing either 2 (occupied) or 0 (unoccupied)
e.g.

occ=2,2,2,0,0

Example inputfiles are found in subdirectory test_H2O, test_He, test_Ne

