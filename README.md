# Description
Scf program based on the libcint library.
This project is a small leightweight electronic structure code that currently is only able to perform restricted open shell hartree fock calculations for molecules.
Future plans are:
Address open inner-shells that are relevant for X-ray absorption.

# Building instructions
to build the code you need 
*cmake
*a lapack and blas library

to build the code, follow these steps:
```
mkdir build
cd build 
cmake ..
make
```
# Usage
to run the code you need an input file that contains

 * Specification of the molecular geometry:
  ```geom=(
  atom1 x1 y1 z1
  atom2 x2 y2 z2
  ...
  atomn xn yn zn
  )
  ```
  where the coordinates (xi,yi,zi) are given in atomic units (bohr radius)

* Specification of the basis set by
  ```
  basis_set=basis_set_file
  ```
  where basis_set_file contains the specification of the basis set (Gaussian format).
  There are some basis sets specifications in subdirectory basis_set/

 * Specifications of the number of basis functions:
  ```
  nmo=x
  ```
  where x is the number of basis functions.
  If x is not correct the program will abort.

 * Specifications of the occupation of the molecular orbitals.
  This is a comma seperated list containing either 2 (occupied) or 0 (unoccupied)
  e.g.
  ```
  occ=2,2,2,0,0
  ```

Example inputfiles are found in subdirectory test_H2O, test_He, test_Ne

