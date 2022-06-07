## Nano Transport
A code to evaluate linear response transport in nanostructures. 

It was originally developed as a post-processing suite of a Lánczos exact diagonalization solver \
for dynamical mean-field theory (DMFT) - GitHub repo: https://github.com/aamaricci/dmft-ed

The software is now an independent package. Its structure is as follows:

* code itself
    * `main.f90` is the actual driver
    * `INPUT_VARS.f90` is a module with the definition of all relevant variables
    * `COMMON.f90` is a module containing utility variables and routines with global scope
    * `GREENS_FUNCTIONS.f90` implements the electronic propagator and self-energy routines
    * `TRANSPORT.f90` implements the transmission function and electric current calculation
    * `BUILD_H.f90` implements the real-space one-body Hamiltonian of the scattering region
    
* external libraries
    * `blas` & `lapack` folders contain a subset of routines from the corresponding libraries,\
which are required by the package to compile

*Extensive documentation to come.*

--

**COPYRIGHT & LICENSING** \
Copyright 2015 - (c), Angelo Valli. \
All rights reserved.

The software is provided with no license, as such it is protected by copyright. \
The software is provided as it is and can be read and copied, in agreement with the Terms of Service of GitHub. \
Use of the code is constrained to author's agreement.
