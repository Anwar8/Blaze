# Blaze: Requirements
Requirements for the Blaze project

## Summary
Blaze is a finite element (FE) program developed specifically for structures in fire. It is explicitly designed from the ground up for scalability on high performance computing (HPC) facilities. The codebase is written in C++ and provides an interface for user customisation.

## Requirements
### Must have
- A geometrically nonlinear shell element - preferrably triangular (3 nodes) for modelling flexibility
- A nonlinear 2-dimensional metallic material model with temperature dependency
- Efficient I/O performance using MPI-IO, HDF5, or NetCDF
- At least one static solver
- Reliance only on very common numerical libraries such as BLAS, LAPACK, and ScaLAPACK
- Thorough unit testing
- Thorough automated verification problems that can be run for different elements and solvers
- Interface to efficiently apply thermal loading
- Detailed error messages that specify exactly the reason for any divergence or crash
### Should have
- At least one dynamic solver
- Decoupled shape function and elements: an element can be constructed by combining a shape function object with an element object
- No reliance on heavy frameworks such as PetSc that may inflate the executable size and complicate the build process
- Thorough automated validation problems that can be run for different elements and solvers
- Interface for the user to specify material degradation properties
- Powerful numerical debugging tools that, for example, show shape of matrix, its classification (positive definie, negative definite, etc.), and rank

### Could have
- Breadcrumbs: the ability to restart a simulation from a specific milestone
