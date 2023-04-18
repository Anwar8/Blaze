# Blaze: Requirements
Requirements for the Blaze project

## Summary
Blaze is a finite element (FE) program developed specifically for structures in fire. It is explicitly designed from the ground up for scalability on high performance computing (HPC) facilities. The codebase is written in C++ and provides an interface for user customisation.

## Requirements
### Must have
- A geometrically nonlinear shell element - preferrably triangular (3 nodes) for flexibility in modelling
- A nonlinear 2-dimensional metallic material model with temperature dependency
- Efficient I/O performance using MPI-IO, HDF5, or NetCDF
- At least one static solver
- Reliance only on very common numerical libraries such as BLAS, LAPACK, and ScaLAPACK
- Thorough unit testing
- Thorough automated verification problems that can be run for different elements and solvers
- Interface to efficiently apply thermal loading
- Detailed error messages that specify exactly the reason for any divergence or crash
- Wiki page from the start that documents the architecture, reasoning, and theory implemented
- Detailed Read The Docs page for API and function documentation
- Containerisation support so `Blaze` can run on any machine or cloud service
- A way to enforce boundary conditions  
### Should have
- At least one dynamic solver
- Decoupled shape function and elements: an element can be constructed by combining a shape function object with an element object
- No reliance on heavy frameworks such as PetSc that may inflate the executable size and complicate the build process
- Thorough automated validation problems that can be run for different elements and solvers
- Interface for the user to specify material degradation properties
- Powerful numerical debugging tools that, for example, show shape of matrix, its classification (positive definie, negative definite, etc.), and rank
- Full support for addition of imperfections
- Full support for consideration of pre-stress
- Traceability of matrix errors such as pointing out which elements and/or materials are introducing zeros into the diagonal
### Could have
- Breadcrumbs: the ability to restart a simulation from a specific milestone
- Support for the new Eurocode EN 1993-1-14 terminology (LA, LBA, MNA, GNA, GMNA, and GMNIA)
- Implicit-explicit automatically jumping solver
- An array of "try" conditional solver properties to attempt upon convergence issues
- Intricate customisable conditonal-solver-conditioner switcher that would change the solver and conditioner based on state of the matrix
## Ideas
- I can use Eigen3 as my BLAS and LAPACK library. It is covered in the EPCC [Modern C++ for Computational Scientists](https://youtube.com/playlist?list=PLB4tvLCynFjShf7VLy-1gL1g9RKDhAYfY) online course. It also appears to be quite performant as shown in their [benchmark page](http://eigen.tuxfamily.org/index.php?title=Benchmark). Please note that the benchmarks appear to be quite out of date!
- It may be a good idea to create specialised `Blaze` types for integers and floating point numbers. For example, `Blazefloat` would be either a single, double, or quad precision floating point number. Likewise, there must be a `Blazeresidue` which is a floating point number that is of twice the precision as the typical `Blazefloat` and is used for calculating residues and errors.
- Division between the geometric and material stiffenesses. If an element has geometric nonlinearity turned on, then the geometric stiffness matrix is calculated, perhaps **inplace** for the stiffness matrix. However, it may be necessary to separate the geometric and material stiffnesses in order to perform LBA analyses. 
- Does preconditioning with part of the geometric or material stiffness matrix help in convergence? What about preconditiong with the buckling Eigen Vector matrix? Does that help?
- It is perhaps sensible that the first thing the program does is establish the global **distributed** matrices it needs to represent and solve the system. It can do that by checking the size of the local matrix contributed by the elements and also by checking the boundary conditions and how they change the number of rows and columns. 
- It appears that the `essential boundary condition method` is the best for the FEM, but is problematic because it requires modifying the size of the system of equations. However, if this is considered a-priori by the condition enforcer before creating the matrix, then I think we can make this easier. For example, it is perhaps a good idea to make the enforcement of a boundary condition applied directly to the nodes, and then the element checks, at the begining of a calculation interval, if any of its nodes has a Dirichlet boundary condition applied to it. If yes, then it will know to send a different number of contributions to the stiffness matrix because of this.
- I am currently thinking to distribute the calculation across compute nodes and perhaps NUMA regions rather than MPI processes. The data in each NUMA region is shared using MPI shared memory, NOT OpenMP. The reason for this is that MPI localises the memory to the cores, which makes OpenMP memory utilisation rather poor and performance suffers as a result. The arrays should be distributed in block-cyclic fashion so that load is balanced. This is a bit tricky as ideal load balance will need NUMA regions to contain parts of the array that are far away, which will require inter-NUMA and inter-node communication.
- Regular old substructring may help with allowing each system to only need to perform solves on the same node without having to communicate to perform linear algebra operations. Only communication necessary would be that for getting the boundary conditions from the adjacents parts of the strucutre - essentially we are solving a number of full small-scale problems at each iteration. We can update boundary conditions each *n* iterations, perhaps, to reduce communication, and maybe even have more than one substructure corresponding to far-apart parts of the structure on the same NUMA region. Again, it is an issue of load balance vs communication. Especially for a structure in fire, the parts under fire will need much more computation compared to the other parts. Same with the parts of the structure that are no longer elastic, or those that are becoming less stable. Perhaps calculating the condition number of the tangent matrix of the substructures and element contributions may help in finding how to distribute load. Same with elements which are experiencing plasticity.
- It is a good idea to have a header file that contains variables that are especially used for optimising performance on different HPC systems. Things like blocking factor (how big each block in in the cyclic distribution is), number of processes per `compute unit` (where a `compute unit` is a collection of processes that share memory together - can range from 1 process to the number of processes on one full node).