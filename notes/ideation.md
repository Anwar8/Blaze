# Ideation
Contains thoughts and ideas from the development of prototype 0. The requirements have now been deprecated and are replaced with those in the feasibility report.

## Requirements (deprecated)
### Must have
- A geometrically nonlinear shell element - preferrably triangular (3 nodes) for flexibility in modelling
- A geometrically nonlinear beam-column element
- A nonlinear 2-dimensional metallic material model with temperature dependency
- A nonlinear 1-dimensional metallic material model with temperature dependency
- Efficient I/O performance using MPI-IO, HDF5, or NetCDF
- At least one static solver
- Reliance only on very common numerical libraries such as BLAS, LAPACK, and ScaLAPACK
- Thorough unit testing
- Thorough automated verification problems that can be run for different elements and solvers
- Interface to efficiently apply thermal loading
- Detailed error messages that specify exactly the reason for any divergence or crash
- Wiki page from the start that documents the architecture, reasoning, and theory implemented
- Detailed Read The Docs page for API and function documentation
- Containerisation support so `XBlaze` can run on any machine or cloud service
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
- Support for element deletion
- Breadcrumbs: the ability to restart a simulation from a specific milestone
- Support for the new Eurocode EN 1993-1-14 terminology (LA, LBA, MNA, GNA, GMNA, and GMNIA)
- Implicit-explicit automatically jumping solver
- An array of "try" conditional solver properties to attempt upon convergence issues
- Intricate customisable conditonal-solver-conditioner switcher that would change the solver and conditioner based on state of the matrix
## Ideas and thoughts
- I can use Eigen3 as my BLAS and LAPACK library. It is covered in the EPCC [Modern C++ for Computational Scientists](https://youtube.com/playlist?list=PLB4tvLCynFjShf7VLy-1gL1g9RKDhAYfY) online course. It also appears to be quite performant as shown in their [benchmark page](http://eigen.tuxfamily.org/index.php?title=Benchmark). Please note that the benchmarks appear to be quite out of date (from 2011)!
- It may be a good idea to create specialised `XBlaze` types for integers and floating point numbers. For example, `Blazefloat` would be either a single, double, or quad precision floating point number. Likewise, there must be a `Blazeresidue` which is a floating point number that is of twice the precision as the typical `Blazefloat` and is used for calculating residues and errors.
- Division between the geometric and material stiffenesses. If an element has geometric nonlinearity turned on, then the geometric stiffness matrix is calculated, perhaps **inplace** for the stiffness matrix. However, it may be necessary to separate the geometric and material stiffnesses in order to perform LBA analyses. 
- Does preconditioning with part of the geometric or material stiffness matrix help in convergence? What about preconditiong with the buckling Eigen Vector matrix? Does that help?
- It is perhaps sensible that the first thing the program does is establish the global **distributed** matrices it needs to represent and solve the system. It can do that by checking the size of the local matrix contributed by the elements and also by checking the boundary conditions and how they change the number of rows and columns. 
- It appears that the `essential boundary condition method` is the best for the FEM, but is problematic because it requires modifying the size of the system of equations. However, if this is considered a-priori by the condition enforcer before creating the matrix, then I think we can make this easier. For example, it is perhaps a good idea to make the enforcement of a boundary condition applied directly to the nodes, and then the element checks, at the begining of a calculation interval, if any of its nodes has a Dirichlet boundary condition applied to it. If yes, then it will know to send a different number of contributions to the stiffness matrix because of this.
- I am currently thinking to distribute the calculation across compute nodes and perhaps NUMA regions rather than MPI processes. The data in each NUMA region is shared using MPI shared memory, NOT OpenMP. The reason for this is that MPI localises the memory to the cores, which makes OpenMP memory utilisation rather poor and performance suffers as a result. The arrays should be distributed in block-cyclic fashion so that load is balanced. This is a bit tricky as ideal load balance will need NUMA regions to contain parts of the array that are far away, which will require inter-NUMA and inter-node communication.
- Regular old substructring may help with allowing each system to only need to perform solves on the same node without having to communicate to perform linear algebra operations. Only communication necessary would be that for getting the boundary conditions from the adjacents parts of the strucutre - essentially we are solving a number of full small-scale problems at each iteration. We can update boundary conditions each *n* iterations, perhaps, to reduce communication, and maybe even have more than one substructure corresponding to far-apart parts of the structure on the same NUMA region. Again, it is an issue of load balance vs communication. Especially for a structure in fire, the parts under fire will need much more computation compared to the other parts. Same with the parts of the structure that are no longer elastic, or those that are becoming less stable. Perhaps calculating the condition number of the tangent matrix of the substructures and element contributions may help in finding how to distribute load. Same with elements which are experiencing plasticity.
- It is a good idea to have a header file that contains variables that are especially used for optimising performance on different HPC systems. Things like blocking factor (how big each block in in the cyclic distribution is), number of processes per `compute unit` (where a `compute unit` is a collection of processes that share memory together - can range from 1 process to the number of processes on one full node).
- It is a good idea to study other parallel packages and FEM tools like PetSc, FireDrake, and [Elmer](http://www.elmerfem.org/blog/documentation/). Elmer appears to be especially powerful, is parallel using domain decomposition, and appears to have **super linear** scaling and **over 150% parallel efficiency** on up to 720 compute nodes as shown in [10.5194/gmdd-6-1689-2013](https://gmd.copernicus.org/articles/6/1299/2013/gmdd-6-1689-2013.pdf). Of course, this is for a particular type of problem (ice-sheet model). Others include: [**FEniCS**](https://fenicsproject.org/) (apprently one of the most widely used), Code_Aster, CalculiX, GetFEM++, Deal.II, Hermes, Hermes++, MOOSE, Kratos Multiphysics, SU2, Nektar++, PetIGA, MFEM, DUNE, SfePy, and FreeFem++. Doing a review on how these are parallelised may be helpful in my own implementation.
- I looked at `gmsh` and `Metis`. The former is written in C++, while the latter is written in C. `gmsh` specialises in mesh operations and has a rich API for retrieving various information about the mesh. It also has numberers such the Cuthill-McKee and nested dissection node numbering (good for HPC apparently). `Metis` is more of a graph/map library with similar capabilities. `gmsh` has it own format, `.msh`, and also has a visualiser and an adaptive mesher in addition for support for element deletion. I can use `gmsh` to create a mesh easily, and then use it to relate my local and global dofs efficiently by performing a loop over the mesh elements as shown in the [sample code file](sample_code.md).
- Element removal can be facilitated by removing their stiffness contribution. Problems arise when enough elements are removed leaving a node, or several nodes, as free bodies. In that case, we would have more freedoms than equations and we would end up with a problem. To resolve this, the bare nodes need to be constrained to zero displacement. As they are not connected to any elements, it does not matter what their value is. At predefined intervals or after a predefined number of elements has been removed, it is not a terrible idea to consider re-assembling the global system of equations. At Tsinghua, they reassemble the system of equations after each element deletion which is madly expensive.