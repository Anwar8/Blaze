# XBlaze: Requirements
Requirements for the XBlaze project

## Summary
XBlaze is a finite element method (FEM) program developed specifically for structures in fire. It is explicitly designed from the ground up for scalability on high performance computing (HPC) facilities. The codebase is written in C++ and provides an interface for user customisation.

## Journal
#### 19 Feb 24
 Need to figure out how to use Kokkos, and make sure that graders of project prep are able to see that I am doing that. Copied the Kokkos tutorials to the repo and tried to build the solution for exercise 01 but failed. Tutorials were prepared for a preconfigured system which I have to do for my own system. Building Kokkos with Spack, finding its directory with `spack find --paths`, and putting this directory as the Kokkos directory did not succeed as the `Makefile` for the exercises did not find `Makefile.kokkos`. Downloading Kokkos source and placing it in the directory indicated in the tutorial makefile did not succceed as running `make` for the tutorial results in the error message `g++-13: error: unrecognized command-line option '-mrtm'`. Maybe `gcc@13.2` has an issue with building Kokkos on my Mac? Trying to install `gcc@12.2` using Spack now, maybe that helps. Although I will always be working on my Mac to write code, perhaps I should be connecting to Cirrus and doing all my build work there rather than locally. Even if things work on my Mac, it does not mean it will work on Cirrus or Archer2 which are my target devices. Better email JP and get some advice. Perhaps working on Cirrus directly with remote connection from VSCode could be a good option.

## Requirements
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

## Some reviews
- [Elmer uses Metis](https://youtu.be/84K6OxEKEjQ?t=1358) for graph partitioning to discretise the domain for the processes. Elmer is open source and its code can be found [here](https://github.com/ElmerCSC/elmerfem). 
- FEniCS has an entire [book](https://launchpadlibrarian.net/83776282/fenics-book-2011-10-27-final.pdf). This book also includes information about what algorithms are used, and how they implemented in parallel.
- Code Aster appears to be less well developed and maintained than the other two, but has been [used for structures in fire](https://www.code-aster.de/project-cases/project-cases-detail/analysis-of-steel-reinforced-concrete-exposed-to-fire-mfpa-leipzig-gmbh-copy.html). The website has much more of a commercial bend to it with the product being seminars and tutorials.
- [Calculix](http://www.dhondt.de/) uses [PreProMax](https://prepomax.fs.um.si/), which is an open source pre and post processor. From looking around the internet it may be that this program was not really made for HPC.
## Some tips
- You need to specify the friends of class A in class A, so that the friend class can access class A members. "John is my friend and what's mine is his" vs. ~~"I am friends with Alex so I will access his stuff!"~~
## TO DO
- [x] ~~Create a subtype (possibly templated) inheriting from `std::vector`. This subtype is for storing nodes and elements. It combines the abilities of a vector and a sorted map. The sorted map functionality is pretty bare only concerned with finding a node or an element with a specific id. It does this by simply sorting the vector by the ids and then using the index. before returning the member, it checks that the id indeed matches. If the id is higher or lower, it uses the difference to find it. if it's not found it raises an error.~~ 
    - Apprently, inheriting from `std::vector` is a (bad idea)[https://stackoverflow.com/questions/16812606/adding-custom-methods-to-stdvector-or-typdef]. In stead, it is best to simply write a function to perform the functionality that I need. This may make it easier to write a function that works for both vectors of elements or vectors of nodes.
- [x] Each node keeps a set of unique element IDs of those elements to which it belongs. Each node also has 6 DoFs **by the definition: everything is 3D, and DoFs have a strict ordering**. Each element keeps a vector or an array containing numbers that map its stiffness matrix to these DoFs. For example, the `Basic2DBeam` would have a DoF array of `{0, 2, 5}` which corresponds to its freedoms of $U_x$, $U_y$, and $\theta_z$. This DoF array helps in assembling the global stiffness matrix. ~~For assembly, we loop over all the nodes, and for each node we retrieve the stiffness components for the elements connected to them, and place them in the right place in the global stiffness matrix. This allows us to more efficiently fill up the global stiffness matrix as we are filling the global matrix progressively. This is different from if we were assembling the matrix by looping over the elements. In that case, we would be filling bits and pieces of the global matrix everywhere.~~
- [x] ~~Use test-lead development to~~ build and correct the assembler.
- [x] Update the Transformation matrix to account for DoFs that are required at the node but not included in the element stiffness.
- [x] Update transformation matrix to consider offset by utilising translation matrix multiplied with the rotation matrix to correspond to offset.
- [x] Consider the following: add a list of the indices for where the DoFs of each element goes. That is, each element is given a map for where its local stiffnesses need to go in the begining of the analysis. If, for example, we have a constraint that forces the DoFs of a node to be reflected in the DoFs of a nother node, then this is reflected in this original map that is not needed to be updated too often or at all. This would allow us to tell where each local stiffness goes in the global matrix in a much simpler way.
- [x] The `Assembler` should be separated from the `GlobalMesh`.
- [x] The `BasicSolver` should be separated from the `GlobalMesh`.
- [ ] The assembler should have a way to map nodes to the load and displacement vectors. Perhaps using `std::map`?
- [ ] Nodes should have a container for loads, reaction forces, and for nodal dispalcements. 
- [ ] Add functionality to add nodal load.
- [ ] The command to load a node adds the loaded node to a `std::set` of loaded nodes so each time we want to map the loaded nodes, we only need to loop over the loaded nodes.
- [ ] The assembler should fill the load vector based on nodal loads.
- [ ] After each successful analysis step, nodal displacements are updated by the assembler which maps the nodal displacement vector back to the nodes.
- [ ] The solution procedure should include calculating the element internal forces and strains, and each element should have these saved.
- [ ] Fixing a node adds it to a `std::set` that corresponds to fixed nodes. After each successful analysis step, these nodes calculate their reaction forces. All nodes that do not have a constraint just have reaction forces of zero.
- [ ] Add loggers to retrieve and log certain displacements or reaction forces from the nodes, and internal forces/stresses/strains of elements.
