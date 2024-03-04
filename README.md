# XBlaze: Requirements
Requirements for the XBlaze project

## Summary
XBlaze is a finite element method (FEM) program developed specifically for structures in fire. It is explicitly designed from the ground up for scalability on high performance computing (HPC) facilities. The codebase is written in C++ and provides an interface for user customisation.

## Project Schedule
#### Week of 19 Feb (10)
- [x] Build `Kokkos` and `Kokkos` Exercises
- [x] Request new Cirrus and Archer2 accounts
- [x] Clone Blaze into Project Prep. Repo
- [x] `Kokkos` first exercise
- [x] Separate Blaze into more files
- [x] Blaze build with `CMake`
- [x] Blaze build tests with `CMake` and `gtest`
#### Week of 26 Feb (09)
- [x] Add geometric nonlinearity to Blaze (1/6)
  - Looked through different textbooks, and then found the lieterature that I want to implement based on (Felippa's NLFEA).
- [x] `Kokkos` Lectures Module 2: Views and Spaces
#### Week of 04 Mar (08)
- [ ] Add geometric nonlinearity to Blaze (2/6)
- [x] `Kokkos` Lectures Module 3: Data Structures + MultiDimensional Loops
#### Week of 11 Mar (07)
- [ ] Add geometric nonlinearity to Blaze (3/6)
- [ ] `Kokkos` Lectures Module 8: Kernels: Sparse and Dense Linear Algebra
#### Week of 18 Mar (06)
- [ ] Add geometric nonlinearity to Blaze (4/6)
#### Week of 25 Mar (05)
- [ ] Add geometric nonlinearity to Blaze (5/6)
#### Week of 01 Apr (04)
- [ ] Add geometric nonlinearity to Blaze (6/6)
#### Week of 08 Apr (03)
- [ ] Report writing (1/3)
#### Week of 15 Apr (02)
- [ ] Report writing (2/3)
#### Week of 22 Apr (01)
- [ ] Report writing (3/3)
#### Report Deadline: Monday 29 Apr


## Journal
#### 4 Mar 24
Finished `Kokkos` lecture 3 which covered parallelising multidimensional loops with the MD policy, subviews, unmanaged views, dual views, and thread safety. All these topics are useful for `Blaze` except for, possibly, dual views which are meant for porting larger pieces of code into `Kokkos`. `MDRange` policy is simply a policy that can be used with many `Kokkos` functions such as those for loop parallelism and reduction operators. The policy has its own access pattern, but can default to that of the memory space. It is also capable of working nicely with tiling and allows for separate access pattern within and between tiles. All policies in `Kokkos` also have scheduling policies that mirror `OpenMP`. Dynamic scheduling allows for *true* work stealing and is supposed to be better than what *OpenMP* could offer at some point. This, however, does nothing on GPUs as they have work stealing by default. There was some talk about using tags for OOP about 35 minutes in, but I did not quite understand that. 

Subviews are a fantastic datatype that can point to a subview part of a view using slicing operators. It, however, has some caveats with what datatype it can be when slicing and so Christian said "use auto. Please!!!". Unmanaged views are another datatype but one that is meant for use with external libraries such as those for IO. It is defined by three NOs: No reference counting, no deallocation when losing scope, and no memory space checks. As such, it is on the programmer to ensure, and strictly so, that they pass absolutely the correct layout and memory space. No label can be added to the unamnaged view, either, but "No label" would have made that four NOs rather than the quintessential 3 (all of which are far more important).

Thread safety is really governed by atomic operations in `Kokkos` because the other solutions such as locks are not scalable and are not compatible with the `Kokkos` programming model. These atomics perform much better on GPUs than CPUs. Performance can be okay on CPUs if there is low conflict, but not otherwise. "scatter contribute" is where, for example, each particle in a molecular dynamic code contributes forces to its neighbours. This could be potentially useful for `Blaze` where elements may need to do something to the nodes. `Kokkos` provides scatterview which is able to choose whether to replicate data or use atomic depending on the architecture. It would use atomics for GPUs, and replicate data for CPUs.

#### 3 Mar 24
Had other work that needed to be done.

#### 2 Mar 24
Finished the second `Kokkos` lecture. This was a very heavy set of content. I found out from this lecture that `CudaUVMSpace` allows you to access the view data from both the CPU and the GPU, however it comes with performance issues as it requires paging every time the view is accessed from a different device resulting in worse performance than explicit copying, even. In stead, it is recommended to use a `mirror` of a view. This mirror is created by `Kokkos::create_mirror_view(view)` and allows us to have a version of the same view that can be accessed on a different device - for example, can be accessed on the host (CPU) if the original view is on the GPU. Note that `Kokkos` would never do a hidden deep copy, and the mirror does not actually copy the data - you have to.

View data in `Kokkos` can be set to be either `LayoutRight` or `LayoutLeft`, which corresponds to row and colum-major data access, respectively. It is also possible to create a custom layout, or use a strided or tiled layout. Layouts are very important because they govern the data access pattern. For a CPU, we want the data access pattern to be "caching", while for a GPU we want it to be "coalescing", and if this is flipped then you can expect a 10-fold performance penalty. The reason for this is that each CPU thread wants to access its own cache line, while GPUs usually run via thread groups called "warp" in `Cuda` or "wavefront" in `HIP`. Given threads $t_1$, $t_2$, and $t_3$, and an array `a` with 9 members, you want the access pattern to be:

**CPU**: caching, give each CPU thread contiguous data
|$t_1$|   |   |$t_2$|   |   |$t_3$|   |   |
|---|---|---|---|---|---|---|---|---|
|`a(0)`|`a(1)`|`a(2)`|`a(3)`|`a(4)`|`a(5)`|`a(6)`|`a(7)`|`a(8)`|

**GPU**: coalescing, give each GPU thread strided data (more at 01:25:45)
|$t_1$|$t_2$|$t_3$|$t_1$|$t_2$|$t_3$|$t_1$|$t_2$|$t_3$|
|---|---|---|---|---|---|---|---|---|
|`a(0)`|`a(1)`|`a(2)`|`a(3)`|`a(4)`|`a(5)`|`a(6)`|`a(7)`|`a(8)`|

More details about the above are provided at 01:22:18. 


If defaults are used: "`Kokkos` index mapping and default layouts provide efficient access if **iteration indices** correspond to the **first index** of an array". One caveat of the views, however, is you are unable to deep-copy between two separate arrays that have different views (`LayoutRight` and `LayoutLeft`), which is why we have to use mirroring. Finally, there are advanced reductions that you can perform in `Kokkos` and you can even define your own. They have all of the reductions available in `MPI` and then some. 

I also now know that if we want to operate across different nodes we will have to use `MPI`, which is covered in Lecture 6 - Internode: MPI and PGAS. I am starting to formulate the idea of the parallel implementation of `Blaze` in my mind. It would be similar to `MPI+OpenMP` but would actually be `MPI+Kokkos`
#### 1 Mar 24
Watching the second `Kokkos` lecture. There are some very important findings from this lecture. `Kokkos` views is the primary data structure that I should be using for `Kokkos`. It behaves like a shared pointer so should be very careful about copying and allocating with it - very easy to just reallocate it by accident. I now know for certain that I should be able to compile `Kokkos` for both `OpenMP` and `Cuda` at the same time. By specifying the Execution space I could make a parallel section run in either execution backend. The execution space can be specified in the execution policy which was an integer only before. The integer is a shortcut for the RangePolicy with its default arguments. Like the execution spaces, there are memory spaces which define where data is kept. Memory spaces need to be defined in a more explicit manner (e.g. `CudaUVMSpace`) because these constructs are quite detailed and some may not have parallels in other architectures (no parallel for that `CudaUVMSpace` in `HIP`). Typedefs would be used in real code to specify which memory space is used. Every function that is called in a parallel region must be marked with a `Kokkos` macro to let the compiler know it needs to compile it for which device. At the 1 hour mark, there is a very good example and set of figures about accessing memory spaces across devices - often you can access the accelerator metadata but not the data itself which is unavailable to the host.

#### 29 Feb 24
Today I read a part of Chapter 13: Corotational formulation overview 1. This chapter has made it clear that the main objective of the corotational formulation is to isolate the rigid body motions from the element deformations. What makes it most useful has been introduced in page 13-8: "...adding and removing rigidy body motions can be visualised as a *front end filter* that lies between the assembler/solver and the element library". That is, the corotational formulation can be decoupled form elements and then be used with the library of linear elements to allow them to become geometrically linear. This is likely how `OpenSees` separates the "Transformation" object (including the *corotational* transformation object) from the rest of element implementation. However, looking at the next chapter, I am afraid things will get much more complicated. It might be wise to skip for a while and move to chapters 16, and 20 through 24.

Other than the reading, I note that I should probably include details on how to build `Blaze` in this `README.md`.
**29 Feb Meeting Notes**
- Agree to narrow down scope of project prep to only cover geometric nonlinearity.
- It is acceptable if I do not achieve the objectives of this project prep after trying. This can be justified in the report if an approporiate reason is given, as well as a plan for the future of `Blaze`.

#### 28 Feb 24
From Felippa's book today I read: 
- Chapter 3 "Residual Force Equations" - the chapter introduces the control variable $\boldsymbol{\kappa}$, state variable $\boldsymbol{u}$, and residual vector $\boldsymbol{r}$. It also, most importantly, establishes the tangent stiffness matrix $\boldsymbol{K}=\frac{\partial \boldsymbol{r}}{\partial \boldsymbol{u}}$. Did not need to go beyond section 3.4.
- Chapter 11 "The Total Lagrangian Plane Beam Element: Formulation" - this chapter introduces the Euler-Bernouli and Timoshenko beam elements, but only formally goes through the Timoshenko beam. It establishes that the Timoshenko beam is an easier beam to work with if corrected with the *residual bending flexibility* (RBF) correction that prevents shear-locking. This beam is only 2D, and I was hoping for a 3D beam-column element. Felippa states that 3D beam-column elements are still research areas. The matrices derived have many variables and look intimidating. There is a flowchart or rather 'roadmap' showing the steps to get the tangent stiffness. Felippa notes that this element has a worse-performing geometric stiffness than the BE beam-column element.
- Chapter 12 "The Total Lagrangian Plane Beam Element: Implementation" - a detailed implementation in `Mathematica` of the previously introduced beam element is presented here. The code is given inside a figure and is very hard to follow without a bigger screen where I would be able to break down the congested and compacted code. The implementation of the beam-column stiffness matrices (both material and geometric) actually only take half a page. It may be possible to directly implement this code myself. This chapter then has many validation exercises. **I am starting to think that just the introduction of geometric nonlinearity (and validation exercises) is enough for the project prep and I should forget about the dynamic solver.** 
 
The following chapters are necessary: 
- Chapter 13: Corotational Formulation Overview I
- Chapter 14: Corotational Formulation Overview II
- Chapter 16: The CR Formulation: BE Plane Beam *(CR = Corotational, BE = Bernouli-Euler)*
- Chapter 20: Overview of Solution Methods
- Chapter 21: Continuation Under Load Control
- Chapter 22: Continuation Under General Control: Predictor
- Chapter 23: Continuation Under General Control: Corrector
- Chapter 24: Incremental-Corrective Methods: Implementation *(there are mistakenly two chapters labelled as 24)*

I have also considered looking beyond Felippa's notes, perhaps at the beam-column elemenet used in `OpenSees` or the one used in `Abaqus`. The reason for that is I want a 3D implementation of a beam-column element not a 2D one. `OpenSees` has a new reference for its beam-column element to a 2021 paper; might be worth checking out. However, my experience is that papers in this area tend to be difficult to follow.

#### 27 Feb 24
Carlos Felippa emailed me his notes! I was too tired to do much work except some very brief reading on the bus on the way to work.

#### 26 Feb 24
Went over NLSA notes by Izzuddin, Crisfield's Nonlinear FEM, IM Smith's Programming the Finite Element Method. After that, I spent the evening collecting additional references to help with building the geometrically nonlinear portion of `Blaze`. Currently trying to get my hands on Carlos Felippa's Nonlinear Finite Elements notes - I even emailed him for them. I was able to get, and had from before, a lot of his previous notes. I collected them in a folder, but they will need some time to go through them. I also went through some of `OpenSees` source code, but will need some additional time with it. Perhaps Frank McKenna's thesis would prove helpful. 

Learned that geometric stiffness is sometimes called "initial stiffness" and "stress stiffness". This is going to be difficult.

#### 25 Feb 24
Did some additional modifications to the `Cmakelists.txt` to better present what `CMake` is doing and make it easier to update it in future. Currently working to build the tests using `CMake` and gtest as well. The issue with the basic build method with `CMake` is that it requires rebuilding the source files for `Blaze` and `BlazeTest`. To overcome this, used the `add_library` command to build a common `BLAZE_LIB` library of sources, and then built `Blaze` and `BlazeTest` by linking to this library. This is the correct way to do this.
Also, completed the first exercise from `Kokkos` and copied the completed solution into the repo. It is worth mentioning that in `Kokkos` reductions, we have to explicitly define a variable which is used by each "thread" to store the partial reduction. While under the hood this is exactly what `OpenMP` does, the implementation is unlike in `OpenMP`. In `Kokkos`, we must explicitly define the intermediate reduction variable that is used by each thread.

#### 24 Feb 24
Successfully built Blaze and its tests using a modified `Makefile` and `config.mk`. The reason I went back to the `Makefile` is because `CMake` was throwing the same error as `make`, so I needed to figure out the issue at the original `Makefile` level first as it is easier and more transparent. I had to separate the `GlobalMesh`, `Assembler`, `BasicSolver`, and `main` from the rest of the build targets and make sure to build all of them with access to `gmsh` header files! Likewise, anything dependent on `maths_defaults.hpp` needed to include the header files from `Eigen3`. Finally, while separating the project objects into files that correspond to only that object, I had missed that the definition of `map_dofs` in `beam_element.hpp` did not match the implementation in `beam_element.cpp`. It was simply an issue of whether or not `map_dofs` belonged to the object `Basic2DBeamElement` or not (it did, but forgot to add `Basic2DBeamElement::` prefix in the `.cpp` file). With this, I figured out how to modify the `Cmakelists.txt` to correctly build and install Blaze! Added the `build.sh` and aliased `bash ./bash.sh` to the word build to make it easier to rebuild and install things from scratch.

#### 23 Feb 24
Friday. Too tired.

#### 22 Feb 24
Tried to configure Blaze to build with `CMake`. Started `Cmakelists.txt`, but currently facing problems in linking to `gmsh`. I think I correctly linked to `Eigen3` by using the command `find_package(Eigen3 REQUIRED NO_MODULE)`. I don't know why I needed the `NO_MODULE`, but that's what is provided on the `Eigen` documentation. There is no such information for `gmsh`. All I found was a message board from 12 years ago where someone suggested to add `gmsh` as a subproject in `CMake`. They referenced a depracated directory in the `gmsh` source code.

#### 21 Feb 24
Did the housekeeping tasks from yesterday:   
  1) Request a new Cirrus **and Archer2** account**s** for PP and dissertation
  2) Prepare a timeline for the work with the report hand-in date in mind.
  3) Move everything to the repo made the Project Prep team
  4) Separate Blaze into a sensible file structure, and updated `config.mk` and `Makefile` albit incorrectly as Blaze fails to compile during the linking stage:
```
Making xblaze:
---------------------------------------------
g++ -std=c++20  -o xblaze basic_utilities.o maths_defaults.o global_coords.o node.o basic_section.o basic_orientation.o basic_shape_function.o beam_element.o assembler.o basic_solver.o main.o  global_mesh.o -L/opt/homebrew/Cellar/gmsh/4.11.1_1/lib -lgmsh 
ld: Undefined symbols:
  Basic2DBeamElement::map_dofs(std::__1::vector<int, std::__1::allocator<int>>, std::__1::set<int, std::__1::less<int>, std::__1::allocator<int>>), referenced from:
      Basic2DBeamElement::create_dof_map() in beam_element.o
clang: error: linker command failed with exit code 1 (use -v to see invocation)
```
I will deal with this when I configure the project to build with CMake.

#### 20 Feb 24
- Built `Kokkos` from scratch with cmake. Must make sure to use the `-S`, `-B`, and `-DCMAKE_INSTALL_PREFIX=` flags to get it to work properly. I can now confidently build `Kokkos` with the `CMake` file they provide with the repository.
- For the exercises, you should download the entire tutorial repo, although as was discovered the `BuildScripts` don't work rendering the need to download everything moot. You will build all of `Kokkos` for each exercise. You will need to have `Kokkos` in `~/Kokkos/kokkos`, and you must modify the `Makefile` in the exercise `Solution` or `Begin` directories to mention your architecture (ARM vs BDW). The Makefile will import another much larger `Makefile.kokkos` from the `Kokkos` directory, and you need to modify some variables there too such as architecture (arm). I had to delete some lines that were causing an issue with my Mac M1 processor. If statements should have taken care of the problem lines but they did not. I set the architecture to `arm-v80` altough the Mac M1 is actually 8.5.
- The `Spack` shell script for the `Kokkos` tutorials simply does not work. It uses the diy command which `Spack` does not recognize.
- When the exercise builds successfully, it will have a `.host` extension which is by design. I don't know why that choice was made. It confused me and made question whether  I built successfully. The final line after a successful build is confusingly "Start Build"
- ⁠To `make` the exercise, use the command `Make -j 24`. I am not sure if the flag is needed, I just realised now that the **second** `Kokkos` lecture shows how to do the exercise including building it with the `Makefile`. Be careful, however, as the `Makefile` was built by design for the architecture they are using for the exercises and so it is an easier process for them.
- TODO: 
  1) Request a new Cirrus account for the PP and dissertation
  2) Prepare a timeline for the work with the report hand-in date in mind.
  3) Move everything to the repo made the Project Prep team
  4) Do the first `Kokkos` exercise
  5) Separate Blaze into more logical files with proper file documentation
  6) Create a `CMakeLists.txt` for Blaze
  7) Include geometric nonlinearity in Blaze
  8) Include dynamic solvers in Blaze

#### 19 Feb 24
 Need to figure out how to use `Kokkos`, and make sure that graders of project prep are able to see that I am doing that. Copied the `Kokkos` tutorials to the repo and tried to build the solution for exercise 01 but failed. Tutorials were prepared for a preconfigured system which I have to do for my own system. Building `Kokkos` with Spack, finding its directory with `spack find --paths`, and putting this directory as the `Kokkos` directory did not succeed as the `Makefile` for the exercises did not find `Makefile.kokkos`. Downloading `Kokkos` source and placing it in the directory indicated in the tutorial makefile did not succceed as running `make` for the tutorial results in the error message `g++-13: error: unrecognized command-line option '-mrtm'`. Maybe `gcc@13.2` has an issue with building `Kokkos` on my Mac? Trying to install `gcc@12.2` using Spack now, maybe that helps. Although I will always be working on my Mac to write code, perhaps I should be connecting to Cirrus and doing all my build work there rather than locally. Even if things work on my Mac, it does not mean it will work on Cirrus or Archer2 which are my target devices. Better email JP and get some advice. Perhaps working on Cirrus directly with remote connection from VSCode could be a good option.






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
