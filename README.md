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
- [x] Add geometric nonlinearity to Blaze (2/6)
- [x] `Kokkos` Lectures Module 3: Data Structures + MultiDimensional Loops
#### Week of 11 Mar (07)
- [x] Add geometric nonlinearity to Blaze (3/6)
  - [x] Load application and restraint for nodes.
  - [x] Additional reading particularly on implementation of $\boldsymbol{T}$ and $\boldsymbol{K}_G$
- [ ] ~~`Kokkos` Lectures Module 6: Internode: MPI and PGAS (approximately 40 minutes out of the total is about MPI)~~
#### Week of 18 Mar (06)
- [x] Add geometric nonlinearity to Blaze (4/6)
  - [x] Read on Stress and element resistance forces
  - [x] Read on Geometric and tangent stiffness?
  - [x] Read on Nonlinear corotational transformation?
- [x] `Kokkos` Lectures Module 6: Internode: MPI and PGAS (approximately 40 minutes out of the total is about MPI)
- [x] `Kokkos` Lectures Module 8: Kernels: Sparse and Dense Linear Algebra
#### Week of 25 Mar (05)
- [x] Add geometric nonlinearity to Blaze (5/6)
  - [x] Implement the calculation of element strains, stresses and nodal forces
  - [x] Implement the geometric stiffness for the beam element
  - [x] Map-out the geometrically nonlinear transformation object
#### Week of 01 Apr (04)
- [ ] Add geometric nonlinearity to Blaze (6/6)
  - [x] Final assembly and load-incrementation
#### Week of 08 Apr (03)
- [ ] Report writing (1/3)
#### Week of 15 Apr (02)
- [ ] Report writing (2/3)
#### Week of 22 Apr (01)
- [ ] Report writing (3/3)
#### Report Deadline: Monday 29 Apr


## Journal
#### 6 April
1. Upon entering the second iteration, the nonlinear transform appears to reset and calculate non-existing displacement as if the corotated and base configurations are identical although in the first iteration that was not the case.
2. It turns out that this is happening because the nodal reaction forces are not being considered which is resulting in out-of-balance forces equal to the load applied. What this is doing is giving a dU equal and opposite the resolved displacement fields thus returning the structure to its base configuration.
3. There was an error where out-of-balance was calculated as $\boldsymbol{G} = \boldsymbol{K}\boldsymbol{U} - \boldsymbol{R}$, but has been updated to $\boldsymbol{G} = \boldsymbol{R} - \boldsymbol{P}$.
4. It appears that the resistance forces are not being extracted from the elements correctly. Should try $\boldsymbol{G} = \boldsymbol{K}\boldsymbol{U} - \boldsymbol{P}$ for out-of-balance.
5. That did not work. I am now able to continue, but solution is  oscillating between negative and positive values. Do not know what is happening.
6. Tried the Basic Beam as well, and it fails to proceed after 0.88 knowing that buckling should take place at 0.86.
7. I give up for now. Note that I tried 4 KU - P, and it was not of great help. It does, however, produce different end contraction (0.35 vs 0.3) which I found a bit unsettling.
   
#### 3 April 
How to resolve the nonlinear effects? Well, first I will need the trigonometric functions needed. That means I will need to establish a relationship between the base configuration and the corotated and current configurations. How do I define a configuration?
#### 2 April 24
Added load factor, iterations, load incrementing, response tracking, and silenced all annoying logging with the `VERBOSE` flag. However, the current `Blaze` cannot capture geometric nonlinearity at all. The beam is supposed to have a buckling load of 2.58096e7 N, but I am loading with up to -30.0e7 without any sign of buckling. Tracking the end-node y-displacement results in a purely linear response. Note that I have added some small (-1000 N) node on the y-axis of the last node without luck. No buckling. It is clear to me that the current system is completely unable to capture geometric nonlinearity. I will need to implement `NonlinearTransform` first and foremost, and then follow that up quickly with `Izzuddin2DNonlinearBeam`. I am running out of time, I only have this week remaining to work on the code before I have to start writing the thesis.

#### 31 Mar and 1 April 24
I spent a long time attempting to map out and build the new nonlinear transformation object under `NonlinearTransform`. This object will have a rather complex interaction with the `BasicOrientation` object, which may become obselete. In fact, upon attempting to build the nonlinear solver I am finding that my implementation of the geometric stiffness is very lacking. The element-level contributions to the geometric stiffness can actually be ignored, and it is the part that I have been ignoring - what is called "external geometric stiffness" in Du and Hajjar (2020), "principal contribution" in Izzuddin's NLSA notes, and "varying moments through the transverse shear force in $\mathcal{C}$" (16.27) in Felippa's notes is what actually matters. Looking at Izzuddin's in-class notes, he makes that clear as well. I have made extensive hand-written notes on the implementation yesterday in OneNote and I think I am almost ready to implement geometric nonlinearity, although I remain rather confused about how I should structure my program classes. I have decided, for now, to implement Izzuddin's nonlinear element under its own nonlinear class rather than implement everything within the `Basic2DBeamElement` class. Thus, I have introduced `Izzuddin2DNonlinearBeam`, which inherits from `Basic2DBeamElement` but seeks to re-implement the way state is assessed using the nonlinear equations in Izzuddin's notes. I am also likely to implement Felippa's or Izzuddin's corotational transform in `NonlinearTransform` to isolate the deformation-inducing displacements from the rigid body motions. Felippa's implementation is a little more clear because of the better figures, and because I spent a rather long time on 31 March taking hand-written notes on it and drawing it for myself. I need to do this, even though I had decided not to implement corotational transformation at first, because I need the trigonometric relationships from the transformation to calculate the element response such as the principal contribution to global geometric stiffness. I have taken some note on both Izzuddin and Felippa's implementation below:

*Why am I unhappy with Felippa’s element?*
1. Linear generalised stresses where axial force does not play a role in the moments.
2. Because it uses a linear interpolation for the nodal forces and for the stiffness matrix.
3. ⁠it seems to have me do some extra work where phi is simply a function of the original base configuration and current configuration that I can solve for easily already.

*Why am I happy with Felippa's element?*
1. The interpolation separates the rigid body modes and I can probably use this approach where I separate the transform from the element state calculation as all element freedoms are not related to the rigid body modes.

*What do I not like about Izzuddin's element?*
1. Appears to be no shear but actually I can use Felippa’s expression for shear to get the generalised shear from the moments.
2. ⁠No evaluation for the tangent stiffness matrix so I have to multiply it myself.
  
*What do I like about Izzuddin's element?*
1. The transformation matrix appears to be quite easy to interpolate actually.
2. And the element includes nonlinearity at every stage including axial force term in the moment calculation.
   
Finally, I have prepared a draft interface for performing the iterative procedure to calcualte the nonlinear response. I have placed this in `main.cpp` as commented-out code which is meant to guide my future implementation.
**References**
Du and Hajjar (2020) Three-dimensional nonlinear displacement-based beam element for members with angle and tee sections. Engineering Structures.

Additional work on night of 1 April:
Added iterations for calculating resistance forces $\boldsymbol{R}$, out of balance forces $\boldsymbol{G}$, and $\Delta \boldsymbol{U}$. Then used them to update the system state, and recalculate the out-of-balance. Out-of-balance appears to decrease with iterations. TODO: **Still need to apply convergence criteria, and load factor $\lambda$.** Remember that my beam is poking to the right, which is the opposite of my imagination for some reason. The load, to apply compression, should be in the negative x-direction. This has been updated. The output values all appear to be more or less correct, but I will need to check again.

#### 27 Mar 24
It might be that I am not yet done with the reading. I have found the book by Belytschko to include many sections that I must read, particularly on transformations! This book also covers some meshless methods which I will need to revisit. I added geometric stiffness and tangent stiffness calculations. Note that the geometric stiffness is incomplete as it only contains contributions of the axial load.

#### 26 Mar 24
Implemented element matrices for nodal forces, nodal displacements, constitutive matrix, local strains, and local stresses. These matrices are updated via the command `update_state`. Currently, the state update is done outside its right place in the solver and is actually done in `main.cpp`. I am facing a problem, however, as the $\boldsymbol{B}$ matrix needs an `x` value to be evaluated for. I implemented locations for Gauss Points, but should I have two copies of the element stress one for each Gauss point? Most likely not,  since I do not need to carry out any integration operations. In cases where I do need to perform such operations, however, I *will* have a value for each Gauss point as is the case for Abaqus.

Moreover, implementing most of the calculations as matrix operations is okay, but as noted by Felippa's notes, may be wasteful. I could just implement them as direct calculations - after all, I only have $P$ and $M$ corresponding to $\boldsymbol{\varepsilon} = \begin{bmatrix}\varepsilon_c & \kappa\end{bmatrix}^T$. $P$ is also constant throughout the element, which will make the estimation of the tangent stiffness much easier. So, in stead of doing $\boldsymbol{\sigma} = \boldsymbol{D}\boldsymbol{\varepsilon}$, I could simply do something like: $P = EA\varepsilon_c$, where $\varepsilon_c = \frac{-1}{L} u_1 + \frac{1}{L} u_2$ where $u_1$ and $u_2$ are the axial elemental freedoms from $\boldsymbol{d}$, and the coefficients are the corresponding element from the derivative of the shape function matrix $\boldsymbol{B}$.

Calculated element mid-point moment correctly: -285000 N.m for 0.15 away from support - this corresponds to $(3 - 0.15) * -1e5 = -2.85e5$ for a total length of 3 m, load of -100 kN, and 10 elements for the beam. Element nodal forces also equilibrate one-another! Nodal force at the last element is equal to 100 kN, which is exactly right for equilibrium with the 

TODO: calculation of fixed node reaction forces.

#### 25 Mar 24
Began to implement the necessary expressions and element members to calculate the stress state. 

#### 24 Mar 24 
Completed `Kokkos` lectures needed for Project Prep. Entries for 20 through 24 will be added later. Busy week, and limited work, but finished all needed reading to the extent that I have nothing left to read, and I have finished all `Kokkos` lectures. Most likely no explicit corotational transformation will be implemented as we do not need to nonlinearise a large library. Different types of element formulation (total Lagrangian, updated Lagrangian) will be used to corotational formulation not needed as its own thing. A bit too complicated to implement right now to be honest.

#### 21 Mar 24
1. The orientation object must only calculate the orientation angles - it should not create nor control transformation matrices.
2. Constraining DoFs - for example saying lateral displacement of node 4 (so, if all are ordered from 0, it would be global DoF 3 previous nodes * 6 DoFs per nodes = 18 considering we start from global DoF 0) is dependent on vertical displacement of nodes 2 and 3 (global DoFs 2 and 8), then there should be a constraint object that would apply after $\boldsymbol{T}^T \boldsymbol{k}' \boldsymbol{T}$ is done, and its job would be to access the global stiffness triplets and transform their numbers to correspond to the interpolation rules between 2 and 8; this should also be done for nodal forces and displacements, etc. This is why I am considering that there should be a separate object whose job is to handle transformations and perhaps even constraints. we  need to remove global DoF 18 from the global matrices and order things accordingly, but also keep the node-object itself. 
3. I completely forgot about calculating reaction forces - how is that done for constrained nodes? do we go over them one by one and then set their reactions as the sum of the element nodal resistances at that node? 
4. All nodes and elements have these local arrays. These arrays will need to be copied around a LOT between compute nodes. It might be worth it making it such that these arrays are actually Kokkos subviews of a global Kokkos view that contains contiguous data for all of these local quantities. This might potentially allow for efficient data access patterns (and copying, and `MPI` communications)

#### 19 Mar 24
Here are my notes on `Kokkos` lecture 6. The idea behind mixing `MPI` with `Kokkos` is rather simple: send and receive views as normal, and be very certain to overlap communication and computation. `MPI` can be made aware of GPUs, and can communicate with them. Therefore, it is possible to exchange GPU data directly with `MPI`, but need to be careful about how that is done. Doing this across different nodes with the `UVMSpace` will destroy performance. The exchange of data via `MPI` can be done via subviews or views if you are sending the entire data, but one must be careful that only contiguous layouts could be sent so do not send `LayoutStride` data without first creating an appropriate `MPI` datatype.  When overlapping communication and computation, be sure that the compute kernel does not touch the data of the communication kernel. They have an exercise for "packing and unpacking" which is a concept, to be honest, that I did not understand. There is also an exercise about overlapping computation and communication where a 3D heat conduction problem is solved. It involves some information about exchanging boundaries and I believe it is going to be quite a useful exercise for me. A final very important consideration is resource affinity when using `Kokkos`+`MPI`. They recommend not over-subscribing the CPU cores by requesting more than physically available knowing that when mixing with `MPI`, each `MPI` rank will have some `OpenMP` or `Devices` associated with it using the bind command. There are ways to set these numbers but were not updated in the slides uploaded. There is also a way for telling `Kokkos` how many devices there are by using a flag or a system variable. 

#### 18 Mar 24
Reading Rankin (1986) paper. The paper starts off rather clear, but then the notation quickly becomes cumbersome. It is important to also note some limitations: this paper only covers rotations not deflections, and the focus is primarily on shell elements. I still need to continue with this paper, but it is not looking good for me. Right now, I have reached the implementation part of the paper and I am hoping this will be easier to go through.

I have also started going through `Kokkos` lecture 6 about MPI, but will need to continue that tomorrow.


Rankin, C. C., & Brogan, F. A. (1986). An element independent corotational procedure for the treatment of large rotations

#### 17 Mar 24
I chose to take a break today. My efficiency has not been amazing, and I needed some rest.

#### 16 Mar 24 
The transformation matrix, with offset $b$, is currently: 
$$\boldsymbol{T} \boldsymbol{U} = \boldsymbol{d}$$

$$ (6 \times 12) (12 \times 1) = (6 \times 1)$$

$$\begin{bmatrix} \cos\theta & 0 & \sin\theta & 0 & 0 & b & 0 & 0 & 0 & 0 & 0 & 0\\ -\sin\theta & 0 & \cos\theta & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & \cos\theta & 0 & \sin\theta & 0 & 0 & b\\ 0 & 0 & 0 & 0 & 0 & 0 & -\sin\theta & 0 & \cos\theta & 0 & 0 & 1\end{bmatrix} \begin{bmatrix} U_{1,x} \\ U_{1,xx} \\ U_{1,y} \\ U_{1,yy} \\ U_{1,z} \\ U_{1,zz} \\ U_{2,x} \\ U_{2,xx} \\ U_{2,y} \\ U_{2,yy} \\ U_{2,z} \\ U_{2,zz}\end{bmatrix} = \begin{bmatrix} \cos\theta U_{1,x} + \sin\theta U_{1,y} + bU_{1,zz}\\ -\sin\theta U_{1,x} + \cos\theta U_{1,y} \\ U_{1,zz} \\ \cos\theta U_{2,x} + \sin\theta U_{2,y} + bU_{2,zz}\\ -\sin\theta U_{2,x} + \cos\theta U_{2,y} \\ U_{2,zz} \end{bmatrix}= \begin{bmatrix} u_1 \\ v_1 \\ \theta_1 \\ u_2 \\ v_2 \\ \theta_2\end{bmatrix}$$

Which works well if we assume $U_{i,zz}$ is small such that $bU_{i,zz} = \delta x_{i, zz}$. I think the following should work, too, but I need to check that $\boldsymbol{T}^T = \boldsymbol{T}^{-1}$.
$$\boldsymbol{T}^T \boldsymbol{d} = \boldsymbol{U}$$

$$ (12 \times 6) (6 \times 1) = (12 \times 1)$$
Perhaps I will need a small Python toy program to convince myself of this, especially given that I am including the offset in this calculation. I should read Chapter 8 Coordinate Transformation and Selected Analysis Options from Concepts and Applications of Finite Element Analysis (2001).

#### 15 Mar 24
Today I attempted to add the functionality to calculate the element resistance vector, which I don't actually have a variable name for yet, and then use it to calculate the resistance vector $\boldsymbol{R}$. I quickly ran into some issues as the object design, particularly for `BasicShapeFunction` is such that all meaningful element matrices, such as $\boldsymbol{B}$, $\boldsymbol{N}$, and $\boldsymbol{k}$ live within the `BasicShapeFunction` object. The element should *use* the shape function and transformation objects, but right now each of them is its own entity. Additionally, I seem to have confused what properties are constant throughout the steps of teh analysis such as the state of shape function $\boldsymbol{N}$ and $\boldsymbol{B}$ at the Gauss points, which are only affected by the displacements of the element nodes when they are used such as in the equation to calculate the element strains via $\boldsymbol{\varepsilon} = \boldsymbol{B}\boldsymbol{d}$. There are currently no Gauss Points in the element and these should be added, and then I should figure out how I will be using them for my beam-column elements.

Then comes the issue of constraints and how they are applied. I mapped the state vector $\boldsymbol{U}$ back to the nodes, but now these need to be mapped to the elements. That requires using the transformation matrix $\boldsymbol{T}$. I am also using this transformation matrix to enforce boundary conditions, and so I really need to remind myself of the implementation details that I used for the `Basic2DBeamElement` and ALL associated objects: `BasicShapeFunction`, and the inaptly-named `Orientation` class. Perhaps the transformation matrix $\boldsymbol{T}$ should be decoupled from the constraint handling and substituted by another matrix which job is to enforce constraints in the same way? (i.e. by using only the number of rows and columns necessary to transfer only the active coordinates). However, that would add another set of matrix multiplication operations each solution step and iteration, and for each element. Perhaps I should not do this. 
It is additionally confusing that in nonlinear analysis, I should be using a nonlinear transformation matrix: $\boldsymbol{T} = \partial \boldsymbol{d}/ \partial \boldsymbol{U}$. Based on this, it is essential for me to understand the corotational transformation, and although introduced in Felippa's notes, I should probably look at it in the original source by Rankin and Brogan (1986) which seems to do an even better job in explaining it and its implementation. At least that's how I feel about the paper from a very quick and cursory scan.

I completed the remainder of chapter 20, and finished chapter 21 as well of Felippa's NLFEA. Useful chapters, I have to say, especially in that they introduced (such as the explicit midpoint rule) new integration schemes and error assessing equations. I will not write my notes on them here just yet. I have a lot of other things to think about right now.
Also, note that I will not commit changes to the code today as I want the files that I am modifying to be clearly highlighted in VSCode.

#### 14 Mar 24
Added the functionality to map global nodal displacements back from the state vector $\boldsymbol{U}$ to the nodes. Modifications made to `Assembler` and `Node` to make it happen. Functionality appears to work correctly, although the `print_info` function from the `GlobalMesh` is a bit...much (prints too much stuff). 

TODO: Improve logging/`print_info` functions.

#### 13 Mar 24
Tested a 3 m cantilever beam with E = 206 GPa, I = 0.000457 $m^4$, and 100 kN tip point load. The actual end displacement for this problem is -0.00956003 m. The results from `Blaze` with 10 elements and 11 nodes is exactly -0.00956003 m. Few things to note: 
1. I don't currently have a way to output results other than printing to output stream
2. There is no mapping between state vector (displacements) and nodes; have to map it manually in my head
3. There is no calculation of stresses and internal resistances based on displacements. This step is missing.

#### 12 Mar 24
Corrected the bugs and added functionality to assembler and global mesh to apply load and properly create the P vector. Seems to work correctly, but should really consider creating some tests, and should also create an ideal model that would allow to verify the solver results (which is more about the stiffness assembly, really...). 

#### 11 Mar 24
Began adding commands to load nodes. Currently extremely tired so code is buggy and I am unable to debug it properly.

#### 6 to 10 Mar 24
This week has been quite tough for me at work. I have been reading the chapters from Felippa's Nonlinear Finite Element Analysis (NLFEA)book at a rate of about an hour/day or a bit more. However, I am finding the implementation of geometrically nonlinear analysis to be a little bit difficult. I have also gone through Bassam Izzuddin's Nonlinear Structural Analysis (NLSA), and found them to be a bit easier to follow but leave out a lot for interpretation. Going through the past exams, I have tried solving some of the problems by hand which I found to be very helpful.

Here's a summary of my reading of Felippa's NLFEA:
- Chapter 13: Corotational Formulation Overview I\
The corotational formulation chapter introduces the concept of using local axes to separate the rigid body motions and element deformations from the rest of the element transformations. The most important idea in this chapter, however, is that this process allows us to decouple the rigid body motions from the rest of the calculation process by applying it as a filtering step. So, the "Corotational filter" would act as an intermediary step between Assembler, element library, and solver. There is a very nice figure for this on page 13-8. As per Felippa's book, and as per his usual style, the idea is that we want to develop as many element types as possible and make everything in between sort of an add-on to enable "element reuse". This is not central to `Blaze`, which is likely to only employ 10 or so elements. The problem with this chapter, and the next one, is that they quickly become overwhelming (see page 13-13 for example), and I can no longer tell what is needed and what is not. 
- Chapter 14: Corotational Formulation Overview II\
Skipped in favour of moving on to other parts of the book that may be more useful for now.
- Chapter 16: The CR Formulation: BE Plane Beam *(CR = Corotational, BE = Bernouli-Euler)*\
This element is exactly what I want to implement in `Blaze` for 2D analysis. It shares many features in common with the beam-column element taught by Izzuddin in his NLSA notes, except it also has the additional y-direction (z-direction in `Blaze`; shear) freedom which is missing from Izzuddin's notes. Very interestingly as I revisit the chapter to write these notes, Felippa's notes introduce the second partial derivatives of the local displacements with respect to the global displacements in equations 16.13 through 16.15. When I read this chapter I did not understand why. Now, after having gone through Izzuddin's notes again, I know why - it is because we will need those derivatives in order to compute the geometric stiffness from the formula $\frac{\partial\boldsymbol{d}}{\partial \boldsymbol{u}\partial \boldsymbol{u}}\boldsymbol{f}$, where $\boldsymbol{d}$, $\boldsymbol{u}$, and $\boldsymbol{f}$ are the element freedoms, global freedoms, and element forces respectively. Additionally, Felippa has expressions for the element "stresses", internal strain energy, material stiffness, force vector, and both material and geometric stiffnesses. What is exceptionally nice is that Felippa divides the parts that depend on axial, bending, shear, and pure geometry into their separate components allowing an easier understanding of the matrix. The corotational transformations are embedded in all of these equations. Izzuddin's notes contain something very similar for the beam-column element (which I think is exactly the same) there. This will allow me to cross reference the two and understand Felippa's implementation in the same way the Rosetta stone allowed cross-translation of different languages.
- Chapter 20: Overview of Solution Methods\
Solution methods is one of the most important chapters to go through. It starts by first covering that the hierarchy of NLFEA solution procedure is threefold: stages > increments > iterations. Stages are like solution stages and are often completely forgotten about (Izzuddin's notes never mentioned them). What's really important is *incremental steps* or *increments* for short, and the iterations. Increments allow us to perform the solution procedure in steps rather than applying the full load. Thus, we can trace the equilibrium path of the structure by controlling the state vector $\boldsymbol{u}$ and the stage control parameter $\kappa$, which when together are called the *state control pair*. There are two primary ways to "march along" (really, better called trace equilibrium path): 1. Predictor-only (purely incremental) methods, and 2. Predictor-corrector (corrective) methods. The *forward-Euler method* $\Delta \boldsymbol{u}_n^0 = \boldsymbol{K}_n^{-1}\boldsymbol{q}_n^0=\bold{v}_n\Delta\kappa_n^0$ where the subscript refers to the increment (step) and the superscript refers to the iteration, and $\boldsymbol{q}$ is the *incremental load vector*. $\bold{v}$ is the incremental velocity vector. In a predictor-only method, the solution would drift off from the correct solution as we get away from equilibrium due to the fact that there will be a residual vector. 
- Izzuddin's NLSA Notes\
Bassam Izzuddin's notes are a very important source for me as they contain the description of NLFEA as I first learned it. Izzuddin is also the developer of `Adaptic`, a finite element program used at Imperial College for blast and fire, which makes it especially topical. Izzuddin's notes highlight the importance of *work*, and use clear notation to describe NLFEA in terms of it. Perhaps it is best to number my notes below rather than explain them in narrative style:
1. $\boldsymbol{T}$ is a nonlinear function of system state (global displacements) $\boldsymbol{U}$. That is: $\boldsymbol{T} = \partial \boldsymbol{d}/\partial\boldsymbol{U}$ or $T_{ij} = \partial d_i/ \partial U_j$. Note that $\boldsymbol{d}$ is a nonlinear function of $\boldsymbol{U}$
2. Straight equilibrium calculations are difficult to use due to constraints so instead we use the conservation of energy to achieve equilibrium: "A structural system is in equilibrium in its deflected configuration if the external work performed by teh applied loading over any possible infinitesimal displacement mode is equal to the internal work performed by the component forces over the corresponding compatible infinitesimal deformations."
3. We usually assemble the tangent stiffness from the element level tangent stiffnesses rather than assemble it using global equations.
4. The forces $\boldsymbol{P}$ would represent a work-conjugate load: $\boldsymbol{P} = \partial W/\partial\boldsymbol{U}$. Forces corresponding to constrained freedoms are not included here, as these would be included directly in the nonlinear relationship between $\boldsymbol{d}$ and $\boldsymbol{U}$. What I think Izzuddin is saying is that constraints are directly applied via the transformation, which is exactly what I have been doing.
5. The most important equation presented is likely the one formalising the concept of stationary potential energy $V$. It states that the change in total potential energy is equal to the change in internal strain energy minus the change in potential energy due to change in work, which should equal zero: 
$$\frac{\partial V}{\partial \boldsymbol{U}} = \left[\frac{\partial \boldsymbol{d}}{\partial \boldsymbol{U}}\right]^T \boldsymbol{f} - \frac{\partial W}{\partial \boldsymbol{U}} = \boldsymbol{0}$$
6. The following definitions are made also for equilibrium, where $\boldsymbol{G}$ is the out-of-balance force, $\boldsymbol{R}$ the *resistance forces*, and where $\boldsymbol{K}_T$ is the tangent stiffness matrix:
$$\boldsymbol{G} = \frac{\partial V}{\partial \boldsymbol{U}} = \boldsymbol{0}$$
$$\boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}$$
$$\boldsymbol{K}_T = \frac{\partial \boldsymbol{G}}{\partial \boldsymbol{U}}; K_{T,ij} = \partial G_i/\partial U_j$$

7. The first term of the geometric stiffness $\boldsymbol{K}_G$ is related to the forces of the strained system, while the second part is proportional to the applied work-conjugate loading:
  $$\boldsymbol{K}_G=\frac{\partial^2 \boldsymbol{d}}{\partial \boldsymbol{U} \partial{\boldsymbol{U}}}\boldsymbol{f} - \frac{\partial^2 W}{\partial\boldsymbol{U}\partial\boldsymbol{U}}$$
8. The tangent stiffness matrix is symmetric if (1) the system is elastic, and (2) the loading is conservative. Also, $\boldsymbol{d}$ is uniquely defined by $\boldsymbol{U}$, and so is $W$. Thus, $\boldsymbol{K}_G$ is also symmetric.
9. We define the loading in terms of the load factor $\lambda$, which was the system control parameter $\kappa$ in Felippa's notes. So, the loading is: $\boldsymbol{P} = \boldsymbol{P}_i + \lambda \boldsymbol{P}_n$
10. The tangent stiffness matrix is positive definite if, when transformed into an upper triangular matrix using Gaussian elimination, all its diagonals are positive. This indicates a stable structure. If any diagonals are negative, then the structure is unstable, and if none are negative but there are zeroes then the matrix is critically stable or critically unstable. If $(\partial \boldsymbol{U}/\partial\lambda \boldsymbol{P}_n < 0)$ then the system is unstable, but this does not need to be the case. 
11. For solving our nonlinear problems, we compute the out of balance $\boldsymbol{G}$ from $\boldsymbol{G} = \boldsymbol{R} - \boldsymbol{P}_i - \lambda\boldsymbol{P}_n$, where $\boldsymbol{R} = \partial U/\partial \boldsymbol{U}$ ($U$ being the internal strain energy). The out of balance forces due to the change in state are given by: 
$$\delta \boldsymbol{G} = \boldsymbol{K}_T \delta \boldsymbol{U} = \frac{\partial^2 U}{\partial \boldsymbol{U}} \partial{\boldsymbol{U}} \delta \boldsymbol{U}$$
12. Let us do load-control as per the final exam of NLSA from 2008:
    1.  $\boldsymbol{P} = \lambda \boldsymbol{P}_n = \lambda (\partial W_n/\partial \boldsymbol{U})$ where we are careful to capture the load from the expression of work as the load is dependent on the state. However, in many FEM systems, as likely will be the case for `Blaze`, the load will be independent of state (top of page 3 from Izzuddin's Nonlinear Solution Procedures for Traversing Equilibrium Paths)
    2.  Calculate out-of-balance forces: $\boldsymbol{G} = \boldsymbol{R} - \boldsymbol{P}$. If this is the first iteration of an unstressed system, then we start with resistance forces equal to 0: $\boldsymbol{G} = - \boldsymbol{P}$.
    3.  Calculate the change in state that results from the out-of-balance forces: $\Delta \boldsymbol{U} = - \boldsymbol{K}_T^{-1} \boldsymbol{G}$
    4.  Update the state based on the calculate "trial" change: $\boldsymbol{U} += \Delta \boldsymbol{U}$
    5.  Check the new out-of-balance forces again by first calculating the resistance forces as $\boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}$ for a newly calculated $\boldsymbol{T}$ and $\boldsymbol{f}$ based on the updated $\boldsymbol{U}$. Repeat steps 1 and 2 from above to recalculate $\boldsymbol{P}$ and $\boldsymbol{G}$, check out-of-balance for convergence, and then repeat steps 3, 4, and 5 if necessary. 
13. I could implement the procedure in step 12 above right now, if needed, but note that it requires inversing the tangent stiffness matrix which I am not sure is an *efficient* way to do this. The *Newton-Raphson* method requires we recalculate the tangent stiffness every iteration, while the *Modified Newton Raphson* does not. What is it that I need to add to `Blaze` to implement the basic load-controlled solver above?
    1.  Load application and restraint for nodes.
    2.  Nonlinear transformation $\boldsymbol{T}$ matrix for the element level.
    3.  Geometric and tangent stiffness matrix for the beam element.
    4.  Iterative solver.
14. I think I should probably read Frank McKenna's thesis and early papers to follow how he implements the solution procedure in OpenSees. **I should pause reading Felippa's book for the time being as I have enough for a super basic load-stepping solver.**

**Continue:**
- Chapter 20: Overview of Solution Methods (from page 20-9)
- Chapter 21: Continuation Under Load Control

**On hold:**
- Chapter 22: Continuation Under General Control: Predictor
- Chapter 23: Continuation Under General Control: Corrector
- Chapter 24: Incremental-Corrective Methods: Implementation *(there are mistakenly two chapters labelled as 24)*


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
- [x] The assembler should have a way to map nodes to the load and displacement vectors. Perhaps using `std::map`?
- [x] Nodes should have a container for loads, reaction forces, and for nodal dispalcements. 
- [x] Add functionality to add nodal load.
- [ ] The command to load a node adds the loaded node to a `std::set` of loaded nodes so each time we want to map the loaded nodes, we only need to loop over the loaded nodes.
- [x] The assembler should fill the load vector based on nodal loads.
- [x] After each successful analysis step, nodal displacements are updated by the assembler which maps the nodal displacement vector back to the nodes.
- [x] The solution procedure should include calculating the element internal forces and strains, and each element should have these saved.
- [ ] Fixing a node adds it to a `std::set` that corresponds to fixed nodes. After each successful analysis step, these nodes calculate their reaction forces. All nodes that do not have a constraint just have reaction forces of zero.
- [ ] Add loggers to retrieve and log certain displacements or reaction forces from the nodes, and internal forces/stresses/strains of elements.
