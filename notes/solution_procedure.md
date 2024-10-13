# Solution Procedure
## definitions
- $\boldsymbol{P}$ a sparse dynamic `Eigen` matrix that captures the nodal loads. Of size $n_{DoFs}$.
- $\boldsymbol{R}$ a sparse dynamic `Eigen` matrix that captures the nodal loads. Of size $n_{DoFs}$.
- $\boldsymbol{G}$ sparse dynamic `Eigen` matrix for storing out of balance forces: $\boldsymbol{G} = \boldsymbol{R} - \boldsymbol{P}$. Of size $n_{DoFs}$.
- $\boldsymbol{U}$ a sparse `Eigen` vector for the nodal displacements. Of size $n_{DoFs}$.
- $d\boldsymbol{U}$ a sparse `Eigen` vector for the increment in nodal displacements. Of size $n_{DoFs}$.
- $\boldsymbol{K}$ a sparse dynamic `Eigen` matrix for the global tangent stiffness of the model. Of size $n_{DoFs} \times n_{DoFs}$.
## initialisation of solution parameters
- Set `load_factor` to 0 and `step` to 1.
- Define a uniform `dLF` equal to $dLF = LF_{max}/n_{steps}$.
- Define a uniform maximum number of iterations.
## solution loop
- LF Loop: while $LF < LF_{max}$:
  - Increment `load_factor` by `dLF`, and increment all loads by `dLF`
  - Increment all loads by `dLF`: load manager iterates over its `nodal_loads` vector of `NodalLoad` objects. Each of these objects tells its loaded node to increment its `nodal_loads` (in the `Node` object) for each loaded DoF by `dLF * nodal_loads[dof]` where `nodal_loads` is a container handled by the `NodalLoad` object.
  -  Set the iteration counter `iter` to 1, and convergence boolean `converged` to false.
  -  Iterations Loop: while `iter` < `max_iter` and `converged == false`
     - Loop over $n_{nodes}$ to map global displacements vector $\boldsymbol{U}$ to each `Node` object's `std::array<real, 6> nodal_displacements`. <span style="color:lightgray">The mapping is done from the `U` vector over the node objects container `std::vector<std::shared_ptr<Node>> node_vector` in `GlobalMesh`. This is a one-to-one mapping where each nodal DoF corresponds directly to the the DoFs in the $\boldsymbol{U}$ vector. In this loop over, each node retrieves its own displacements for its active DoFs.</span>
     - Loop over $n_{elements}$ to update the state of each element. For `Nonlinear2DPlasticBeamElement` this procedure is as follows:
     - Loop over both $n_{elements}$ and $n_{nodes}$ to ask elements to calculate the global stiffness triplets (contribution and location of element stiffness contribution) that will populate the global stiffness, and ask the nodes to compute the load triplets that correspond to the global force vector $\boldsymbol{P}$. <span style=color:red>The nodal contributions to the force vector $\boldsymbol{P}$ do not change during the iterations loop, only during the LF loop. Should not be doing this every iteration except for stiffness of elements which depend on iteration state.</span>
     - Loop over both $n_{elements}$ and $n_{nodes}$ to retrieve their triplets and assemble them into $\boldsymbol{P}$ and $\boldsymbol{K}$. <span style="color:red"> As before, there is no need to reassemble the nodal load vector $\boldsymbol{P}$ since it does *not* change in the iterations loop.</span> This operation also involves creating 4 temporary `std::vector` to keep the triplets, collecting all contributions and *remaking* a compressed $\boldsymbol{P}$ and $\boldsymbol{K}$.
     - Loop over $n_{elements}$ and assemble their resistance forces into $\boldsymbol{R}$ in a similar way to the previous step creating a temporary `std::vector` and remaking $\boldsymbol{R}$.
     - Calculate the L2 norm of the out-of-balance vector $\boldsymbol{G}$, and compare to tolerance thus setting the convergence state.
     - If the solution is not converged, solve the system of equations: $\boldsymbol{G} = \boldsymbol{K}d\boldsymbol{U}$ for $d\boldsymbol{U}$. This is done in three steps: `analyzePattern`, `factorize`, and then `solve` all of which are handled by `Eigen3`. 
     - Increase the iteration counter `iter`
  - if we are here, it means that either all iterations finished or the solution converged. If the solution converged, increase our step counter `step` and loop over $n_{elements}$ to update their sections' starting stress states.

---
### Plastic Element update procedure 
- Retrieve nodal displacements from the `Node` objects of the elements and store them in the `global_ele_U` vector of size 12 where we have 6 DoFs for the element $\left[ U^1_{1} \ U^1_{2} \ U^1_{33} \ U^2_{1} \ U^2_{2} \ U^2_{33}\right]^T$, and 6 extra DoFs for global compatibility

    <div style="border:1px solid black; padding: 10px;">
    
     - Retrieve nodal displacements from the `Node` objects of the elements and store them in the `global_ele_U` vector of size 12 where we have 6 DoFs for the element $\left[ U^1_{1} \ U^1_{2} \ U^1_{33} \ U^2_{1} \ U^2_{2} \ U^2_{33}\right]^T$, and 6 extra DoFs $\left[ U^1_{11} \ U^1_{22} \ U^1_{3} \ U^2_{11} \ U^2_{22} \ U^2_{3}\right]$ for global compatibility.
     - Calculate $\boldsymbol{d} = \left[ \Delta \ \theta_1\  \theta_2 \right]^T$ from the nonlinear corotational transformation object.
     - Retrieves the updated element length from the nonlinear transformation object.
     - For each Gauss Point, calculate the $\boldsymbol{B}$ matrix: $\boldsymbol{B} = \frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{d}^T} = \begin{bmatrix} \frac{1}{L_0} & \frac{2\theta_1}{15} - \frac{\theta_2}{30} & -\frac{\theta_1}{30} + \frac{2\theta_2}{15} \\ 0 & -\frac{4}{L_0} + \frac{6x}{L_0 ^2} & -\frac{2}{L_0} + \frac{6x}{L_0 ^2}\end{bmatrix}$ 
     - For each Gauss Point, calculate the strain vector $\boldsymbol{\varepsilon} = \begin{bmatrix} \frac{\Delta}{L_0} + \left(\frac{2\theta_1^2 - \theta_1 \theta_2 + 2 \theta_2^2}{30}\right) \\ \left(- \frac{4}{L_0} + \frac{6x}{L_0^2} \right)\theta_1 + \left(- \frac{2}{L_0} + \frac{6x}{L_0^2} \right)\theta_2 \end{bmatrix}$ </div>
