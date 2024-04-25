# Notes
Contains notes and some other things I want to remind myself of and I want access to in one place.

# Transformations
## Translations
### Graphical
There is no rotational DoF, so we don't need to worry about the effects of `element` rotation on translation.
$$

   \begin{bmatrix}
   u'_1 \\
   v'_1 \\
   1 \\
   \end{bmatrix}
   =
   \begin{bmatrix}
   1 & 0 & \Delta_x \\
   0 & 1 & \Delta_y \\
   0 & 0 & 1 \\
   \end{bmatrix}

   \begin{bmatrix}
   u_1 \\
   v_1 \\
   1 \\
   \end{bmatrix}
   =

   \begin{bmatrix}
   u_1 + \Delta_x\\
   v_1 + \Delta_y \\
   1 \\
   \end{bmatrix}
$$
### Physical (FEM)
We have an angular DoF, so our main concern is that the angular DoF will contribute to the lateral displacements.
$$

   \begin{bmatrix}
   u'_1 \\
   v'_1 \\
   \theta'_{z,1} \\
   \end{bmatrix}
   =
   \begin{bmatrix}
   1 & 0 & \Delta_y \\
   0 & 1 & 0 \\
   0 & 0 & 1 \\
   \end{bmatrix}

   \begin{bmatrix}
   u_1 \\
   v_1 \\
   \theta_{z,1} \\
   \end{bmatrix}
   =

   \begin{bmatrix}
   u_1 + \theta_{z,1}\Delta_y\\
   v_1 \\
   \theta_{z,1} \\
   \end{bmatrix}
$$
Where we have used the following small angle assumptions regarding $\theta_{z,1}$:
$$
\sin(\theta\rightarrow0) = \theta
\\
\cos(\theta\rightarrow0) = 1
$$
And we had considered that the change in $u'_1$ due to rotation is $\Delta_y\sin(\theta_{z,1})$, and in $v'_1$ is $\Delta_y\cos(\theta_{z,1}) - \Delta_y$. The caveat here is that this only works for small rotation angles, not for large ones: <span style="color:red;">this will not work very well for large deformations!</span>

When we combine the translation transformation with the usual rotation transformation, we get:
$$
   \begin{bmatrix}
   1 & 0 & \Delta_y \\
   0 & 1 & 0 \\
   0 & 0 & 1 \\
   \end{bmatrix}
   \begin{bmatrix}
   c & s & 0 \\
   -s & c & 0 \\
   0 & 0 & 1 \\
   \end{bmatrix}
   =
   \begin{bmatrix}
   c & s & \Delta_y \\
   -s & c & 0 \\
   0 & 0 & 1 \\
   \end{bmatrix}
$$
The order is **very important** here. We have to recall that the offset is applied to the $\theta'_{z,1}$ local rotation. That means that we need to first get that rotation from the global DoFs whatever they may be. As such, we must first perform the rotational transformation followed by the translation transformation. As in, we first multiply $\{d\}$ by $[T_{\theta}]$ then by $[T_{\Delta}]$: $\{d'\} = [T_\Delta][T_\theta]\{d\}$.



## Review of some other software
- [Elmer uses Metis](https://youtu.be/84K6OxEKEjQ?t=1358) for graph partitioning to discretise the domain for the processes. Elmer is open source and its code can be found [here](https://github.com/ElmerCSC/elmerfem). 
- FEniCS has an entire [book](https://launchpadlibrarian.net/83776282/fenics-book-2011-10-27-final.pdf). This book also includes information about what algorithms are used, and how they implemented in parallel.
- Code Aster appears to be less well developed and maintained than the other two, but has been [used for structures in fire](https://www.code-aster.de/project-cases/project-cases-detail/analysis-of-steel-reinforced-concrete-exposed-to-fire-mfpa-leipzig-gmbh-copy.html). The website has much more of a commercial bend to it with the product being seminars and tutorials.
- [Calculix](http://www.dhondt.de/) uses [PreProMax](https://prepomax.fs.um.si/), which is an open source pre and post processor. From looking around the internet it may be that this program was not really made for HPC.
## Some tips
- You need to specify the friends of class A in class A, so that the friend class can access class A members. "John is my friend and what's mine is his" vs. ~~"I am friends with Alex so I will access his stuff!"~~
