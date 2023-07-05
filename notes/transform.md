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
