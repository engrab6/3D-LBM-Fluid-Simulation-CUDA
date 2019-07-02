==================================
D3Q27 Lattice Boltzmann BGK Method
==================================

:author: CSE 6230
:date:   Fall 2017

Definition
----------

We consider a periodic 3D lattice of points :math:`(i,j,k)` for :math:`0\leq
i,j,k < N`.  At every such point there are densities :math:`f_\alpha` of flow
in 27 quantized velocities :math:`v_\alpha`, where :math:`\alpha \in
\{-1,0,1\}\times \{-1,0,1\} \times \{-1,0,1\}` corresponds to either the point
itself (:math:`\alpha = (0,0,0)`) or one of its neighbors.  The following steps going into update the densities in a timestep of a *Lattice Boltzmann Method*.

First the bulk characteristics of the velocity densities are calculated:

.. math::

    \begin{aligned}
      \rho &= \sum_{\alpha} f_\alpha & [\text{density}], \\
      u_x  &= \frac{1}{\rho}\sum_{j,k} f_{1,j,k} - f_{-1,j,k} &
      [x\text{ velocity}], \\
      u_y  &= \frac{1}{\rho}\sum_{i,k} f_{i,1,k} - f_{i,-1,k} &
      [y\text{ velocity}],\\
      u_z  &= \frac{1}{\rho}\sum_{i,j} f_{i,j,1} - f_{i,j,-1} &
      [z\text{ velocity}].
    \end{aligned}

These velocities are then spread between the quantized directions:

.. math::

    u_{(i,j,k)} = i * u_x + j * u_y + k * u_z.

These values determine the equilibrium densities :math:`f^0_{\alpha}`:

.. math::

    f^0_{\alpha} = \rho w_{\alpha} (1 + 3 u_{\alpha} + \frac{9}{2} u_{\alpha}^2 - \frac{3}{2} (u_x^2 + u_y^2 + u_z^2)),

Finally, the new densities at the point downstream of :math:`v_\alpha` are set to a weighted average of the original and equilibrium densities:

.. math::

    f_{\alpha}((i,j,k) + \alpha,t + 1) = \frac{1}{\tau} f^0_{\alpha}((i,j,k),t) + (1 - \frac{1}{\tau})f_{\alpha}((i,j,k),t)

The parameters :math:`w_{\alpha}` and :math:`\tau` determine the
characteristics of the flow.
