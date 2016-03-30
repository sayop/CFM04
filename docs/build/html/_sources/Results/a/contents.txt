=============
 Problem1 - a
=============

Describe the essential steps of the solution method. Include the discretized equations and implementation of boundary conditions.

In this problem set, we are supposed to solve the Navier-Stokes equations having continuity and momentum conservation equations together. Tensor forms of continuity and momentum equations are given below:

- Continuity (incompressible)

  .. math::

     \frac{\partial u_{i}}{\partial x_{i}} = 0 

- Momentum equation:

  .. math::

     \frac{\partial u_{i}}{\partial t} + \frac{\partial u_{i}u_{j}}{\partial x_{j}} = -\frac{1}{\rho}\frac{\partial p}{\partial x_{i}} + \nu \frac{\partial}{\partial x_{j}}\left ( \frac{\partial u_{i}}{\partial x_{j}} \right )


- Non-dimensionalization of the Navier-Stokes equations

  In some cases, it is beneficial to non-dimensionalize the given transport equation because it eases the analysis of problem of interest, and also may reduce the number of parameters. The non-dimensionalized form of the Navier-Stokes equation can be achieved by first normalizing the primitive variables as followings:

  .. math::

     \tilde{u_{i}} = \frac{u_{i}}{U_{\text{ref}}},\;\;  \tilde{x_{i}} = \frac{x_{i}}{L_{\text{ref}}},\;\; \tilde{\rho}=\frac{\rho}{\rho_{\text{ref}}},\;\;\tilde{P} = \frac{P}{\rho_{\text{ref}}\, U^{2}_{\text{ref}}},\;\; \tilde{t}=\frac{t}{L/U_{\text{ref}}}

  For the final form of non-dimensionalized Navier-Stokes equation, tilda, :math:`\tilde{}`, will be dropped out for brevity and a new non-dimensional physical parameter :math:`Re` that represents the flow intertia against the fluid viscosity is introduced. Now we got:

  .. math::

     \frac{\partial u_{i}}{\partial t} + \frac{\partial u_{i}u_{j}}{\partial x_{j}} = -\frac{1}{\rho}\frac{\partial p}{\partial x_{i}} + \frac{1}{\text{Re}} \frac{\partial}{\partial x_{j}}\left ( \frac{\partial u_{i}}{\partial x_{j}} \right )

  where the Reynolds number is defined as:

  .. math::

     \text{Re} = \frac{U_{\text{ref}}L_{\text{ref}}}{\nu}


- Projection method

  Given partial differention equation set is composed of continuity and momentum conservation equations. By contrary to the artificial compressiblity method, the projection method is to solve this equation set with incompressible solution method, so-called, pressure-based method. In this type of method, the continuity equation is not directly solved but used for divergence-free constraint, which is an outcome of incompressible flow continuity feature.

  The momentum equation given in the tensor form can be recast in terms of time differential form as shown below. And the convective and diffusion flux terms are combined together to express a simple form of spatially differenciated flow quantities in terms of :math:`u` and :math:`v`. And the pressure derivative term is left alone for a particular purpose of pressure correction.


  .. math::

     \frac{\partial u_{i}}{\partial t} = H_{i} -\frac{\partial p}{\partial x_{i}}

  where the convection terms are represented by :math:`H_{i}`:

  .. math::

     H_{i} = -\frac{\partial u_{i}u_{j}}{\partial x_{j}} + \frac{1}{\text{Re}}\frac{\partial^{2} u_{i}}{\partial x_{j} \partial x_{j}}


  Here, the above expression can be again discretized in two separate time steps having different right hand side quantities. Two separate the temporal difference, we employs an intermediate time indicated by :math:`*`. The intermediate time level of velocity is updated without pressure derivative term. So this time level is incomplete so it is not divergence free.

  .. math::

     \frac{u^{*}_{i} - u^{n}_{i}}{\Delta t} = H_{i}

  .. math::

     \frac{u^{n+1}_{i} - u^{*}_{i}}{\Delta t} = -\frac{\partial p^{n+1}}{\partial x_{i}}

  Summing above two equations gives a complete form of temporal difference equation of momentum conservation equation. Then to find the divergence free and next time level of velocity components, we need to first evaluate the intermediate velocity quantities and pressure at the next time level. But, the given form of equation is not well posed because we haven't yet solved next time level of pressure which also satisfies the divergence free constraint. To do this, we need to take divergence of the equation then it gives a elliptic equation form what is called Poisson's equation.

  .. math::

     \frac{1}{\Delta t} \frac{\partial u^{*}}{\partial x_{i}} = \frac{\partial^{2} p^{n+1}}{\partial x_{i} \partial x_{i}}

  Since the equation is characterized with elliptic equation so it becomes boundary condition problem in the given domain. Because of two dimensional space and non-linearity, it needs to be solved with point-iterative method such as Jacobi method and Gauss-Seidel method. In this project, this solution for pressure is resolved with succesive over-relaxation method to have faster convergence.

  
- Vector form of transport equations

  Rewriting the previously drived non-dimensionalized momentum equations in vector form generates a simple format that eases implementation of the numerical method. The above transport equation can be newly formed as shown below:

  .. math::

     \frac{\partial \vec{U}}{\partial t} + \frac{\partial \vec{E}}{\partial x} + \frac{\partial \vec{F}}{\partial y} + \vec{\triangledown }p = \frac{1}{\text{Re}} \left ( \frac{\partial^{2}}{\partial x^{2}} + \frac{\partial^{2}}{\partial y^{2}} \right ) \vec{U}

  where the each of vector elements are summarized below:

  .. math::

     \vec{U} = \begin{bmatrix}u\\ v \end{bmatrix}, \;\; \vec{E} = \begin{bmatrix} uu \\ uv\end{bmatrix}, \;\; \vec{F} = \begin{bmatrix} uv\\ vv \end{bmatrix}, \;\; \vec{\triangledown} = \begin{bmatrix} \frac{\partial}{\partial x}\\ \frac{\partial}{\partial y} \end{bmatrix}

  Compared to the previous homework that constructs the flux vector including pressure force, the current flux vectors :math:`E` and :math:`F` contain only convection terms in x and y directions. This is the purpose that the numerical method separates two different steps: evaluating intermediates velocities and projection of this velocities onto the divergence free satisfied velocity field.

  Now this is good to go further for descritization because the given task is to solve explicit form of discretization equation. Even though the derived form of transport equation is not linearized, each of vectors above are easily discretized in terms of their elements that are combinations of each primitive variables. Thus, in this project, actual discretization has been doen form the driven transport equation above.

- Finding time step algorithm

  Contrary to the previous homework, the system of partial differential equation is composed of only momentum conservation equations in x and y direction, that is, this system does not contain continuity equation. Note that the mass conservation only applies as a divergence free constraint. So the system convection velocity is straighforward rather than having to find eigenvalues of coefficient matrices. So the Courant number relation for the two-dimensional system of equation can be defined by:

  .. math::

     \text{Courant} = \frac{u dt}{dx} + \frac{v dt}{dy}

  The numerical time step is evaluted as a minimum time step that satifies the given relation at every computational node points.
