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



  
- Vector form of transport equations

  Rewriting the previously drived non-dimensionalized continuity and momentum equation in vector form generates a simple format that eases implementation of the numerical method. The above transport equation can be newly formed as shown below:

  .. math::

     \frac{\partial \vec{U}}{\partial t} + \frac{\partial \vec{E}}{\partial x} + \frac{\partial \vec{F}}{\partial y} = \frac{1}{\text{Re}} \left ( \frac{\partial^{2}}{\partial x^{2}} + \frac{\partial^{2}}{\partial y^{2}} \right ) \vec{U}

  where the each of vector elements are summarized below:

  .. math::

     \vec{U} = \begin{bmatrix}P\\ u\\ v \end{bmatrix}, \;\; \vec{E} = \begin{bmatrix} \frac{u}{\beta}\\ uu + P\\ uv\end{bmatrix}, \;\; \vec{F} = \begin{bmatrix} \frac{v}{\beta}\\ uv\\ vv + P\end{bmatrix}

  Now this is good to go further for descritization because the given task is to solve explicit form of discretization equation. Even though the derived form of transport equation is not linearized, each of vectors above are easily discretized in terms of their elements that are combinations of each primitive variables. Thus, in this project, actual discretization has been doen form the driven transport equation above.

- Finding time step algorithm

  In order to find time step that may stabilize the numerical solution, we need to know system convecting velocity as we pick the coefficient of spatial derivative terms in Burger's and Euler equations as the convection velocity. The driven system of equation is not a single partial different equation but a set of three different partial different equation. To find the convection speed of numerical information in the time and space domains, we need to first linearize the given system of equations and find the Eigen values. The linearization can be obatained by following process. The driven system of PDE should be reformulated in linearized set of equations:

  .. math::

     \frac{\partial \vec{U}}{\partial t}  + \left [ A \right ] \frac{\partial \vec{U}}{\partial x} + \left [ B \right ] \frac{\partial \vec{U}}{\partial y} = \frac{1}{\text{Re}} \left ( \frac{\partial^{2}}{\partial x^{2}} + \frac{\partial^{2}}{\partial y^{2}} \right ) \vec{U}

  Now we have found two coefficient matrices of convection terms and the spatial derivatives is now taken with respect to :math:`\vec{U}` only. Despite the vector form, the PDE form is a identical with Burger's equation. The coefficient matrices are below listed:

  .. math::

     \left [ A \right ] = \begin{bmatrix} 0 & \frac{1}{\beta} & 0 \\ 1 & 2u & 0\\ 0 & v & u \end{bmatrix}, \;\; \left [ B \right ] = \begin{bmatrix} 0 & 0 & \frac{1}{\beta} \\ 0 & v & u\\ 1 & 0 & 2v \end{bmatrix}

  The resolved Eigen values of :math:`\left [ A \right ]` and :math:`\left [ B \right ]` matrices are :math:`u, u+a, u-a` and :math:`v, v+a, v-a`, respectively. Taking :math:`\left [ A \right ]` for example, the maximum convection velocity that transmit the numerical information can then be :math:`\left | u  \right | + a`. Therefore, the Courant number for this case can also be determined by:

  .. math::

     


     
