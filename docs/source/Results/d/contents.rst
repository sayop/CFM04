=============
 Problem1 - d
=============

Examine the method stability with different grids. Determine the maximum time step that leads to a stable solution and compare it to the stability criteria.


------------------
 Grid spacing test
------------------

In this project, variable time stepping method is employed such that temporal numerical instability is avoided. Thus, given conditions of Re=100 and Re=500 are numerically stable with relevant choice of Courant number. For this reason, the grid spacing test is performed with much higher Reynolds number condition such that physical viscousity effect is not significantly taken into account. Note that we already know Euler equation having no viscousity is unconditionally unstable with central finite difference method. 

On the other hand, having visous terms will attenuate the possible growth of numerical instability and leads to stable solution with properly chosen Courant number. All the preliminary tests with given Re conditions gave stable solutions. This is because viscous terms play a role of smoothing out the stiff solutions so as to result in numerically stable solution.

Here, the test has conducted with Re = 10,000 such that the equation almost goes to Euler equation. We can expect that numerical solution will be unconditionally unstable because our employed finite difference is done with central differencing. The test shown below proves this. On the other hand, the effect of grid spacing on the unstable solution's temporal evolution is clarified. Even though all the cases go to unstable solution after all, the onset of amplified numerical instability appears at different location in iteration number.

- Effect of grid spacing on the temporal evolution of axial velocity residual. (Re = 10,000)

  .. figure:: ./images/gridSpacing.png
     :scale: 80%


------------------
 Maximum time step
------------------

In this code, the variable time step method is used to maintain stable numerically. Therefore, the code does not run with constant time step. The maximum time step test is performed with different set of *Courant* number condition. The grid spacing is fixed with 20x20 to have fast running of simulation.


  .. figure:: ./images/timeStep.png
     :scale: 80%

  +-----------+-----------------------+
  | Courant # | dt at 100th iteration |
  +===========+=======================+
  |  0.5      | 0.038791              |
  +-----------+-----------------------+
  |  0.8      | 0.061655              |
  +-----------+-----------------------+
  |  1.0      | 0.064686              |
  +-----------+-----------------------+
  |  1.1      | 0.063479              |
  +-----------+-----------------------+
  |  1.2      | 0.072822              |
  +-----------+-----------------------+


