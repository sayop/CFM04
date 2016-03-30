=================
 Code Instruction
=================

The present project is aimed to develop a computer program for solving a steady solution with projection correction method. The code being used for answering all the question here is written with Python language. The current version of Python code is an extended version from the previous homework. Therefore, the current Python code has also a feature of pressure correction method as well as artificial compressibility method. This program is to run with simple command::
 
  $ python main.py

Quick instruction for running the simulation
--------------------------------------------

The Python code used for this project can be cloned from *github.com* repository::

  $ git clone https://github.com/sayop/CFM04

You can also see the code directly by visiting the website: https://github.com/sayop/CFM04 If you clone the code, you will see the following set of files and directories::

  $ sayop@reynolds:~$ ls CFM04/
  docs  README.md  src

*docs* contains the document files set for the current project using *Sphinx* software. This *pdf* document is online available at: http://cfm04-gatech.readthedocs.org. The Python script for this simulation is stored in *src* folder.

Before running the simulation, you need to open the file named *input.in* using editor for example, VI on unix system::
 
  $ vi input.in

Then, you should be able to see the following set of simulation parameters::

  #grid dimension
  iDim            20
  jDim            20
  xmin            0
  xmax            10
  ymin            0
  ymax            10
  #boundary conditions
  uLeft           0.0
  vLeft           0.0
  uRight          0.0
  vRight          0.0
  uBottom         0.0
  vBottom         0.0
  uUp             10.0
  vUp             0.0
  # fluid kinematic properties
  nu              1.0
  pInit           10.0
  #simulation setup
  pCorr           1
  alpha           1.0
  pResidual       0.01
  maxIter         100000000000
  Courant         0.5
  dtInit          0.0001
  Beta            0.5
  residualMin     0.00005
  #Post-Process
  nIterWrite      100


The parameter's name above will literally tell you what every single variables indicates in the simulation. For the post-processing as requested in this project, *nIterWrite* will write a solution plot and CSV file at speicifed interval of time integration number.

The most important feature here is to set the pressure correction method. In order to have this goal, you will need to set *pCorr* to 1 to switch on its feature. Otherwise, you will run your simulation with artificial compressibility method. 

Also, to set the specified Reynolds number, you need to change the *uUp* that will show a correponding Reynolds number at the beginning of your simulation on screen. Here, 10.0 of *uUp* will maintain the Reynolds number 100.0. All the input parameters are dimensional quantities. And these variables will be non-dimensionalized when they are transitioned to the main loop of the simulation.
