import numpy as np
from solutionMethod import *
from pCorrection import pressureCorrect
from post import *
import time

def timeIntegrate(inputDict):
   tStart  = 0.0
   Cr      = float(inputDict['Courant'])
   imax    = int(inputDict['iDim'])
   jmax    = int(inputDict['jDim'])
   maxIter = int(inputDict['maxIter'])
   pCorr   = int(inputDict['pCorr'])
   beta    = float(inputDict['Beta'])
   residualMin = float(inputDict['residualMin'])
   nIterWrite  = int(inputDict['nIterWrite'])


   # initialize residual variables for p, u, and v
   residualInit = np.zeros(3)
   residual     = np.zeros(3)
   UresidualLog  = []

   # start to count time for calculting computation performance
   start = time.clock()

   # Non-dimensionalize flow and domain variables
   nondimensionalize(inputDict, 1, 1, 1)

   #
   # Time Marching:
   #
   print '=============================================='
   print '# Time integration starts at t = %s' % tStart
   print '=============================================='
   t = tStart
   nIter = 0
   if pCorr != 1:
      nSubIter = 1
   elif pCorr == 1:
      # If pressure correction method, sub-iteration needs to update u^* and u^(n+1), sequentially.
      nSubIter = 2
   
   while True:
      nIter += 1
      # populate flux vectors' elements with updated u, v, p
      populateFluxVectors(pCorr,inputDict)
      
      # find dt with stability condition
      dt = updateTimeStep(inputDict,nIter)

      for n in range(nSubIter):
         if n == 1:
            # Second step for projection method
            # pressure correction: Run laplacian equation for pressure
            pressureCorrect(inputDict, dt)

         # if it is projection method, u* and v* will be updated at first and
         # u^(n+1) and v^(n+1) will be later updated at second step.
         # Q vector represents RHS terms when time-integrating U vector.
         # update Q vector for explicit time integration
         updateQvector(inputDict, dt, n, pCorr)

         # update new time level of p, u, v in primitive variable form
         updatePrimitiveVars(pCorr,imax,jmax,dt)

      # update boundary condition for pressure only
      updatePressureBC(imax, jmax)

      # compute residual value to verify convergence
      residual = computeResidual(imax, jmax, dt, FDM.Q)
      if nIter == 1:
         residualInit = residual
      # trace residual for u-velocity
      resNorm = residual[1] / residualInit[1]
      UresidualLog.append(resNorm)

      # dimensionless time increment because all the time variables have been non-dimensionalized above.
      t += dt
      if pCorr != 1:
         MachX, MachY = computeMaximumMach(imax, jmax, beta)
         print "|- nIter = %s" % nIter, ", t = %.6f" % t, ", dt = %.6f" % dt, ", Maximum Mach_x = %.4f" % MachX, ", Maximum Mach_y = %.4f" % MachY, ", u-residual = %.5f" % resNorm
      elif pCorr == 1:
         print "|- nIter = %s" % nIter, ", t = %.6f" % t, ", dt = %.6f" % dt, ", u-residual = %.5f" % resNorm

      if (nIter % nIterWrite == 0):
         dimensionalize(inputDict, 1, 1, 0)
         plotStreamLine(domainVars.x, domainVars.y, flowVars.u, flowVars.v, nIter)
         plotContour(domainVars.x, domainVars.y, flowVars.u, flowVars.v, nIter)
         nondimensionalize(inputDict, 1, 1, 0)

      if (nIter >= maxIter or resNorm <= residualMin): break

   #
   # time elapsed:
   elapsedTime = (time.clock() - start)
   print "## Elapsed time: ", elapsedTime

   # Dimensionalize flow and domain variables
   dimensionalize(inputDict, 1, 1, 1)

   # plot contour of artificial pressure
   #pltFile = 'p_contour.png'
   #phi = flowVars.p
   #phiMin = phi.min()
   #phiMax = phi.max()
   #plotContour(domainVars.x, domainVars.y, phi, phiMin, phiMax, pltFile)

   # plot contour of velocity magnitude
   #pltFile = 'Umag_contour.png'
   #flowVars.Umag = np.sqrt(flowVars.u ** 2 + flowVars.v ** 2)
   #phi = flowVars.Umag
   #phiMin = phi.min()
   #phiMax = phi.max()
   #plotContour(domainVars.x, domainVars.y, phi, phiMin, phiMax, pltFile)

   # plot streamline of velocity
   plotStreamLine(domainVars.x, domainVars.y, flowVars.u, flowVars.v, nIter)

   # trace center-line data to be compared to the Ghia's paper data
   nondimensionalize(inputDict, 1, 1, 0)
   traceCenterLineData('x', flowVars.v, 'v-velocity_in_x.csv')
   traceCenterLineData('y', flowVars.u, 'u-velocity_in_y.csv')
   dimensionalize(inputDict, 1, 1, 0)

   # write a log file for u-residual
   csvFile = 'u-residualLog.csv'
   logResidual(UresidualLog, csvFile)
