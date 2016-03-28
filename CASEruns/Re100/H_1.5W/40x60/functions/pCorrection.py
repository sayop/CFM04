import numpy as np
#from solutionMethod import centralFiniteDifference
from solutionMethod import *

def pressureCorrect(inputDict, dt):
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])
   alpha = float(inputDict['alpha'])
   convergeCrit = float(inputDict['pResidual'])

   # update Right Hand-Side terms
   RHS = np.zeros((imax,jmax))

   RHS += centralFiniteDifference(flowVars.u,'x',1,True) / dt
   RHS += centralFiniteDifference(flowVars.v,'y',1,True) / dt

   # Run point iteration method with Jacobi method
   iterMax = 1000 # maximum iteration for point-iterative method
   for niter in range(iterMax):
      #if niter == 0: residualInit = max(1e-99,computeResidual(flowVars.p, RHS, imax, jmax))
      # run with point-iteration method: will update new pressure at new time level
      #flowVars.p, residual = pointIterJacobi(flowVars.p, RHS, imax, jmax, alpha)
      # Jacobi method
      #flowVars.p, residual = pointIterJacobi(flowVars.p, RHS, imax, jmax, alpha)
      # Gauss-Seidel method
      flowVars.p, residual = pointIterGS(flowVars.p, RHS, imax, jmax, alpha)
      #updatePressureBC(imax, jmax)

      #residual = computeResidual(flowVars.p, RHS, imax, jmax)
      # update boundary condition for pressure only
      # compute residual at initial sate
      if niter == 0: residualInit = max(1e-99,residual)
      residual = residual / residualInit
      #print niter, residual
      #print flowVars.p
      if residual <= convergeCrit or niter == iterMax-1:
         print '| Pressure Correction: converged at residual = %.6f' % float(residual), ' and iterations of %s' % niter
         break

def pointIterGS(phi, RHS, imax, jmax, alpha):
   residual = 0.0
   dx = domainVars.dx
   dy = domainVars.dy
   dxSqr = dx * dx
   dySqr = dy * dy

   for j in range(jmax):
      #if j == 0: continue
      for i in range(imax):
         #if i == 0: continue
         # coefficient matrix
         # neighborL: values at points left to (i,j) point
         # neighborR: values at points right to (i,j) point
         # neighborN: values at points north to (i,j) point
         # neighborS: values at points south to (i,j) point
         #myCoef = 0.0
         #if i == 0:
         #   myCoef += 1.0 / dxSqr
         #   neighborL = 0.0
         #   neighborR = (phi[i+2,j] - 2.0 * phi[i+1,j]) / dxSqr
         #elif i == imax-1:
         #   myCoef += 1.0 / dxSqr
         #   neighborL = (phi[i-2,j] - 2.0 * phi[i-1,j]) / dxSqr
         #   neighborR = 0.0
         #else:
         #   myCoef += -2.0 / dxSqr
         #   neighborL = phi[i-1,j] / dxSqr
         #   neighborR = phi[i+1,j] / dxSqr

         #if j == 0:
         #   myCoef += 1.0 / dySqr
         #   neighborS = 0.0
         #   neighborN = (phi[i,j+2] - 2.0 * phi[i,j+1]) / dySqr
         #elif j == jmax-1:
         #   myCoef += 1.0 / dySqr
         #   neighborS = (phi[i,j-2] - 2.0 * phi[i,j-1]) / dySqr
         #   neighborN = 0.0
         #else:
         #   myCoef += -2.0 / dySqr
         #   neighborS = phi[i,j-1] / dySqr
         #   neighborN = phi[i,j+1] / dySqr
         myCoef    = -2.0 * (1.0 / dxSqr + 1.0 / dySqr)
         if i == 0:
            neighborL = 1.0 / dxSqr * phi[i+1,j]
         else:
            neighborL = 1.0 / dxSqr * phi[i-1,j]
         if i == imax-1:
            neighborR = 1.0 / dxSqr * phi[i-1,j]
         else:
            neighborR = 1.0 / dxSqr * phi[i+1,j]
         if j == 0:
            neighborS = 1.0 / dySqr * phi[i,j+1]
         else:
            neighborS = 1.0 / dySqr * phi[i,j-1]
         if j == jmax-1:
            neighborN = 1.0 / dySqr * phi[i,j-1]
         else:
            neighborN = 1.0 / dySqr * phi[i,j+1]

         tmp = (RHS[i,j] - neighborL - neighborR - neighborS - neighborN) / myCoef
         phiIncrement = alpha * (tmp - phi[i,j])
         residual += phiIncrement ** 2
         phi[i,j] = phi[i,j] + phiIncrement

   numberOfInnerPoints = imax * jmax
   residual = residual / numberOfInnerPoints
   residual = np.sqrt(residual)

   return phi, residual

def pointIterJacobi(phi, RHS, imax, jmax, alpha):
   newPhi = phi
   #newPhi = np.zeros((imax,jmax))
   residual = 0.0
   dx = domainVars.dx
   dy = domainVars.dy
   dxSqr = dx * dx
   dySqr = dy * dy

   for j in range(jmax):
      #if j == 0: continue
      for i in range(imax):
         #if i == 0: continue
         # coefficient matrix
         # neighborL: values at points left to (i,j) point
         # neighborR: values at points right to (i,j) point
         # neighborN: values at points north to (i,j) point
         # neighborS: values at points south to (i,j) point
         myCoef = 0.0
         if i == 0:
            myCoef += 1.0 / dxSqr
            neighborL = 0.0
            neighborR = 1.0 / dxSqr * (phi[i+2,j] - 2.0 * phi[i+1,j])
         elif i == imax-1:
            myCoef += 1.0 / dxSqr
            neighborL = 1.0 / dxSqr * (phi[i-2,j] - 2.0 * phi[i-1,j])
            neighborR = 0.0
         else:
            myCoef += -2.0 / dxSqr
            neighborL = 1.0 / dxSqr * phi[i-1,j]
            neighborR = 1.0 / dxSqr * phi[i+1,j]

         if j == 0:
            myCoef += 1.0 / dySqr
            neighborS = 0.0
            neighborN = 1.0 / dySqr * (phi[i,j+2] - 2.0 * phi[i,j+1])
         elif j == jmax-1:
            myCoef += 1.0 / dySqr
            neighborS = 1.0 / dySqr * (phi[i,j-2] - 2.0 * phi[i,j-1])
            neighborN = 0.0
         else:
            myCoef += -2.0 / dySqr
            neighborS = 1.0 / dySqr * phi[i,j-1]
            neighborN = 1.0 / dySqr * phi[i,j+1]

         tmp = (RHS[i,j] - neighborL - neighborR - neighborS - neighborN) / myCoef
         phiIncrement = alpha * (tmp - phi[i,j])
         residual += phiIncrement ** 2
         newPhi[i,j] = phi[i,j] + phiIncrement

   numberOfInnerPoints = imax * jmax
   residual = residual / numberOfInnerPoints
   residual = np.sqrt(residual)

   return newPhi, residual


def computeResidual(phi, RHS, imax, jmax):
   residual = 0.0
   dx = domainVars.dx
   dy = domainVars.dy
   dxSqr = dx * dx
   dySqr = dy * dy
   for j in range(jmax-1):
      if j == 0: continue
      for i in range(imax-1):
         if i == 0: continue
         # This is coefficient for my PHI at location (i,j).
         me        = -2.0 * (1.0/dxSqr + 1.0/dySqr) * phi[i,j]
         # coeff * T_{i-1,j} in finite difference equation
         neighborL = 1.0 / dxSqr * phi[i-1,j]
         # coeff * T_{i+1,j} in finite difference equation
         neighborR = 1.0 / dxSqr * phi[i+1,j]
         # coeff * T_{i,j-1} in finite difference equation
         neighborS = 1.0 / dySqr * phi[i,j-1]
         # coeff * T_{i,j+1} in finite difference equation
         neighborN = 1.0 / dySqr * phi[i,j+1]
         residual += ( RHS[i,j] - me - neighborL - neighborR - neighborS - neighborN ) ** 2

   numberOfInnerPoints = (imax-2) * (jmax-2)
   residual = residual / numberOfInnerPoints
   residual = np.sqrt(residual)
   return residual
