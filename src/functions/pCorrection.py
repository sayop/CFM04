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
   iterMax = 200 # maximum iteration for point-iterative method
   for niter in range(iterMax):
      if niter == 0: residualInit = max(1e-99,computeResidual(flowVars.p, RHS, imax, jmax))
      # run with point-iteration method: will update new pressure at new time level
      #flowVars.p, residual = pointIterJacobi(flowVars.p, RHS, imax, jmax, alpha)
      # Jacobi method
      #flowVars.p, res = pointIterJacobi(flowVars.p, RHS, imax, jmax, alpha)
      # Gauss-Seidel method
      flowVars.p, res = pointIterGS(flowVars.p, RHS, imax, jmax, alpha)
      #updatePressureBC(imax, jmax)

      residual = computeResidual(flowVars.p, RHS, imax, jmax)
      # update boundary condition for pressure only
      # compute residual at initial sate
      #if niter == 0: residualInit = max(1e-99,residual)
      residual = residual / residualInit
      #print niter, residual
      if residual <= convergeCrit or niter == iterMax-1:
         print '| Pressure Correction: converged at residual = %.6f' % float(residual), ' and iterations of %s' % niter
         break

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
         newPhi[i,j] = phi[i,j] + phiIncrement

   numberOfInnerPoints = (imax-2) * (jmax-2)
   residual = residual / numberOfInnerPoints
   residual = np.sqrt(residual)

   return newPhi, residual

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

   numberOfInnerPoints = (imax-2) * (jmax-2)
   residual = residual / numberOfInnerPoints
   residual = np.sqrt(residual)

   return phi, residual

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
