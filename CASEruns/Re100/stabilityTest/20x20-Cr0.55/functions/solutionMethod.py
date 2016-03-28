from variables import flowVars, domainVars, FDM
import numpy as np

def nondimensionalize(inputDict, idomain, iflow):
   # idomain: nondimensionlize only domain variables
   if idomain == 1:
      Lref = domainVars.Lref
      domainVars.x = domainVars.x / Lref
      domainVars.y = domainVars.y / Lref
      domainVars.dx = domainVars.dx / Lref
      domainVars.dy = domainVars.dy / Lref

   # iflow: nondimensionalize only flow variables
   if iflow == 1:
      Uref = flowVars.Uref
      flowVars.u = flowVars.u / Uref
      flowVars.v = flowVars.v / Uref
      flowVars.p = flowVars.p / Uref ** 2

def dimensionalize(inputDict, idomain, iflow):
   # idomain: nondimensionlize only domain variables
   if idomain == 1:
      Lref = domainVars.Lref
      domainVars.x = domainVars.x * Lref
      domainVars.dx = domainVars.dx * Lref
      domainVars.y = domainVars.y * Lref
      domainVars.dy = domainVars.dy * Lref
   # iflow: nondimensionalize only flow variables
   if iflow == 1:
      Uref = flowVars.Uref
      flowVars.u = flowVars.u * Uref
      flowVars.v = flowVars.v * Uref
      flowVars.p = flowVars.p * Uref ** 2

def populateFluxVectors(inputDict):
   # Make sure that all the flux vectors' elements should be evaluated 
   # with non-dimensionalized domain and flow variables.
   beta = float(inputDict['Beta'])
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])

   # Populate E vector
   FDM.E[0] = flowVars.u / beta
   FDM.E[1] = flowVars.u ** 2 + flowVars.p
   FDM.E[2] = flowVars.u * flowVars.v
   
   # Populate F vector
   FDM.F[0] = flowVars.v / beta
   FDM.F[1] = flowVars.u * flowVars.v
   FDM.F[2] = flowVars.v ** 2 + flowVars.p

   # Populate D vector
   FDM.D[0] = np.zeros((imax,jmax))
   FDM.D[1] = flowVars.u
   FDM.D[2] = flowVars.v


def updateQvector(inputDict, dt):
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])

   for n in range(3):
      if n == 0: 
         # if True, the flux vector will have non-zero value at boundary
         # in order to update the primative variable at such boundary.
         #updateBoundary = True
         updateBoundary = False
      else:
         updateBoundary = False
      # clean Q vector to store new values
      FDM.Q[n] = np.zeros((imax,jmax))
      # Add convective flux in x-direction with E vector
      FDM.Q[n] += centralFiniteDifference(-FDM.E[n],'x',1,updateBoundary)
      # Add convective flux in y-direction with F vector
      FDM.Q[n] += centralFiniteDifference(-FDM.F[n],'y',1,updateBoundary)
      # Add diffusive flux in x- and y-directions with D vector
      # continuity equation doesn't incude second derivative, so skip!
      if n == 0: continue
      FDM.Q[n] += centralFiniteDifference(FDM.D[n],'x',2,updateBoundary) / flowVars.Re
      FDM.Q[n] += centralFiniteDifference(FDM.D[n],'y',2,updateBoundary) / flowVars.Re


def centralFiniteDifference(phi, direction, nOrder, updateBoundary):
   imax = len(phi[:,0])
   jmax = len(phi[0,:])
   dx = domainVars.dx
   dy = domainVars.dy
   f = np.zeros((imax,jmax))
   # first order spatial derivative in central difference
   if nOrder == 1:
      # x-derivative
      if direction == 'x':
         for j in range(jmax-1):
            if updateBoundary == False and j == 0: continue
            for i in range(imax-1):
               if updateBoundary == False and i == 0: continue
               # if current node is located at boundary, forward or backward difference is used. (First order)
               if i == 0:
                  # forward difference
                  f[i,j] = (phi[i+1,j] - phi[i,j]) / dx
               if i == imax-1:
                  # backward difference
                  f[i,j] = (phi[i,j] - phi[i-1,j]) / dx
               else:
                  f[i,j] = 0.5 * (phi[i+1,j] - phi[i-1,j]) / dx
      # y-derivative
      if direction == 'y':
         for i in range(imax-1):
            if updateBoundary == False and i == 0: continue
            for j in range(jmax-1):
               if updateBoundary == False and j == 0: continue
               # if current node is located at boundary, forward or backward difference is used. (First order)
               if j == 0:
                  # forward difference
                  f[i,j] = (phi[i,j+1] - phi[i,j]) / dy
               if j == jmax-1:
                  # backward difference
                  f[i,j] = (phi[i,j] - phi[i,j-1]) / dy
               else:
                  f[i,j] = 0.5 * (phi[i,j+1] - phi[i,j-1]) / dy
         
   # second order spatial derivative in central difference
   if nOrder == 2:
      # x-derivative
      if direction == 'x':
         for j in range(jmax-1):
            if updateBoundary == False and j == 0: continue
            for i in range(imax-1):
               if updateBoundary == False and i == 0: continue
               f[i,j] = (phi[i+1,j] - 2.0*phi[i,j] + phi[i-1,j]) / dx ** 2
      # y-derivative
      if direction == 'y':
         for i in range(imax-1):
            if updateBoundary == False and i == 0: continue
            for j in range(jmax-1):
               if updateBoundary == False and j == 0: continue
               f[i,j] = (phi[i,j+1] - 2.0*phi[i,j] + phi[i,j-1]) / dy ** 2

   return f

def updateTimeStep(inputDict):
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])
   Cr   = float(inputDict['Courant'])
   beta = float(inputDict['Beta'])
   # Artificial speed of sound
   a    = 1.0 / np.sqrt(beta)
   dx   = domainVars.dx
   dy   = domainVars.dy
   dt   = 999999.9
   for j in range(jmax-1):
      if j == 0: continue
      for i in range(imax-1):
         if i == 0: continue
         Uconv = abs(flowVars.u[i,j]) + a
         Vconv = abs(flowVars.v[i,j]) + a
         tmp   = (Uconv / dx + Vconv / dy)
         dt    = min(dt, Cr / tmp)

   return dt


def computeResidual(imax, jmax, dt, Q):
   res = np.zeros(3)
   # compute residual for pressure, u, and v
   for n in range(3):
      for j in range(jmax):
         for i in range(imax):
            res[n] += (dt * Q[n][i,j]) ** 2
   
      res[n] = np.sqrt( res[n] / (imax*jmax) )

   return res


def computeMaximumMach(imax, jmax, beta):
   Usqr = flowVars.u ** 2
   Vsqr = flowVars.v ** 2
   Asqr    = 1.0 / beta
   MachX = np.sqrt( Usqr / Asqr )
   MachX = MachX.max()
   MachY = np.sqrt( Vsqr / Asqr )
   MachY = MachY.max()
   return MachX, MachY


def updatePressureBC(imax, jmax):
   # left boundary
   i = 0
   for j in range(jmax-1):
      if j == 0: continue
      flowVars.p[i,j] = flowVars.p[i+1,j]

   # right boundary
   i = imax - 1
   for j in range(jmax-1):
      if j == 0: continue
      flowVars.p[i,j] = flowVars.p[i-1,j]

   # bottom boundary
   j = 0
   for i in range(imax-1):
      if i == 0: continue
      flowVars.p[i,j] = flowVars.p[i,j+1]

   # upper boundary
   j = jmax - 1
   for i in range(imax-1):
      if i == 0: continue
      flowVars.p[i,j] = flowVars.p[i,j-1]

   # update corner ponints
   flowVars.p[0,0] = flowVars.p[1,1]
   flowVars.p[imax-1,0] = flowVars.p[imax-2,1]
   flowVars.p[0,jmax-1] = flowVars.p[1,jmax-2]
   flowVars.p[imax-1,jmax-1] = flowVars.p[imax-2,jmax-2]
