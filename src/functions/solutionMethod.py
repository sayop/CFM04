from variables import flowVars, domainVars, timeVars, FDM

import numpy as np

def nondimensionalize(inputDict, idomain, iflow, itime):
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

   # itime: nondimensionalize only time variables
   if itime == 1:
      Lref = domainVars.Lref
      Uref = flowVars.Uref
      timeVars.t = timeVars.t * Uref / Lref
      timeVars.dt = timeVars.dt * Uref / Lref
      timeVars.dtInit = timeVars.dtInit * Uref / Lref
      

def dimensionalize(inputDict, idomain, iflow, itime):
   # idomain: dimensionlize only domain variables
   if idomain == 1:
      Lref = domainVars.Lref
      domainVars.x = domainVars.x * Lref
      domainVars.dx = domainVars.dx * Lref
      domainVars.y = domainVars.y * Lref
      domainVars.dy = domainVars.dy * Lref
   # iflow: dimensionalize only flow variables
   if iflow == 1:
      Uref = flowVars.Uref
      flowVars.u = flowVars.u * Uref
      flowVars.v = flowVars.v * Uref
      flowVars.p = flowVars.p * Uref ** 2

   # itime: dimensionalize only time variables
   if itime == 1:
      Lref = domainVars.Lref
      Uref = flowVars.Uref
      timeVars.t = timeVars.t * Lref / Uref
      timeVars.dt = timeVars.dt * Lref / Uref
      timeVars.dtInit = timeVars.dtInit * Lref / Uref

def populateFluxVectors(pCorr,inputDict):
   # Make sure that all the flux vectors' elements should be evaluated 
   # with non-dimensionalized domain and flow variables.
   beta = float(inputDict['Beta'])
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])

   # For artificial compressibility method (ACM)
   if pCorr != 1:
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

   if pCorr == 1:
      # Populate E vector
      FDM.E[0] = np.zeros((imax,jmax))
      FDM.E[1] = flowVars.u ** 2
      FDM.E[2] = flowVars.u * flowVars.v

      # Populate F vector
      FDM.F[0] = np.zeros((imax,jmax))
      FDM.F[1] = flowVars.u * flowVars.v
      FDM.F[2] = flowVars.v ** 2

      # Populate D vector
      FDM.D[0] = np.zeros((imax,jmax))
      FDM.D[1] = flowVars.u
      FDM.D[2] = flowVars.v


def updateQvector(inputDict, dt, nSub, pCorr):
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])

   for n in range(3):
      
      if n == 0:
         # if True, the flux vector will have non-zero value at boundary
         # in order to update the primative variable at such boundary.
         # At boundary, no BC needed for continuity equation.
         updateBoundary = True
         #updateBoundary = False
      else:
         updateBoundary = False
      # clean Q vector to store new values
      FDM.Q[n] = np.zeros((imax,jmax))
      if nSub == 0:
         # pCorr == 1: in first step of updating intermediate time level of velocities,
         # u*, v* at boundaries are evaulated with FDM.
         if pCorr == 1: updateBoundary = True
         # Add convective flux in x-direction with E vector
         FDM.Q[n] += centralFiniteDifference(-FDM.E[n],'x',1,updateBoundary)
         # Add convective flux in y-direction with F vector
         FDM.Q[n] += centralFiniteDifference(-FDM.F[n],'y',1,updateBoundary)
         # Add diffusive flux in x- and y-directions with D vector
         # continuity equation doesn't incude second derivative, so skip!
         if n == 0: continue
         FDM.Q[n] += centralFiniteDifference(FDM.D[n],'x',2,updateBoundary) / flowVars.Re
         FDM.Q[n] += centralFiniteDifference(FDM.D[n],'y',2,updateBoundary) / flowVars.Re
      elif nSub == 1:
         # this is for pressure correction method
         # corrected pressure has been resolved, and future time level of primitive variables need to be updated here.
         if n == 0: continue
         if n == 1: FDM.Q[n] = centralFiniteDifference(-flowVars.p,'x',1,updateBoundary)
         if n == 2: FDM.Q[n] = centralFiniteDifference(-flowVars.p,'y',1,updateBoundary)
         #print 'project u* onto divergence free space'


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
         for j in range(jmax):
            if updateBoundary == False and j == 0: continue
            if updateBoundary == False and j == jmax-1: continue
            for i in range(imax):
               if updateBoundary == False and i == 0: continue
               if updateBoundary == False and i == imax-1: continue
               # if current node is located at boundary, forward or backward difference is used. (First order)
               if i == 0:
                  # forward difference
                  f[i,j] = (phi[i+1,j] - phi[i,j]) / dx
               elif i == imax-1:
                  # backward difference
                  f[i,j] = (phi[i,j] - phi[i-1,j]) / dx
               else:
                  f[i,j] = 0.5 * (phi[i+1,j] - phi[i-1,j]) / dx
      # y-derivative
      if direction == 'y':
         for i in range(imax):
            if updateBoundary == False and i == 0: continue
            if updateBoundary == False and i == imax-1: continue
            for j in range(jmax):
               if updateBoundary == False and j == 0: continue
               if updateBoundary == False and j == jmax-1: continue
               # if current node is located at boundary, forward or backward difference is used. (First order)
               if j == 0:
                  # forward difference
                  f[i,j] = (phi[i,j+1] - phi[i,j]) / dy
               elif j == jmax-1:
                  # backward difference
                  f[i,j] = (phi[i,j] - phi[i,j-1]) / dy
               else:
                  f[i,j] = 0.5 * (phi[i,j+1] - phi[i,j-1]) / dy
         
   # second order spatial derivative in central difference
   if nOrder == 2:
      # x-derivative
      if direction == 'x':
         for j in range(jmax):
            if updateBoundary == False and j == 0: continue
            if updateBoundary == False and j == jmax-1: continue
            for i in range(imax):
               if updateBoundary == False and i == 0: continue
               if updateBoundary == False and i == imax-1: continue
               # if current node is located at boundary, forward or backward difference is used. (Second order)
               if i == 0:
                  # forward difference
                  f[i,j] = (phi[i,j] - 2.0*phi[i+1,j] + phi[i+2,j]) / dx ** 2
               elif i == imax-1:
                  # backward difference
                  f[i,j] = (phi[i-2,j] - 2.0*phi[i-1,j] + phi[i,j]) / dx ** 2
               else:
                  f[i,j] = (phi[i+1,j] - 2.0*phi[i,j] + phi[i-1,j]) / dx ** 2
      # y-derivative
      if direction == 'y':
         for i in range(imax):
            if updateBoundary == False and i == 0: continue
            if updateBoundary == False and i == imax-1: continue
            for j in range(jmax):
               if updateBoundary == False and j == 0: continue
               if updateBoundary == False and j == jmax-1: continue
               # if current node is located at boundary, forward or backward difference is used. (Second order)
               if j == 0:
                  # forward difference
                  f[i,j] = (phi[i,j] - 2.0*phi[i,j+1] + phi[i,j+2]) / dy ** 2
               elif j == jmax-1:
                  # backward difference
                  f[i,j] = (phi[i,j-2] - 2.0*phi[i,j-1] + phi[i,j]) / dy ** 2
               else:
                  f[i,j] = (phi[i,j+1] - 2.0*phi[i,j] + phi[i,j-1]) / dy ** 2

   return f

def updateTimeStep(inputDict,nIter):
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])
   Cr   = float(inputDict['Courant'])
   pCorr = float(inputDict['pCorr'])
  
   if pCorr != 1:
      # Artificial speed of sound
      beta = float(inputDict['Beta'])
      a    = 1.0 / np.sqrt(beta)
   else:
      a = 0.0

   dx   = domainVars.dx
   dy   = domainVars.dy
   dt   = 0.25 * flowVars.Re * min(domainVars.dx, domainVars.dy) ** 2
   #if nIter == 1:
   #   dt = timeVars.dtInit
   for j in range(jmax-1):
      if j == 0: continue
      for i in range(imax-1):
         if i == 0: continue
         Uconv = max(1e-99,abs(flowVars.u[i,j])) + a
         Vconv = max(1e-99,abs(flowVars.v[i,j])) + a
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

def updatePrimitiveVars(pCorr,imax,jmax,dt,updateBoundary):

   # update primative vector U in interior points
   for j in range(jmax):
      if updateBoundary == False and j == 0: continue
      if updateBoundary == False and j == jmax-1: continue
      for i in range(imax):
         if updateBoundary == False and i == 0: continue
         if updateBoundary == False and i == imax-1: continue
         if pCorr != 1: flowVars.p[i,j] += dt * FDM.Q[0][i,j]
         flowVars.u[i,j] += dt * FDM.Q[1][i,j]
         flowVars.v[i,j] += dt * FDM.Q[2][i,j]



