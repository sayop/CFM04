import numpy as np
from variables import *

def initSimulationVars(inputDict):
   print '# Initializing flow variables...'
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])
   dtInit = float(inputDict['dtInit'])

   flowVars.nu   = float(inputDict['nu'])
   flowVars.Uref = float(inputDict['uUp'])

   pInit = float(inputDict['pInit'])
   flowVars.p    = pInit * np.ones((imax,jmax))
   flowVars.u    = np.zeros((imax,jmax))
   flowVars.v    = np.zeros((imax,jmax))

   # initial time step
   timeVars.dtInit = dtInit

   # Evaluate reference Re number
   flowVars.Re = flowVars.Uref * domainVars.Lref / flowVars.nu
   print "# - Reference Reynolds number: Re = ", flowVars.Re

   # Populate boundary values that is specified by user in inputs.in
   updateVelocityBC(inputDict,imax,jmax)


def updateVelocityBC(inputDict,imax,jmax):
   # update boundary values of u and v in dimensionalized form
   # Left boundary
   flowVars.u[0,:] = float(inputDict['uLeft'])
   flowVars.v[0,:] = float(inputDict['vLeft'])
   # Right boundary
   flowVars.u[imax-1,:] = float(inputDict['uRight'])
   flowVars.v[imax-1,:] = float(inputDict['vRight'])
   # Bottom boundary
   flowVars.u[:,0] = float(inputDict['uBottom'])
   flowVars.v[:,0] = float(inputDict['vBottom'])
   # Upper boundary
   flowVars.u[:,jmax-1] = float(inputDict['uUp'])
   flowVars.v[:,jmax-1] = float(inputDict['vUp'])

