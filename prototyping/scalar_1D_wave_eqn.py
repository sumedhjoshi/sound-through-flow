#!/usr/local/bin/python

import numpy as np
import os
import sys
#from matplotlib import pyplot as plot
import math

# Set some inputs.
n  = 512
L  = 1000.0
x  = np.linspace(-L,L,n)
dx = x[1] - x[0]
w  = 2.0*L/25.0

# Set some physical constants.
rho0 = 1000.0
c    = np.array( [1500.0]*len(x) )
beta = rho0*c**2


# Time stepping parameters.
dt     = dx/max(c)/10.0
tfinal = 0.3
nt     = round( tfinal / dt )
t      = np.linspace( 0, tfinal, nt )

# Set the initial condition.
u0 = np.exp( - ( x**2.0 ) / w**2.0 )
u0 = 0.0 * u0
s0 = np.exp( - ( x**2.0 ) / w**2.0 )

# Build the finite difference matrices.
two = [-2.0]*len(x)
one = [1.0]*(len(x)-1)
D2 = np.diag( two ) + np.diag( one, -1) + np.diag( one, 1 )
D2 = D2/dx**2.0
D1 = np.diag( one, 1) - np.diag( one, -1 )
D1 = D1 / ( 2.0 * dx )

# Set up the matrices for the time-stepping.
U = np.zeros( ( len(x), len(t) ) )
S = np.zeros( ( len(x), len(t) ) )
U[:,0] = 0.0*u0
S[:,0] = s0
S2 = np.zeros( ( len(x), len(t) ) );
S2[:,0] = s0
S2[:,1] = s0

# Loop over time, explicitly computing the derivative.
for iit in range( 1,len(t) ):

   # Compute the 1st-order hyperbolic terms.
   S[:,iit] = S[:,iit-1] - dt * np.dot( D1, U[:,iit-1] )
   U[:,iit] = U[:,iit-1] - dt * np.dot( D1, S[:,iit-1] ) * c**2

#for iit in range( 2, len(t) ):
#
#   # Compute the 2nd-order time-difference solution.
#   S2[:,iit] = 2.0*S2[:,iit-1] - S2[:,iit-2] + ( beta / rho0 ) * dt**2 * np.dot( D2, S2[:,iit-2] )

# Set the data.
Sdata = S
Udata = U

#
#
# Apply the adjoint method to invert for beta/rho0.

# Set some constants.
Nadj = 50

# Make a bad initial guess.
B0 = beta * ( 1. )
R0 = rho0 * ( 1. )
C = np.array( np.sqrt(B0/R0) - 0.01 * x )
Ck     = np.zeros( (len(x), Nadj) )
Gk     = np.zeros( (Nadj) )
gradGk = np.zeros( (len(x), Nadj) )
Ck[:,0] = C
#C = np.array( [np.sqrt(B0/R0)] * len( x ) )

print " Initial guess for sound speed: " + str(C[0])

# Apply the adjoint method.
for ii in range(1,Nadj):

   print " Adjoint method iteration number " + str(ii) + "\n"

   # Step 1: Compute the field based on the current guess.
   iiU = np.zeros( (len(x),len(t)) )
   iiS = np.zeros( (len(x),len(t)) )
   iiU[:,0] = u0
   iiS[:,0] = s0
   for iit in range( 1, len(t) ):

      # Compute the first order hyperbolic equations.
      iiS[:,iit] = iiS[:,iit-1] - dt * np.dot( D1, iiU[:,iit-1] )
      iiU[:,iit] = iiU[:,iit-1] - dt * np.dot( D1, iiS[:,iit-1] ) * C * C

   # Step 2: Compute the right-hand-side for the adjoint problem.
   diffU = 2 * ( iiU - Udata )
   diffS = 2 * ( iiS - Sdata )

   # Step 3: Solve the adjoint problem.
   iiUa  = np.zeros( (len(x),len(t)) )
   iiSa  = np.zeros( (len(x),len(t)) )
   #iiUa[:,0] = -1.0 * diffU[:,-1]
   #iiSa[:,0] = -1.0 * diffS[:,-1]
   iiUa[:,0] = diffU[:,-1]
   iiSa[:,0] = diffS[:,-1]
   for iit in range( 1, len(t) ):

      # Compute the first order hyperbolic equations.
      #iiSa[:,iit] = iiSa[:,iit-1] - dt * np.dot( np.transpose(D1), iiUa[:,iit-1] )
      #iiUa[:,iit] = iiUa[:,iit-1] - dt * np.dot( np.transpose(D1), iiSa[:,iit-1] ) * C * C
      iiSa[:,iit] = iiSa[:,iit-1] - dt * np.dot( np.transpose(D1), iiUa[:,iit-1] ) * C * C
      iiUa[:,iit] = iiUa[:,iit-1] - dt * np.dot( np.transpose(D1), iiSa[:,iit-1] )

   # Step 3.5: Build the derivative of the operator wrt to the parameters.
   DLu = np.zeros( (len(x), len(t)) )
   for iit in range(0,len(t)):
      DLu[:,iit] = np.dot( D1, iiS[:,iit] )

   # Step 4: Compute the gradient of the data-model mismatch functional.
   #gradG = ( - 2.0 * C ) * np.trapz( np.trapz( DLu * iiUa)) * dx * dt
   gradG = ( 2.0 * C ) * np.trapz( np.trapz( DLu * iiUa)) * dx * dt

   # Step 5: Update the guess for the parameters.
   G = np.trapz(np.trapz((iiU - Udata)**2))*dx*dt
   alpha = 1.
   step_size = 25.0
   C = C - gradG*step_size

   # Step 6: Store some data.
   Ck[:,ii]     = C
   gradGk[:,ii] = gradG
   Gk[ii]       = G

   # Step 7: Print a status report.
   print "     This estimate of the gradient: " + str(gradG[0])
   print "     This guess for sound speed: "    + str(C[0])
   print "     The current residual value: "    + str(np.sqrt(G)) + "\n"
