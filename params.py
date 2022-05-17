import numpy as np
from numpy import pi,sqrt
import json
import sys

R  = int(sys.argv[1])
BN = int(sys.argv[2])
dtconst = float(sys.argv[3])
nu = float(sys.argv[4])

N  = 2**BN
nu = 0.01
Ltotal = 1.
Ttotal = 1.
L = Ltotal/10.
alpha = 0.02

dx = 1./float(N)
dt = dtconst * dx*dx
# linear analysis suggests
# dt < dx^2 / ( pi^2 nu Ltotal^2 )
visc = 4.*pi*pi*nu
sqdx = sqrt(dx)
sqhdx = sqrt(.5*dx)

NT_estimate = int(Ttotal/dt)

N_eval = 10000
N_skip = NT_estimate//N_eval
NT = N_skip * N_eval
t_eval = np.arange(1,NT,N_skip)*dt
