import numpy as np
from numpy import pi,sqrt
import json
import sys

N  = 2**8
nu = 0.01
Ltotal = 1.
Ttotal = 1.
L = Ltotal/10.
alpha = 0.02

dx = 1./float(N)
dt = 0.05*dx*dx
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
