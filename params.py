import numpy as np
from numpy import pi,sqrt
import json
import sys

N  = 2**10
nu = 1e-5
Ltotal = 1.
Ttotal = 1.
L = Ltotal/10.
alpha = 0.1

dx = 1./float(N)
dt = dx*dx
visc = 4.*pi*pi*nu
sqdx = sqrt(dx)

NT_estimate = int(Ttotal/dt)

N_eval = 100
N_skip = NT_estimate//N_eval
NT = N_skip * N_eval
t_eval = np.arange(1,NT,N_skip)*dt