import numpy as np
import sys
from numpy import pi,sqrt,exp
from numpy.fft import fft,ifft,fftfreq

from ext_params import *
from params import *

all_params = {"R": R, "N": N, "nu": nu, "Ltotal": Ltotal, "Ttotal": Ttotal,
              "L": L, "alpha": alpha, "dx": dx, "dt": dt, "NT": NT,
              "visc": visc, "sqdx": sqdx}
              #, "t_eval": t_eval, "N_eval": N_eval, "N_skip": N_skip}

with open(f'data/params_R_{R:06d}.json','w',encoding="utf-8") as file:
    json.dump(all_params, file)

# ### FUNCTIONS ####

# vx: argument given in real space
def NonlinearDriftFunction(vx):

    # 1. linear
    return - alpha * vx

    # 2.
    #return alpha*(1.-vx)/(1.+vx*vx)

def EulerMaruyamaStep(v0,f0):

    # this term holds Fourier[ F(X)], it is an array
    # F(X) is the nonlinear drift contribution
    F_fourier  = NonlinearDriftFunction( v0 )

    v0 -= visc*dt*v0

    v0 += dt * F_fourier

    v0 += sqrt(dt) * f0

    return v0

def JentzenKloedenWinkelStep(v0,f0):

    F_fourier  = NonlinearDriftFunction( v0 )

    # nonlinear drift
    v0 += dt * F_fourier
    # linear drift, exponentiated
    v0 *= exp(-visc*dt)

    v0 += noise_weight * f0

    return v0

# ### INTEGRATION ####

# arrays
v0 = np.zeros((N,),dtype=np.float64)
f0 = np.zeros((N,),dtype=np.float64)

X = fftfreq(N) * Ltotal
K = fftfreq(N) * N
K2 = K*K

kernel = exp(-.5*X**2/L/L) # exponential correlation function
kernel = sqrt(fft(kernel))

# noise
noise_weight  = np.ones(N)
noise_weight -= exp( -2.*visc*dt )
with np.errstate(divide='ignore', invalid='ignore'):
    noise_weight *= .5/visc
noise_weight      = sqrt(noise_weight)
noise_weight[0]   = sqrt(dt)

# velocity fields are saved to these arrays
v_all = np.empty((N_eval,N),dtype=np.float64) # stationary evolution, spatial profile
f_all = np.empty((N_eval,N),dtype=np.float64) # forcing

for ii,t in enumerate(t_eval):

    for _ in range(N_skip):

        # new time, new forcing
        f0 = np.random.normal(size=(N,))

        # CHOOSE algorithm
        # 1. JKW
        #v0 = JentzenKloedenWinkelStep(v0,f0)
        # 2. Euler-Maruyama
        v0 = EulerMaruyamaStep(v0,f0)
        # END CHOOSE algorithm

    v_all[ii,:] = v0
    f_all[ii,:] = f0

fname = f'data/test_sde_R_{R:06d}'
np.savez( fname , v_all=v_all, f_all=f_all )
