import numpy as np
import sys
from numpy import pi,sqrt,exp

from ext_params import *
from params import *

all_params = {"R": R, "N": N, "nu": nu, "Ltotal": Ltotal, "Ttotal": Ttotal,
              "L": L, "alpha": alpha, "dx": dx, "dt": dt, "NT": NT,
              "visc": visc, "sqdx": sqdx}
              #, "t_eval": t_eval, "N_eval": N_eval, "N_skip": N_skip}

with open(f'data/params_R_{R:06d}.json','w',encoding="utf-8") as file:
    json.dump(all_params, file)




#### FUNCTIONS ####

# vx: argument given in real space
def NonlinearDriftFunction(vx):

    return np.full((N,),alpha) * vx

def JentzenKloedenWinkelStep(v0,f0):

    # this term holds Fourier[ F(X)], it is an array
    # F(X) is the nonlinear drift contribution
    F_fourier = NonlinearDriftFunction(np.fft.ifft(v0))
    F_fourier = np.fft.fft(F_fourier)

    # zero mode, linear part
    v0[0] += sqrt(dt) * f0[0]
    # zero mode, nonlinear drift
    v0[0] += dt * F_fourier[0]

    # other modes
    cte = dt*visc;

    # nonlinear drift
    v0[1:] += dt * F_fourier[1:]
  	# linear drift
    v0[1:] *= exp(-cte*K[1:]**2)
    # noise
    weight  = 1.-exp(-2.*visc*dt*K[1:]**2)
    weight *= .5/visc/K[1:]**2
    weight  = sqrt(weight)

    v0[1:] += weight * f0[1:]

    return v0





#### INTEGRATION ####

# arrays
v0 = np.zeros((N,),dtype=np.complex128)
f0 = np.zeros((N,),dtype=np.complex128)

X = np.fft.fftfreq(N) * Ltotal
K = np.fft.fftfreq(N) * N
K2 = K*K

kernel = exp(-.5*X**2/L/L) # exponential correlation function
kernel = sqrt(np.fft.fft(kernel))

# velocity fields are saved to these arrays
v_real = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
v_four = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
f_four = np.empty((N_eval,N),dtype=np.complex128) # forcing

for ii,t in enumerate(t_eval):

    # new time, new forcing
    f0 = np.random.normal(size=(N,)) * sqdx
    f0 = kernel * np.fft.fft(f0)

    for i in range(N_skip):
        v0 = JentzenKloedenWinkelStep(v0,f0)

    v_four[ii,:] = v0
    v_real[ii,:] = np.fft.ifft(v0)
    f_four[ii,:] = f0

fname = f'data/test_jkw_R_{R:06d}'
np.savez( fname , v_real=v_real, v_four=v_four, f_four=f_four )
