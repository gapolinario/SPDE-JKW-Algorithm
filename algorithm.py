import numpy as np
import sys
from numpy import pi,sqrt,exp
from numpy.fft import fft,ifft,fftfreq

from ext_params import *
from params import *

all_params = {"R": R, "N": N, "nu": nu, "Ltotal": Ltotal, "Ttotal": Ttotal,
              "L": L, "alpha": alpha, "dx": dx, "dt/dx^2": dt/dx**2, "NT": NT,
              "visc": visc, "sqdx": sqdx}
              #, "t_eval": t_eval, "N_eval": N_eval, "N_skip": N_skip}

with open(f'data/params_R_{R:06d}.json','w',encoding="utf-8") as file:
    json.dump(all_params, file)

#### FUNCTIONS ####

# vx: argument given in real space
def NonlinearDriftFunction(vx):

    # 1. linear
    #return -np.full((N,),alpha) * vx
    return - alpha * vx

    # 2.
    #return alpha*(1.-vx)/(1.+vx*vx)

# v0: argument given in Fourier space
# return is given in Fourier space too
def NonlinearDriftFunctionFourier(v0):

    # F(u) = - alpha * (1.-vx) / (1.+vx*vx)

    vx  = ifft(v0)
    tmp = ifft(Dealias(v0))

    return fft(-alpha*(1.-vx)/(1.+tmp))

def EulerMaruyamaStep(v0,f0):

    # this term holds Fourier[ F(X)], it is an array
    # F(X) is the nonlinear drift contribution
    F_fourier  = NonlinearDriftFunction( ifft(v0) )
    F_fourier  = fft(F_fourier)

    #F_fourier  = NonlinearDriftFunctionFourier( v0 )

    v0 -= visc*dt*K2*v0

    v0 += dt * F_fourier

    v0 += sqrt(dt) * f0

    return v0

def PredictorCorrectorStep(v0,f0):

    # this term holds Fourier[ F(X)], it is an array
    # F(X) is the nonlinear drift contribution
    F_fourier  = NonlinearDriftFunction( ifft(v0) )
    F_fourier  = fft(F_fourier)

    #F_fourier  = NonlinearDriftFunctionFourier( v0 )

    # predictor step
    vp  = np.copy(v0)
    vp -= visc*dt*K2*v0
    vp += dt * F_fourier
    vp += sqrt(dt) * f0

    # 1st half of corrector step
    vn  = np.copy(v0)
    vn += sqrt(dt) * f0
    vn -= visc*.5*dt*K2*v0
    vn += .5 * dt * F_fourier

    # 2nd half of corrector step
    F_fourier  = NonlinearDriftFunction( ifft(vp) )
    F_fourier  = fft(F_fourier)

    vn -= visc*.5*dt*K2*vp
    vn += .5 * dt * F_fourier

    return vn

def JentzenKloedenWinkelStep(v0,f0):

    # this term holds Fourier[ F(X)], it is an array
    # F(X) is the nonlinear drift contribution
    #F_fourier  = NonlinearDriftFunction( np.fft.ifft(v0) )
    #F_fourier  = fft(F_fourier)

    F_fourier = 0.

    #F_fourier  = NonlinearDriftFunctionFourier( v0 )

    # nonlinear drift
    v0 += dt * F_fourier
    # linear drift, exponentiated
    v0 *= exp(-visc*dt*K2)

    v0 += noise_weight * f0

    return v0

# Gets v in Fourier space
# Returns FFT(u^2) dealised where u = IFFT(v)
def Dealias(v0):

    M = N//2*3

    vpad  = np.zeros(M,dtype=np.complex128)

    vpad[0:N//2] = v0[0:N//2]
    vpad[M-N//2:M] = v0[N//2:]

    vpad  = fft( ifft(vpad)**2 )
    vpad *= M/N

    vf = np.zeros(N,dtype=np.complex128)
    vf[0:N//2] = vpad[0:N//2]
    vf[N//2:]    = vpad[M-N//2:M]

    return vf

# ### INTEGRATION ####

# arrays
v0 = np.zeros((N,),dtype=np.complex128)
f0 = np.zeros((N,),dtype=np.complex128)

X = fftfreq(N) * Ltotal
K = fftfreq(N) * N
K2 = K*K

kernel = exp(-.5*X**2/L/L) # exponential correlation function
kernel = sqrt(fft(kernel))

# noise
noise_weight  = np.ones(N)
noise_weight -= exp( -2.*visc*dt*K2 )
with np.errstate(divide='ignore', invalid='ignore'):
    noise_weight *= .5/visc/K2
noise_weight      = sqrt(noise_weight)
noise_weight[0]   = sqrt(dt)

# velocity fields are saved to these arrays
v_real = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
v_four = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
f_four = np.empty((N_eval,N),dtype=np.complex128) # forcing
f_real = np.empty((N_eval,N),dtype=np.complex128) # forcing

# initial condition

for ii,t in enumerate(t_eval):

    for _ in range(N_skip):

        # new time, new forcing
        f0 = np.random.normal(size=(N,)) * sqhdx + 1j * np.random.normal(size=(N,)) * sqhdx
        f0 = kernel * fft(f0)

        # CHOOSE algorithm
        # 1. JKW
        v0 = JentzenKloedenWinkelStep(v0,f0)
        # 2. Euler-Maruyama
        #v0 = EulerMaruyamaStep(v0,f0)
        # 3. Predictor-Corrector
        v0 = PredictorCorrectorStep(v0,f0)
        # END CHOOSE algorithm

    v_four[ii,:] = v0
    v_real[ii,:] = np.fft.ifft(v0)
    f_four[ii,:] = f0
    f_real[ii,:] = np.fft.ifft(f0)

fname = f'data/test_jkw_R_{R:06d}'
np.savez( fname , v_real=v_real, v_four=v_four, f_four=f_four, f_real=f_real )
