import numpy as np
import sys
from numpy import pi,sqrt,exp
from numpy.fft import fft,ifft,fftfreq

from ext_params import *
from params import *

# ### FUNCTIONS ####

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
    vf[N//2:]  = vpad[M-N//2:M]

    return vf


# arrays
v0 = np.zeros((N,),dtype=np.complex128)

X = fftfreq(N) * Ltotal
K = fftfreq(N) * N
K2 = K*K
K2 = 1.

kernel = exp(-.5*X**2/L/L) # exponential correlation function
kernel = sqrt(fft(kernel))

# velocity fields are saved to these arrays
v_real  = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
v_four  = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
v2_real = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile
v2_four = np.empty((N_eval,N),dtype=np.complex128) # stationary evolution, spatial profile

for ii in range(N_eval):

    # new time, new forcing
    # complex field
    #v0 = np.random.normal(size=(N,)) * sqhdx + 1j * np.random.normal(size=(N,)) * sqhdx
    # real field
    v0 = np.random.normal(size=(N,)) * sqdx
    v0 = kernel * fft(v0)

    # save u=IFFT(v) (in real and fourier space)
    v_four[ii,:] = v0
    v_real[ii,:] = ifft(v0)

    v0 = Dealias(v0)

    # save u^2, u=IFFT(v) (in real and fourier space)
    v2_four[ii,:] = v0
    v2_real[ii,:] = ifft(v0)

fname = f'data/test_dealias_R_{R:06d}'
np.savez( fname , v_real=v_real, v2_real=v2_real, v_four=v_four, v2_four=v2_four )
