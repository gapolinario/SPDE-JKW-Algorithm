import numpy as np
from numpy import pi,sqrt,exp
from scipy.special import erfc
from glob import glob
import matplotlib.pyplot as plt
import os

from params import *

# parent directory
path = os.path.abspath(os.path.join(os.path.abspath('')))

all_files = glob(f'{path}/data/*_R_{R:03d}*_N_{BN:02d}_dtconst_{dtconst:.06f}_nu_{nu:.06f}.npz')
all_files.sort()

print(f"R_{R:03d}_N_{BN:02d}_dtconst_{dtconst:.06f}_nu_{nu:.06f}")
print(f"Ensemble size is {len(all_files)}")

# variance of individual Fourier modes

spec = np.zeros((N_eval,N))

#taxis = np.arange(N_eval) * dt * N_skip

for f in all_files:

    npzfile = np.load(f)

    v_four = npzfile['v_four']

    spec += np.abs(v_four)**2 / N

spec *= 1./len(all_files)

fname = f'results/test_jkw_spectrum_R_{R:03d}_N_{BN:02d}_dtconst_{dtconst:.06f}_nu_{nu:.06f}'
np.save( fname , spec )
