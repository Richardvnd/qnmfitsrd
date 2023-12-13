"""

This calculates the 'spatial mismatch' (not to be confused with mismatch) between the predicted 
spatial pattern and the reconstructed pattern. The purpose of this code is to look at the mismatch
for multiple black holes - particularly those with negative remnant spin - as these do not 
appear to be reconstructed as well by the spatial reconstruction code. 


"""


import numpy as np
import matplotlib.pyplot as plt
import qnmfitsrd as qnmfits
from multiprocessing import Pool
from spatial_reconstruction import *
from qnm_visualisation import qnm_viz
from matplotlib.animation import FuncAnimation

l_max = 6
n_max = 6
t0 = {0:40., 1:18.5, 2:12., 3:8., 4:5.5, 5:3., 6:1.5, 7:0.}[n_max]

sim = qnmfits.SXS(ID=305, zero_time=(2,2))

mapping = [(2,2,0,1)]

QNMs = [(lam,mu,n,p) for lam in np.arange(2, l_max+1)
                        for mu in np.arange(-lam, lam+1)
                           for n in np.arange(0, n_max+1)
                              for p in (-1, +1)]

best_fit_linear = qnmfits.mapping_multimode_ringdown_fit(sim.times, 
                                        sim.h, 
                                        modes=QNMs.copy(),
                                        Mf=sim.Mf,
                                        chif=sim.chif_mag,
                                        t0=t0,
                                        mapping_modes=mapping,
                                        spherical_modes=[(l,m) for l in np.arange(2, l_max+1)
                                                               for m in np.arange(-l,l+1)])

sm = spatial_mismatch(mapping[0], best_fit_linear, sim.chif_mag, l_max)

print(sm) 