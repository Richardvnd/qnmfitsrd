"""

Updated: 13/12/2023

This code uses the qnm_viz class to generate animations of the spatial reconstruction. 

"""

import numpy as np
import matplotlib.pyplot as plt
import qnmfitsrd as qnmfits
from multiprocessing import Pool
from spatial_reconstruction import *
from qnm_visualisation import qnm_viz
from matplotlib.animation import FuncAnimation
import datetime 

sim = qnmfits.SXS(ID=305, zero_time=(2,2))

l_max = 2
n_max = 2
t0 = {0:40., 1:18.5, 2:12., 3:8., 4:5.5, 5:3., 6:1.5, 7:0.}[n_max]
qnm_viz = qnm_viz(sim, l_max=l_max)

mapping = [(2,2,0,1)]

if len(mapping[0])==8:
   QNMs = [(lam,mu,n,p) for lam in np.arange(2, l_max+1)
                            for mu in np.arange(-lam, lam+1)
                               for n in np.arange(0, n_max+1)
                                  for p in (-1, +1)] + mapping
else:
   QNMs = [(lam,mu,n,p) for lam in np.arange(2, l_max+1)
                           for mu in np.arange(-lam, lam+1)
                              for n in np.arange(0, n_max+1)
                                 for p in (-1, +1)]



fig, axs = plt.subplots(nrows=2, ncols=2, 
                     subplot_kw={'projection': 'mollweide'}, 
                     figsize=(12,8))

Lon, Lat = qnm_viz.latlon() 

map = mapping[0]

G = spheroidal(np.pi/2-Lat, Lon, map, l_max, sim.chif_mag)

id = datetime.datetime.now()

def update(step):
   best_fit = qnmfits.mapping_multimode_ringdown_fit(sim.times, 
                                             sim.h, 
                                             modes=QNMs.copy(),
                                             Mf=sim.Mf,
                                             chif=sim.chif_mag,
                                             t0=step,
                                             mapping_modes=mapping,
                                             spherical_modes=[(l,m) for l in np.arange(2, l_max+1)
                                                               for m in np.arange(-l,l+1)])

   F = mode_mapping(np.pi/2-Lat, Lon, best_fit, map, l_max)

   fig.suptitle(f"t_0 = {step}", fontsize=16)

   axs[0,0].title.set_text(str(map)+'_real')
   axs[0,0].pcolormesh(Lon, Lat, np.real(F), cmap=plt.cm.jet)

   axs[1,0].title.set_text(str(map)+'_imag')
   axs[1,0].pcolormesh(Lon, Lat, np.imag(F), cmap=plt.cm.jet)

   axs[0,1].title.set_text('S_real')
   axs[0,1].pcolormesh(Lon, Lat, np.real(G), cmap=plt.cm.jet)

   axs[1,1].title.set_text('S_imag')
   axs[1,1].pcolormesh(Lon, Lat, np.imag(G), cmap=plt.cm.jet)

   return fig 

#ani = FuncAnimation(fig, update, frames=range(-10, 10, 1), interval=50)
#ani.save(f'mapping_animation_{map}_{id}.mp4', writer='ffmpeg')


#ani = qnm_viz.animate_mapping_projection(sim, QNMs, mapping, min_t0 = -100, max_t0 = 200, step = 1, save = True)

# Calculate how fast the lobes are moving

T = 25
print("Velocity of lobes ~ {:.2f}c".format(4*np.pi*sim.Mf/T))