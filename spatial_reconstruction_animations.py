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

l_max = 4
n_max = 4
t_start = 0
t_end = 50
t_step = 1
t0 = {0:40., 1:18.5, 2:12., 3:8., 4:5.5, 5:3., 6:1.5, 7:0.}[n_max]
id = datetime.datetime.now()

qnm_viz = qnm_viz(sim, l_max=l_max)

mapping = [(4,4,0,1)]

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

fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(12,8))

axs[0, 0].axis('off')  
axs[0, 1].axis('off')  
axs[1, 0].axis('off')  
axs[1, 1].axis('off')  
axs[3, 0].axis('off')  
axs[0, 0] = plt.subplot(4, 2, 1, projection='mollweide')
axs[0, 1] = plt.subplot(4, 2, 2, projection='mollweide')
axs[1, 0] = plt.subplot(4, 2, 3, projection='mollweide')
axs[1, 1] = plt.subplot(4, 2, 4, projection='mollweide')


Lon, Lat = qnm_viz.latlon() 

map = mapping[0]

G = spheroidal(np.pi/2-Lat, Lon, map, l_max, sim.chif_mag)

sm_list = []
arg_list = []
z_list = [] 

times_list = [] 

def update(step):

   times_list.append(step)

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
   sm, arg, z = spatial_mismatch(F, G, num_points=100)
   sm_list.append(sm)
   arg_list.append(arg)
   z_list.append(z)

   fig.suptitle(f"Map: {map} \n $t_0$ = {step}", fontsize=16)

   axs[0, 0].title.set_text('Real reconstruction')
   axs[0, 0].title.set_position([0.5, 1.05])  
   axs[0, 0].pcolormesh(Lon, Lat, np.real(F), cmap=plt.cm.jet)

   axs[1,0].title.set_text('Imag reconstruction')
   axs[1,0].title.set_position([0.5, 1.05])  
   axs[1,0].pcolormesh(Lon, Lat, np.imag(F), cmap=plt.cm.jet)

   axs[0, 1].title.set_text('Real calculated')
   axs[0, 1].pcolormesh(Lon, Lat, np.real(G), cmap=plt.cm.jet)

   axs[1,1].title.set_text('Imag calculated')
   axs[1,1].pcolormesh(Lon, Lat, np.imag(G), cmap=plt.cm.jet)

   axs[2, 0].plot(times_list, sm_list, color='black')
   axs[2, 1].plot(times_list, arg_list, color='black')

   axs[2, 0].set_ylabel('Spatial Mismatch')
   axs[2, 1].set_xlabel('Time')
   axs[2, 0].set_xlabel('Time')
   axs[2, 1].set_ylabel('Argument of Z') 

   axs[3, 1].plot(np.real(z_list), np.imag(z_list), color='black')
   axs[3, 1].set_xlabel('Re(Z)')
   axs[3, 1].set_ylabel('Im(Z)') 

   return fig 

ani = FuncAnimation(fig, update, frames=range(t_start, t_end, t_step), interval=100)
ani.save(f'mapping_animation_{map}_{id}.mp4', writer='ffmpeg')

diff_arg = np.diff(arg_list)
# Look for the discontinuity in the argument of Z
jump_indices = np.where(abs(diff_arg) > 2) 
period = np.mean(np.diff(jump_indices)) * t_step
print(period) 

#ani = qnm_viz.animate_mapping_projection(sim, QNMs, mapping, min_t0 = -100, max_t0 = 200, step = 1, save = True)