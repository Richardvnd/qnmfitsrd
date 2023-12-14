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
import datetime 

l_max = 4
n_max = 4
t_start = -50
t_end = 50
t_step = 1
times_list = np.arange(t_start, t_end, t_step)

mapping = [(2,2,0,1)]
map = mapping[0]

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

lon = np.linspace(-np.pi, np.pi, 200)
lat = np.linspace(-np.pi/2, np.pi/2, 200)
Lon, Lat = np.meshgrid(lon, lat)

bh_ids = [8, 13, 21, 34, 55, 89, 64, 114, 115, 207, 305, 1419, 1420, 1421, 1422, 1423, 1425, 1427]
#bh_ids = [305, 1427] 

sm_dict = {}
arg_dict = {}
spin_dict = {} 

for bh in bh_ids:

   try:
      sim = qnmfits.SXS(ID=bh, zero_time=(2,2))
   except Exception as e:
      print(f"Error occurred for bh={bh}: {str(e)}")
      continue

   G = spheroidal(np.pi/2-Lat, Lon, map, l_max, sim.chif_mag)

   sm_dict[bh] = []
   arg_dict[bh] = []

   for step in times_list:

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
      sm_dict[bh].append(sm)
      arg_dict[bh].append(arg)
      spin_dict[bh] = sim.chif[2]


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))

for bh in bh_ids:
   if spin_dict[bh] < 0:
      color = 'r'
      ax2.plot(times_list, sm_dict[bh], label=f'ID={bh}', color=color)
   else:
      color = 'b'
      ax1.plot(times_list, sm_dict[bh], label=f'ID={bh}', color=color)
   ax3.plot(times_list, arg_dict[bh], label=f'ID={bh}', color=color)

ax1.set_ylabel('Spatial Mismatch')
ax2.set_ylabel('Spatial Mismatch')
ax3.set_ylabel('Argument')
ax3.set_xlabel('Time [M]')

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=len(bh_ids)/2)
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=len(bh_ids)/2)
ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=len(bh_ids)/2)

plt.tight_layout()
plt.savefig('spatial_mismatch_spins.png')
plt.show()

