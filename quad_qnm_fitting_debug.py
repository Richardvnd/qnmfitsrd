import matplotlib.pyplot as plt
import numpy as np

import qnmfitsrd as qnmfits

sim = qnmfits.SXS(ID=305, zero_time=(2,2))

modes = [(5,5,0,1,5,5,0,1), 
              (6,6,0,1,6,6,0,1),
              (7,7,0,1,7,7,0,1),
              (8,8,0,1,8,8,0,1),
              (9,9,0,1,9,9,0,1)]

best_fit = qnmfits.multi_multimode_ringdown_fit(
    sim.times,
    sim.h,
    modes,
    Mf=sim.Mf,
    chif=sim.chif_mag,
    t0=0,
)

print(f"Mismatch = {best_fit['mismatch']}")



