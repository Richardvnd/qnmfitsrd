"""

Updated: 13/12/2023

This animates the spherical distribution of the spatial reconstruction.

"""


from matplotlib.animation import FuncAnimation
import numpy as np
import spherical
import matplotlib.pyplot as plt
import quaternionic
import numpy as np
from scipy.special import sph_harm

import qnmfitsrd as qnmfits
from spatial_reconstruction import mode_mapping

l_max = 9

sim = qnmfits.SXS(ID=305, zero_time=(2,2))

modes = [(2,2,n,1) for n in range(7+1)]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update(t): 

    t = t/10

    best_fit = qnmfits.multimode_ringdown_fit(
        sim.times,
        sim.h,
        modes,
        Mf=sim.Mf,
        chif=sim.chif_mag,
        t0=t
    )

    mapping = (2,2,0,1)

    lon = np.linspace(-np.pi, np.pi, 200)
    lat = np.linspace(-np.pi/2, np.pi/2, 200)

    Lon, Lat = np.meshgrid(lon, lat)
    F = mode_mapping(np.pi/2-Lat, Lon, best_fit, mapping, l_max)

    wigner = spherical.Wigner(l_max)

    def sYlm(theta, phi, l, m, s=-2):
        Y = np.empty_like(theta, dtype=complex)
        for i in range(theta.shape[0]):
            for j in range(theta.shape[1]):
                R = quaternionic.array.from_spherical_coordinates(theta[i, j], phi[i, j])
                Y[i, j] = wigner.sYlm(s, R)[wigner.Yindex(l, m)]
        return Y


    def Ylm(theta, phi, l, m):
        Y = sph_harm(m, l, phi, theta)
        return Y

    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)
    theta, phi = np.meshgrid(theta, phi) 

    l = 2
    m = 0

    spherical_harmonics = sYlm(theta, phi, l, m)

    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    ax.plot_surface(x, y, z, facecolors=plt.cm.jet(F.real), rstride=1, cstride=1)

    return ax, 

min_t0 = 0
max_t0 = 10000
step_t0 = 1

ani = FuncAnimation(fig, update, frames=range(min_t0, max_t0, step_t0), interval=0.0001)

plt.show()
