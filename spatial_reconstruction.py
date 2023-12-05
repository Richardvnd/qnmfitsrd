import numpy as np
import pickle
import json
import matplotlib.pyplot as plt
import qnmfitsrd as qnmfits
import spherical
import quaternionic

def mode_mapping(theta, phi, best_fit, mapping, l_max):
    wigner = spherical.Wigner(l_max)
    R = quaternionic.array.from_spherical_coordinates(theta, phi)
    Y = wigner.sYlm(-2, R)

    ans = np.zeros_like(theta, dtype=complex)
    i = 0
    for loop in range(len(best_fit['C'])):
        A = best_fit['C'][loop]
        if best_fit['modes'][loop]==mapping:
            ans += A * Y[:,:,wigner.Yindex(*best_fit['spherical_modes'][i])]
            i += 1

    ans /= np.max(np.abs(ans)) # normalise peak value
    return ans

def spheroidal(theta, phi, mapping, l_max, chif):
    wigner = spherical.Wigner(l_max)
    R = quaternionic.array.from_spherical_coordinates(theta, phi)
    Y = wigner.sYlm(-2, R)

    ans = np.zeros_like(theta, dtype=complex)
    l, m, n, p = mapping
    for lp in np.arange(2, l_max+1):
        ans += qnmfits.qnm.mu(lp, m, l, m, n, p, chif) * Y[:,:,wigner.Yindex(lp, m)]

    ans /= np.max(np.abs(ans)) # normalise peak value
    return ans


def my_function(args):

    mapping, sim, t0_vals, l_max, n_max, num_points = args

    QNMs = [(lam,mu,n,p) for lam in np.arange(2, l_max+1)
                            for mu in np.arange(-lam, lam+1)
                               for n in np.arange(0, n_max+1)
                                  for p in (-1, +1)]

    spatial_mismatch = np.zeros_like(t0_vals, dtype=float)
    for i, t0 in enumerate(t0_vals):
        best_fit = qnmfits.mapping_multimode_ringdown_fit(sim.times,
                                                sim.h,
                                                modes=QNMs.copy(),
                                                Mf=sim.Mf,
                                                chif=sim.chif_mag,
                                                t0=t0,
                                                mapping_mode=mapping,
                                                spherical_modes=[(l,m) for l in np.arange(2, l_max+1)
                                                                       for m in np.arange(-l, l+1)])

        # I'm doing a simple, 2D "rectangle rule" integration
        # should probably switch to using dbl quadrature
        dx, dphi = 2./num_points, 2*np.pi/num_points
        x = np.arange(-1, 1, dx)
        phi = np.arange(-np.pi, np.pi, dphi)
        Theta, Phi = np.meshgrid(np.arccos(x), phi)

        f = mode_mapping(Theta, Phi, best_fit, mapping, l_max)
        g = spheroidal(Theta, Phi, mapping, l_max, sim.chif_mag)

        numerator = np.abs(np.sum(f*np.conj(g)))
        denominator = np.sqrt( np.abs(np.sum(f*np.conj(f))) * np.abs(np.sum(g*np.conj(g))) )

        match = numerator/denominator

        spatial_mismatch[i] = 1-match

    return spatial_mismatch
