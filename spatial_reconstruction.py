"""

Updated: 13/12/2023

These are functions Chris wrote that compute the spatial distribution of the mapped mode (mode_mapping) 
and the expected distribution (spheroidal) based on the mixing coefficient. The final function (my_function)
computes the 'spatial mismatch'. 


"""

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
        if best_fit['modes'][loop]==mapping:
            A = best_fit['C'][loop]
            ans += A * Y[:,:,wigner.Yindex(*best_fit['spherical_modes'][i])]
            i += 1
    ans /= np.max(np.abs(ans)) # normalise peak value
    return ans


def spheroidal(theta, phi, mapping, l_max, chif):
    wigner = spherical.Wigner(l_max)
    R = quaternionic.array.from_spherical_coordinates(theta, phi)
    Y = wigner.sYlm(-2, R)

    ans = np.zeros_like(theta, dtype=complex)
    if len(mapping)==4:
        l, m, n, p = mapping
        for lp in np.arange(2, l_max+1):
            ans += qnmfits.qnm.mu(lp, m, l, m, n, p, chif) * Y[:,:,wigner.Yindex(lp, m)]
    elif len(mapping)==8:
        a, b, c, sign1, e, f, g, sign2 = mapping
        j = b + f
        for i in np.arange(2, l_max+1):
            ans += qnmfits.qnm.alpha([(i,j)+mapping], chif) * Y[:,:,wigner.Yindex(i,j)]

    ans /= np.max(np.abs(ans)) # normalise peak value
    return ans


def spatial_mismatch(f, g, num_points=100):

    # TODO: switch to dbl quadrature 

    dx, dphi = 2./num_points, 2*np.pi/num_points
    x = np.arange(-1, 1, dx)
    phi = np.arange(-np.pi, np.pi, dphi)
    Theta, Phi = np.meshgrid(np.arccos(x), phi)

    z = np.sum(f*np.conj(g))

    numerator = np.abs(z)
    denominator = np.sqrt( np.abs(np.sum(f*np.conj(f))) * np.abs(np.sum(g*np.conj(g))) )

    arg = np.angle(z)

    spatial_mismatch = 1 - numerator/denominator

    return spatial_mismatch, arg, z
