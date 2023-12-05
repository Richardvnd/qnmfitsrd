""" 
Code to calculate the integral of three spherical harmonics 
"""

import json
import os
import spherical
import numpy as np 
import quaternionic
import spherical
from scipy.integrate import dblquad as dbl_integrate
from scipy.interpolate import interp1d
from scipy.special import sph_harm as Yml
from sympy.physics.wigner import wigner_3j, wigner_6j


if os.path.isfile('integrals_real.json') and os.path.isfile('integrals_imag.json'):
    with open('integrals_real.json', 'r') as f:
        betas_real = json.load(f)
    with open('integrals_imag.json', 'r') as f:
        betas_imag = json.load(f)

def sYlm(l, m, theta, phi, s=-2, l_max=9):
    wigner = spherical.Wigner(l_max)
    R = quaternionic.array.from_spherical_coordinates(theta, phi)
    Y = wigner.sYlm(s, R)
    return Y[wigner.Yindex(l, m)]

# See appendix D of 3+1 numerical relativity 

def triple_sph_int(l1, m1, l2, m2, l3, m3, s1=-2, s2=-2, s3=-2):
    return (((2*l1+1)*(2*l2+1))/(4*np.pi*(2*l3+1)))**0.5 * \
    spherical.clebsch_gordan(l1, s1, l2, s2, l3, -s3) * \
    spherical.clebsch_gordan(l1, m1, l2, m2, l3, m3) * \
    (-1)**m3 * \
    (-1)**(l1-l2+m3) *\
    (-1)**(l1+l2+s3) * \
    (1/(2*l3 + 1))
    

def calc_integrals(l1, m1, l2, m2, l3, m3):

    f_real = lambda theta, phi: np.real(
            np.sin(theta) * 
            sYlm(l1, m1, theta, phi) * 
            sYlm(l2, m2, theta, phi) * 
            np.conj(sYlm(l3, m3, theta, phi)))

    f_imag = lambda theta, phi: np.imag(
            np.sin(theta) * 
            sYlm(l1, m1, theta, phi) * 
            sYlm(l2, m2, theta, phi) * 
            np.conj(sYlm(l3, m3, theta, phi)))

    beta_real = dbl_integrate(f_real, 0, 2*np.pi, 0, np.pi)[0]
    beta_imag = dbl_integrate(f_imag, 0, 2*np.pi, 0, np.pi)[0]
    return beta_real + 1j * beta_imag

def get_integrals(l1, m1, l2, m2, l3, m3):
    return betas_real[f'{l1}{m1}{l2}{m2}{l3}{m3}'] + 1j * betas_imag[f'{l1}{m1}{l2}{m2}{l3}{m3}'] 

for a in range(2,5):
    for b in range(-a,a+1):
        for c in range(2,5):
            for d in range(-c,c+1):
                for e in range(2,5):
                    for f in range(-e,e+1):
                        print(a,b,c,d,e,f,triple_sph_int(a,b,c,d,e,f))

breakpoint() 
