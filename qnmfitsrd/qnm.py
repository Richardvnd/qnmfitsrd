import json
import math
import os.path
from itertools import product

import numpy as np
import qnm as qnm_loader
import quaternionic
import spherical
from scipy.integrate import dblquad as dbl_integrate
from scipy.interpolate import interp1d
from scipy.special import sph_harm as Yml


class qnm:
    """
    Class for loading quasinormal mode (QNM) frequencies and spherical-
    spheroidal mixing coefficients. This makes use of the qnm package,
    https://arxiv.org/abs/1908.10377
    """
    
    def __init__(self, l_max=9):
        """
        Initialise the class.
        """
        
        # Dictionary to store the qnm functions
        self.qnm_funcs = {}
        
        # Dictionary to store interpolated qnm functions for quicker 
        # evaluation
        self.interpolated_qnm_funcs = {}

        # For integration of spin-weighted spherical harmonics
        self.l_max = l_max
        self.wigner = spherical.Wigner(self.l_max)

        # Load files once when class is initialised
        if os.path.isfile('integrals_real.json') and os.path.isfile('integrals_imag.json'):
            with open('integrals_real.json', 'r') as f:
                self.betas_real = json.load(f)
            with open('integrals_imag.json', 'r') as f:
                self.betas_imag = json.load(f)
        
    def interpolate(self, l, m, n):
        
        qnm_func = self.qnm_funcs[l,m,n]
        
        # Extract relevant quantities
        spins = qnm_func.a
        real_omega = np.real(qnm_func.omega)
        imag_omega = np.imag(qnm_func.omega)
        all_real_mu = np.real(qnm_func.C)
        all_imag_mu = np.imag(qnm_func.C)

        # Interpolate omegas
        real_omega_interp = interp1d(
            spins, real_omega, kind='cubic', bounds_error=False, 
            fill_value=(real_omega[0],real_omega[-1]))
        
        imag_omega_interp = interp1d(
            spins, imag_omega, kind='cubic', bounds_error=False, 
            fill_value=(imag_omega[0],imag_omega[-1]))
        
        # Interpolate mus
        mu_interp = []
        
        for real_mu, imag_mu in zip(all_real_mu.T, all_imag_mu.T):
            
            real_mu_interp = interp1d(
                    spins, real_mu, kind='cubic', bounds_error=False, 
                    fill_value=(real_mu[0],real_mu[-1]))
                
            imag_mu_interp = interp1d(
                spins, imag_mu, kind='cubic', bounds_error=False, 
                fill_value=(imag_mu[0],imag_mu[-1]))
            
            mu_interp.append((real_mu_interp, imag_mu_interp))

        # Add these interpolated functions to the frequency_funcs dictionary
        self.interpolated_qnm_funcs[l,m,n] = [
            (real_omega_interp, imag_omega_interp), mu_interp]
        
    def omega(self, l, m, n, sign, chif, Mf=1):
        """
        Return a complex frequency, :math:`\omega_{\ell m n}(M_f, \chi_f)`,
        for a particular mass, spin, and mode. One or both of chif and Mf can
        be array_like, in which case an ndarray of complex frequencies is
        returned.
        
        Parameters
        ----------
        l : int
            The angular number of the mode.
            
        m : int
            The azimuthal number of the mode.
            
        n : int
            The overtone number of the mode.
            
        sign : int
            An integer with value +1 or -1, to indicate the sign of the real
            part of the frequency. This way any regular (+1) or mirror (-1)
            mode can be requested. Alternatively, this can be thought of as
            prograde (sign = sgn(m)) or retrograde (sign = -sgn(m)) modes.
            
        chif : float or array_like
            The dimensionless spin magnitude of the black hole.
            
        Mf : float or array_like, optional
            The mass of the final black hole. This is the factor which the QNM
            frequencies are divided through by, and so determines the units of 
            the returned quantity. 
            
            If Mf is in units of seconds, then the returned frequency has units 
            :math:`\mathrm{s}^{-1}`. 
            
            When working with SXS simulations and GW surrogates, we work in 
            units scaled by the total mass of the binary system, M. In this 
            case, providing the dimensionless Mf value (the final mass scaled 
            by the total binary mass) will ensure the QNM frequencies are in 
            the correct units (scaled by the total binary mass). This is 
            because the frequencies loaded from file are scaled by the remnant 
            black hole mass (Mf*omega). So, by dividing by the remnant black 
            hole mass scaled by the total binary mass (Mf/M), we are left with
            Mf*omega/(Mf/M) = M*omega.
            
            The default is 1, in which case the frequencies are returned in
            units of the remnant black hole mass.
            
        Returns
        -------
        complex or ndarray
            The complex QNM frequency or frequencies of length 
            max(len(chif), len(Mf)).
        """
        # Load the correct qnm based on the type we want
        m *= sign
        
        # Test if the qnm function has been loaded for the requested mode
        if (l,m,n) not in self.qnm_funcs:
            self.qnm_funcs[l,m,n] = qnm_loader.modes_cache(-2, l, m, n)
            
        if type(chif) in [float, np.float64]:
            omega, A, mu = self.qnm_funcs[l,m,n](chif, store=True)
            
        else:
            # Test if the interpolated qnm function has been created (we create 
            # our own interpolant so that we can evaluate the frequencies for 
            # all spins simultaneously)
            if (l,m,n) not in self.interpolated_qnm_funcs:
                self.interpolate(l,m,n)
                
            omega_interp = self.interpolated_qnm_funcs[l,m,n][0]
            omega = omega_interp[0](chif) + 1j*omega_interp[1](chif)
            
        # Use symmetry properties to get the mirror mode, if requested
        if sign == -1:
            omega = -np.conjugate(omega)
        
        return omega/Mf
    
    def omega_list(self, modes, chif, Mf=1):
        """
        Return a frequency list, containing frequencies corresponding to each
        mode in the modes list (for a given mass and spin).
        
        Parameters
        ----------            
        modes : array_like
            A sequence of (l,m,n,sign) tuples to specify which QNMs to load 
            frequencies for. For nonlinear modes, the tuple has the form 
            (l1,m1,n1,sign1,l2,m2,n2,sign2,...).
            
        chif : float or array_like
            The dimensionless spin magnitude of the final black hole.
            
        Mf : float or array_like, optional
            The mass of the final black hole. See the qnm.omega docstring for
            details on units. The default is 1.
            
        Returns
        -------
        list
            The list of complex QNM frequencies.
        """
        # For each mode, call the qnm function and append the result to the
        # list
        
        # Code for linear QNMs:
        # return [self.omega(l, m, n, sign, chif, Mf) for l, m, n, sign in modes]
        
        # Code for nonlinear QNMs:
        return [
            sum([self.omega(l, m, n, sign, chif, Mf) 
                 for l, m, n, sign in [mode[i:i+4] for i in range(0, len(mode), 4)]
                 ]) 
            for mode in modes
            ]
    
        # Writen out, the above is doing the following:
            
        # return_list = []
        # for mode in modes:
        #     sum_list = []
        #     for i in range(0, len(mode), 4):
        #         l, m, n, sign = mode[i:i+4]
        #         sum_list.append(self.omega(l, m, n, sign, chif, Mf))
        #     return_list.append(sum(sum_list))
        # return return_list
        
    def mu(self, l, m, lp, mp, nprime, sign, chif):
        """
        Return a spherical-spheroidal mixing coefficient, 
        :math:`\mu_{\ell m \ell' m' n'}(\chi_f)`, for a particular spin and 
        mode combination. The indices (l,m) refer to the spherical harmonic. 
        The indices (l',m',n') refer to the spheroidal harmonic. The spin chif
        can be a float or array_like.

        Parameters
        ----------
        l : int
            The angular number of the spherical-harmonic mode.
            
        m : int
            The azimuthal number of the spherical-harmonic mode.
            
        lp : int
            The angular number of the spheroidal-harmonic mode.
            
        mp : int
            The azimuthal number of the spheroidal-harmonic mode.
            
        nprime : int
            The overtone number of the spheroidal-harmonic mode.
            
        sign : int
            An integer with value +1 or -1, to indicate the sign of the real
            part of the QNM frequency. If the mixing coefficient associated
            with a -1 QNM (i.e. a mirror mode) is requested, then symmetry 
            properties are used for the calculation.
            
        chif : float or array_like
            The dimensionless spin magnitude of the final black hole.

        Returns
        -------
        complex or ndarray
            The spherical-spheroidal mixing coefficient.
        """
        # There is no overlap between different values of m
        if mp != m:
            return 0
        
        # Load the correct qnm based on the type we want
        m *= sign
        mp *= sign
        
        # Our functions return all mixing coefficients with the given 
        # (l',m',n'), so we need to index it to get the requested l
        if abs(m) > 2:
            index = l - abs(m)
        else:
            index = l - 2
        
        # Test if the qnm function has been loaded for the requested mode
        if (lp,mp,nprime) not in self.qnm_funcs:
            self.qnm_funcs[lp,mp,nprime] = qnm_loader.modes_cache(
                -2, lp, mp, nprime)
            
        if type(chif) in [float, np.float64]:
            
            # Access the relevant functions from the qnm_funcs dictionary, and 
            # evaluate at the requested spin. Storing speeds up future 
            # evaluations.
            omega, A, mu = self.qnm_funcs[lp,mp,nprime](chif, store=True)
            mu = mu[index]
            
        else:
            
            if (lp,mp,nprime) not in self.interpolated_qnm_funcs:
                self.interpolate(lp,mp,nprime)
                
            mu_interp = self.interpolated_qnm_funcs[lp,mp,nprime][1][index]
            mu = mu_interp[0](chif) + 1j*mu_interp[1](chif)
            
        # Use symmetry properties to get the mirror mixing coefficient, if 
        # requested
        if sign == -1:
            mu = (-1)**(l+lp)*np.conjugate(mu)
            
        return mu
        
    def mu_list(self, indices, chif):
        """
        Return a list of mixing coefficients, for all requested indices. See
        the qnm.mu() docstring for more details.
        
        Parameters
        ----------
        indices : array_like
            A sequence of (l,m,l',m',n',sign) tuples specifying which mixing 
            coefficients to return.
            
        chif : float
            The dimensionless spin magnitude of the final black hole.

        Returns
        -------
        mus : list
            The list of spherical-spheroidal mixing coefficients.
        """
        # List to store the mixing coeffs
        mus = []
        
        # For each mode, call the qnm function and append the result to the 
        # list
        for l, m, lp, mp, nprime, sign in indices:
            mus.append(self.mu(l, m, lp, mp, nprime, sign, chif))
        
        return mus
    
    def sYlm(self, l, m, theta, phi, s=-2):
        R = quaternionic.array.from_spherical_coordinates(theta, phi)
        Y = self.wigner.sYlm(s, R)
        return Y[self.wigner.Yindex(l, m)]
    

    def get_betas(self):
        """
        Calculates the integral of three spherical harmonics to determine 
        'beta' coefficients that can be used to calculate the mixing coefficients
        for quadratic modes. 

        The code calculates these once up to l=9 and saves then as a h5 locally. It 
        checks for this file before attempting to caclulate. 
    
        Returns
        -------
        betas : dictionary
            A dictionary of beta values labelled by their indices (i.e. beta_{dbhfij}).

        """

        if hasattr(self, 'betas_real') and hasattr(self, 'betas_imag'):
            betas_real = self.betas_real
            betas_imag = self.betas_imag
        else:
            betas_real = {}
            betas_imag = {}
            for d, b, h, f, i, j in product(range(2, self.l_max+1), range(-self.l_max, self.l_max+1), repeat=3):

                if (f'{h}{f}{d}{b}{i}{j}') in betas_real:
                    betas_real[f'{d}{b}{h}{f}{i}{j}'] = betas_real[f'{h}{f}{d}{b}{i}{j}']
                    betas_imag[f'{d}{b}{h}{f}{i}{j}'] = betas_imag[f'{h}{f}{d}{b}{i}{j}']
                    continue

                if (b + f - j) != 0:
                    betas_real[f'{d}{b}{h}{f}{i}{j}'] = 0
                    betas_imag[f'{d}{b}{h}{f}{i}{j}'] = 0
                    continue

                f_real = lambda theta, phi: np.real(
                                np.sin(theta) * 
                                self.sYlm(d, b, theta, phi) * 
                                self.sYlm(h, f, theta, phi) * 
                                np.conj(self.sYlm(i, j, theta, phi)))
        
                f_imag = lambda theta, phi: np.imag(
                                np.sin(theta) * 
                                self.sYlm(d, b, theta, phi) * 
                                self.sYlm(h, f, theta, phi) * 
                                np.conj(self.sYlm(i, j, theta, phi)))

                beta_real = dbl_integrate(f_real, 0, 2*np.pi, 0, np.pi)[0]
                beta_imag = dbl_integrate(f_imag, 0, 2*np.pi, 0, np.pi)[0]

                betas_real[f'{d}{b}{h}{f}{i}{j}'] = beta_real
                betas_imag[f'{d}{b}{h}{f}{i}{j}'] = beta_imag
            
            with open('integrals_real.json', 'w') as f:
                json.dump(betas_real, f)

            with open('integrals_imag.json', 'w') as f:
                json.dump(betas_imag, f)

        return betas_real, betas_imag

    def alpha(self, indices, chif):
        """
        Calculates alphas (equivalent to mus for quadratic modes).

        Parameters
        ----------
        indices : array_like
            A sequence of (i,j,a,b,c,sign1,e,f,g,sign2) tuples (corresponding to the spherical indices and two spheroidal indices respectively)
            i.e. (l, m, l', m', n', p', l'', m'', n'', p'').

        chif : float
            The dimensionless spin magnitude of the final black hole.
        
    
        Returns
        -------
        alphas : list
            The list of spherical-spheroidal mixing coefficients for quadratic modes. 

        """

        betas_real, betas_imag = self.get_betas()

        alphas = [sum(self.mu(d, b, a, b, c, sign1, chif) * 
                      self.mu(h, f, e, f, g, sign2, chif) * 
                      (betas_real[f'{d}{b}{h}{f}{i}{j}'] + 1j * betas_imag[f'{d}{b}{h}{f}{i}{j}']) 
                      for d in range(2, self.l_max+1) 
                      for h in range(2, self.l_max+1)) 
                      for i, j, a, b, c, sign1, e, f, g, sign2 in indices]

        return alphas