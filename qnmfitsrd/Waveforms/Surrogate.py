import numpy as np
import spherical

from .Base import BaseClass


class NRSur7dq4(BaseClass):
    """
    A class for the NRSur7dq4 surrogate model presented in Varma et al. 2019,
    arXiv:1905.09300.
    
    Evaluates gravitational waveforms generated by precessing binary black 
    hole systems with generic mass ratios and spins.
    
    This model includes the following spin-weighted spherical harmonic modes:
    2<=ell<=4, -ell<=m<=ell.
    
    The parameter space of validity is:
    q \in [1, 6], and |chi1|,|chi2| \in [-1, 1], with generic directions.
    where q is the mass ratio and chi1/chi2 are the spin vectors of the
    heavier/lighter BH, respectively.
    
    The surrogate has been trained in the range q \in [1, 4] and 
    |chi1|/|chi2| \in [-0.8, 0.8], but produces reasonable waveforms in the
    above range and has been tested against existing NR waveforms in that
    range.

    Parameters
    ----------
    q : float, optional
        The mass ratio of the system. The default is 1.
        
    chi1 : array_like, optional
        Dimensionless spin vector of the heavier black hole. For the 
        reference time we use the earliest available time, ~4000M before
        merger. The default is [0,0,0].
        
    chi2 : array_like, optional
        The dimensionless spin vector of the lighter black hole. The default
        is [0,0,0].
        
    f_ref : float, optional
        Frequency used to set the reference epoch at which the reference frame 
        is defined and the spins are specified. Should be in cycles/M. The 
        default is 0.01. 
 
        The reference frame (or inertial frame) is defined as follows. The +ve 
        z-axis is along the orbital angular momentum at the reference epoch. 
        The separation vector from the lighter BH to the heavier BH at the 
        reference epoch is along the +ve x-axis. The y-axis completes the 
        right-handed triad.
        
    ellMax : int, optional
        Maximum ell index for modes to include. All available m indicies 
        for each ell will be included automatically. The default is None, 
        in which case all available modes will be included.
        
    zero_time : optional
        The method used to determine where the time array equals zero. 
        
        Options are:
            
        - time : float
            The time on the default SXS simulation time array where we set t=0.
        - (l,m) : tuple
            The peak of the absolute value of this mode is where t=0.
        - 'norm'
            The time of the peak of total amplitude (e.g. Eq. 2 of 
            https://arxiv.org/abs/1705.07089 ) is used.
        - 'Edot'
            The time of peak energy flux is used.
          
        The default is 0 (the default simulation zero time).
        
    transform : str or list, optional
        Transformations to apply to the SXS data. Options are:
            
        - 'rotation'
            Reperform the spin-weighted spherical harmonic decomposition in a 
            rotated coordinate system where the z-axis is parallel to the 
            remnant black hole spin.
        - 'dynamic_rotation'
            Reperform the spin-weighted spherical harmonic decomposition in a 
            rotated coordinate system where the z-axis is parallel to the 
            instantanious black hole spin.
            
        The default is None (no tranformations are applied).
    """
    
    def __init__(self, q=1, chi1=[0,0,0], chi2=[0,0,0], f_ref=0.01,
                 ellMax=None, zero_time=0, transform=None):
        """
        Initialize the class.
        """
        self.q = q
        self.chi1 = chi1
        self.chi2 = chi2
        self.f_ref = f_ref
        self.ellMax = ellMax
        self.zero_time = zero_time
        
        # Individual mass values are useful
        self.m1 = self.q/(1+self.q)
        self.m2 = 1/(1+self.q)
        self.M = self.m1 + self.m2
        
        # Load the surrogate
        # ------------------
        
        # The surrogate models
        import gwsurrogate as gws

        # Download the surrogate if it hasn't been already - this saves the 
        # surrogate .h5 file to the 
        # lib/python/site-packages/gwsurrogate/surrogate_downloads folder
        if 'NRSur7dq4' not in dir(gws):
            gws.catalog.pull('NRSur7dq4')
        
        # Load the surrogate
        sur = gws.LoadSurrogate('NRSur7dq4')
        
        # Evaluate the surrogate with the provided parameters. f_low=0 returns 
        # the whole waveform.
        self.times, self.h, self.dyn = sur(
            q=q, chiA0=chi1, chiB0=chi2, f_low=0, f_ref=f_ref, ellMax=ellMax, 
            precessing_opts={'return_dynamics': True})
        
        # The ellMax value will be useful
        if self.ellMax == None:
            self.ellMax = 4
        
        # Estimate the remnant properties
        # -------------------------------
        
        # We use the surfinBH package to calculate the remnant black hole 
        # properties, making sure to use consistant values for the reference
        # epoch.
        import surfinBH as gwsremnant

        # Load the surrogate remnant
        surrem = gwsremnant.LoadFits('NRSur7dq4Remnant')
        
        # The remnant mass and 1-sigma error estimate
        self.Mf, self.Mf_err = surrem.mf(
            q, self.chi1, self.chi2, omega0=np.pi*self.f_ref)
        
        # The remnant spin and 1-sigma error estimate
        self.chif, self.chif_err = surrem.chif(
            q, self.chi1, self.chi2, omega0=np.pi*self.f_ref)
        self.chif_mag = np.linalg.norm(self.chif)
        
        # Angular coordinates of the final spin vector
        chif_norm = self.chif/self.chif_mag
        self.thetaf = np.arccos(chif_norm[2])
        self.phif = np.arctan2(chif_norm[1], chif_norm[0])
        
        # Frame independent flux quantities
        # ---------------------------------
        
        # Calculate waveform mode time derivatives
        self.calculate_hdot()
        
        # Calculate energy flux
        self.calculate_Moft()
        
        # Calculate angular momentum flux
        self.calculate_chioft()
        
        # Frame transformations
        # ---------------------
        
        # Shift the time array to use the requested zero-time.
        self.time_shift()
        
        # Construct a Wigner object
        self.wigner = spherical.Wigner(self.ellMax)
        
        # Apply tranformations
        if type(transform) != list:
            transform = [transform]
        
        for transformation in transform:
            if transformation == 'rotation':
                self.rotate_modes()
            elif transformation == 'dynamic_rotation':
                self.rotate_modes_over_time()
            elif transformation == 'boost':
                pass
            elif transformation is None:
                pass
            else:
                print('Requested transformation not available.')
        
        # Other interesting quantities
        # ----------------------------
        
        # Calculate the approximate frequency evolution
        self.calculate_foft()
        
        
class NRHybSur3dq8(BaseClass):
    """
    A class for the NRHybSur3dq8 surrogate model presented in Varma et al. 
    2018, arxiv:1812.07865.
    
    Evaluates gravitational waveforms generated by aligned-spin binary black 
    hole systems. This model was built using numerical relativity (NR) 
    waveforms that have been hybridized using post-Newtonian (PN) and 
    effective one body (EOB) waveforms.
    
    This model includes the following spin-weighted spherical harmonic modes:
    (2,2), (2,1), (2,0), (3,3), (3,2), (3,1), (3,0), (4,4) (4,3), (4,2) and 
    (5,5). The m<0 modes are deduced from the m>0 modes.
    
    The parameter space of validity is:
    q \in [1, 10] and chi1z/chi2z \in [-1, 1],
    where q is the mass ratio and chi1z/chi2z are the spins of the 
    heavier/lighter BH, respectively, in the direction of orbital angular 
    momentum.
    
    The surrogate has been trained in the range q \in [1, 8] and chi1z/chi2z 
    \in [-0.8, 0.8], but produces reasonable waveforms in the above range and 
    has been tested against existing NR waveforms in that range.

    Parameters
    ----------
    q : float, optional
        The mass ratio of the system. The default is 1.
        
    chi1 : array_like, optional
        Dimensionless spin vector of the heavier black hole. For the 
        reference time we use the earliest available time, ~4000M before
        merger. The default is [0,0,0].
        
    chi2 : array_like, optional
        The dimensionless spin vector of the lighter black hole. The 
        default is [0,0,0].
        
    f_ref : float, optional
        Frequency used to set the reference epoch at which the reference frame 
        is defined and the spins are specified. Should be in cycles/M. The 
        default is 0.01. 
 
        The reference frame (or inertial frame) is defined as follows. The +ve 
        z-axis is along the orbital angular momentum at the reference epoch. 
        The separation vector from the lighter BH to the heavier BH at the 
        reference epoch is along the +ve x-axis. The y-axis completes the 
        right-handed triad.
        
    ellMax : int, optional
        Maximum ell index for modes to include. All available m indicies 
        for each ell will be included automatically. The default is None, 
        in which case all available modes will be included.
        
    zero_time : optional
        The method used to determine where the time array equals zero. 
        
        Options are:
            
        - time : float
            The time on the default SXS simulation time array where we set t=0.
        - (l,m) : tuple
            The peak of the absolute value of this mode is where t=0.
        - 'norm'
            The time of the peak of total amplitude (e.g. Eq. 2 of 
            https://arxiv.org/abs/1705.07089 ) is used.
        - 'Edot'
            The time of peak energy flux is used.
          
        The default is 0 (the default simulation zero time).
        
    transform : str or list, optional
        Transformations to apply to the SXS data. Options are:
            
        - 'rotation'
            Reperform the spin-weighted spherical harmonic decomposition in a 
            rotated coordinate system where the z-axis is parallel to the 
            remnant black hole spin.
        - 'dynamic_rotation'
            Reperform the spin-weighted spherical harmonic decomposition in a 
            rotated coordinate system where the z-axis is parallel to the 
            instantanious black hole spin.
            
        The default is None (no tranformations are applied).
    """
    
    def __init__(self, q=1, chi1=[0,0,0], chi2=[0,0,0], f_ref=0.01, 
                 ellMax=None, zero_time=None, inclination=None, phi_ref=0,
                 transform=None, ):
        """
        Initialize the class.
        """
        self.q = q
        self.chi1 = chi1
        self.chi2 = chi2
        self.f_ref = f_ref
        self.ellMax = ellMax
        self.zero_time = zero_time
        
        # Individual mass values are useful
        self.m1 = self.q/(1+self.q)
        self.m2 = 1/(1+self.q)
        self.M = self.m1 + self.m2
        
        # Load the surrogate
        # ------------------
        
        # The surrogate models
        import gwsurrogate as gws

        # Download the surrogate if it hasn't been already - this saves the 
        # surrogate .h5 file to the 
        # lib/python/site-packages/gwsurrogate/surrogate_downloads folder
        if 'NRHybSur3dq8' not in dir(gws):
            gws.catalog.pull('NRHybSur3dq8')
        
        # Load the surrogate
        sur = gws.LoadSurrogate('NRHybSur3dq8')
        
        # Evaluate the surrogate with the provided parameters. f_low=0 returns 
        # the whole waveform.
        self.times, self.h, self.dyn = sur(
            q=q, chiA0=chi1, chiB0=chi2, f_low=0, f_ref=f_ref, ellMax=ellMax)
        
        # The ellMax value will be useful (although this surrogate includes 
        # the (5,5) mode, no other ell=5 modes are included so we ignore for
        # simplicity)
        if self.ellMax == None:
            self.ellMax = 4
            
        # We need to manually populate the negative m modes (which are related
        # to the positive m modes by symmetry), and also fill the (4,0) mode
        # with zeros, as this isn't given
        for l in range(2, self.ellMax+1):
            for m in range(-l,l+1):
                
                if (l == 4) & (m == 0):
                    self.h[l,m] = np.zeros_like(self.times)
                elif m < 0:
                    self.h[l,m] = (-1)**l * np.conjugate(self.h[l,-m])
                    
        # Estimate the remnant properties
        # -------------------------------
        
        # We use the surfinBH package to calculate the remnant black hole 
        # properties. Note, for this remnant model we have to specify spins at
        # t=-100M from the peak of the waveform, can we be more careful with
        # this?
        import surfinBH as gwsremnant

        # Load the surrogate remnant
        surrem = gwsremnant.LoadFits('NRSur3dq8Remnant')
        
        # The remnant mass and 1-sigma error estimate
        self.Mf, self.Mf_err = surrem.mf(q, self.chi1, self.chi2)
        
        # The remnant spin and 1-sigma error estimate
        self.chif, self.chif_err = surrem.chif(q, self.chi1, self.chi2)
        self.chif_mag = np.linalg.norm(self.chif)
        
        # Angular coordinates of the final spin vector
        chif_norm = self.chif/self.chif_mag
        self.thetaf = np.arccos(chif_norm[2])
        self.phif = np.arctan2(chif_norm[1], chif_norm[0])
        
        # Frame independent flux quantities
        # ---------------------------------
        
        # Calculate waveform mode time derivatives
        self.calculate_hdot()
        
        # Calculate energy flux
        self.calculate_Moft()
        
        # Calculate angular momentum flux
        self.calculate_chioft()
        
        # Frame transformations
        # ---------------------
        
        # Shift the time array to use the requested zero-time.
        self.time_shift()
        
        # Construct a Wigner object
        self.wigner = spherical.Wigner(self.ellMax)
        
        # Apply tranformations
        if type(transform) != list:
            transform = [transform]
        
        for transformation in transform:
            if transformation == 'rotation':
                self.rotate_modes()
            elif transformation == 'dynamic_rotation':
                self.rotate_modes_over_time()
            elif transformation == 'boost':
                pass
            elif transformation is None:
                pass
            else:
                print('Requested transformation not available.')
        
        # Other interesting quantities
        # ----------------------------
        
        # Calculate the approximate frequency evolution
        self.calculate_foft()