"""
galaxy

  currently stealing functions from M. Salem (msalem@astro.columbia.edu)
  will adapt better to personal use later
"""

import numpy as np
import cgs as cgs # units list from Munier

def NFW_mass(r, c = 12, M20 = 1.0E12*cgs.Msun, rho_crit=9.74E-30):
    """
       Spherical NFW DM profile

       Parameters
       ----------
       r  :  nump array or ndarray
             spherical radius from GC
       c  :  float, optional 
             NFW concentration parameter, default : 12
       M200 : float, optional
             NFW halo mass (cgs)
             default : 10^12 Msun
       rho_crit : float, optional 
             critical densit (cgs)
             default : 9.74e-30 g/cm^3
    """

    R200 = (3.0 * M200 / (4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
    fc   = np.log(1.0 + c) - c/(1.0+c)
    x    = r * c /R200

    return M200 / fc * ( np.log(1.0+x) - x/(1.0+x))

def halo_gas_ndensity(r,n0=(.46),rc=(.35)*cgs.kpc , beta=(.71), many_model=False):
    """
        Beta gas density profile for galactic halo in hydrostatic equilibrium
	Default fit from  Miller & Bregman 2013
		
	Parameters
	----------
	r : numpy array or ndarray
		spherical radius from galactic center (cgs)
	n0 : float, optional
		core number density (cgs)
		default : .46+.74
	rc : float, optional
		core radius (cgs)
		default : .35 + .29 kpc
	beta : float, optional
		density falloff exponent
		default : .71 - .14
    """
    
    f = lambda x: n0*(1.0 + (x/rc)**2)**(-3.0*beta/2.0)
        
    if many_model:
    # return several profiles over the range of Millerg and Bregman best 
    # fit parameters
    
        
    
        return all_n0, all_rc, all_beta, f(r) 
    
    else:
        return f(r)
        
def gato_ndensity(r, halo_type='isothermal'):
    """
        Cubic spline interpolation of digitized Fig. 8 from Gato et. al. 2013
        for the halo density. Reads data points from file and returns 
        the values. 
    
        r : numpy array or ndarray
            spherical radius from galactic center (cgs)
        halo_type : string, optional
            type of haly to choose from. Choices are 'isothermal',
            'adiabatic', or 'cooling'. Default: isothermal
    """
    file_dir = "/home/emerick/Research/dwarfs/analysis/halo/gato/"
    
    data_file = file_dir + "gato_" + halo_type + ".dat"
    
    data = np.genfromtxt(data_file, names=True)
    data['r'] = data['r'] * cgs.kpc # kpc to cm
    
    f = interp1d(data['r'], data['n'], kind='cubic')
    
    return f(r)
    
    
def halo_gas_mass(r, n, mu = cgs.mu):
    """
        Enclosed mass of a spherical halo of density n in Msun
        
        Parameters
        ----------
        r : array or ndarray
            radius in cgs units
        n : number density
        mu : float
            mean molecular weight. Default : cgs.mu = 0.61 
        
    """
    
    
    dr = r[1:] - r[:-1]
    density = n[1:] * mu * cgs.mp
    r       = r[1:]
    
    dV = 4.0 * np.pi * r**2 * dr
    
    return dV*density / cgs.Msun
       
def halo_gas_temperature(r, n, M = None, gamma = 1.6667, mu = cgs.mu):
    """
        Gas temperature for galactic halo in hydrostatic equillibrium.
        Calculates the temperature profile given either the number density
        profile, or just the mass profile.
        
        Parameters
        ----------
        r : numpy array or ndarray
            spherical radius from GC, cgs
        n : numpy array or ndarray
            density profile at each r. Mass enclosed is calculated from
            this. 
        M : optional, default none
            enclosed mass at each r. Used in place of calculated version
            if n is 
    """
    
    if M is None:
        M = halo_gas_mass(r, n, mu)
    
    return gamma * cgs.G * mu * cgs.mp * M / (3.0 * r * cgs.kb)
        

