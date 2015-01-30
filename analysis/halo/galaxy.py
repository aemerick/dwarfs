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

def halo_gas_density(r,n0=(.46+.74),rc=(.35+.29)*cgs.kpc ,beta=(.71-.14)):
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
    return cgs.mp*cgs.mu*n0*(1.0 + (r/rc)**2)**(-3.0*beta/2.0)
