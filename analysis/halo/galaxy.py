"""
galaxy

  currently stealing functions from M. Salem (msalem@astro.columbia.edu)
  will adapt better to personal use later
"""

import numpy as np
import cgs as cgs # units list from Munier
from scipy import interpolate as interp

def NFW_mass(r, c = 12, M200 = 1.0E12*cgs.Msun, rho_crit=9.74E-30):
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

def halo_gas_ndensity(r,n0=(.46),rc=(.35)*cgs.kpc , beta=(.71), ambient=1.0E-5, many_model=False):
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
    ambient : float, optional
        ambient medium scaling used in Miller Bregman 2013 to account
        for RPS of dwarf spheroidals at large radii. Default 1.0E-5
    """
    
    f = lambda x: n0*(1.0 + (x/rc)**2)**(-3.0*beta/2.0) + ambient
        
    if many_model:
    # return several profiles over the range of Millerg and Bregman best 
    # fit parameters
    
        
    
        return all_n0, all_rc, all_beta, f(r) 
    
    else:
        return f(r)
        

#####
## DO NOT HAVE PROPER FIT PARAMETERS FOR RHO_O, PHI_O OR A
def gato_ndensity_analytic(r, a = 170.0*cgs.kpc, M=1.9E12*cgs.Msun,
                     gamma=None, halo_type='isothermal',
                     rho_o=1, phi_o=1, A=1, mu = cgs.mu):
    """
        returns the density of gas using the profile obtained 
        by Gato et. al. using the truncated flat potential for 
        the Milky way
    """

    if gamma == None:
        if halo_tye == 'isothermal':
            gamma = 1.0
        elif halo_type == 'adiabatic':
            gamma = 5.0/3.0
        elif halo_type == 'cooling':
            gamma = 1.33333333
    

    if gamma == 1.0:
        rho = rho_o * np.exp(- (TF_potential(r) - phi_o)/A)

    else:
        rho = rho_o * ((1.0 - TF_potential(r) - phi_o) *\
              (gamma-1)/(gamma*A))**(1.0/(gamma-1.0)))


    return rho / (cgs.mp * mu)

def TF_potential(r, M=1.9E12*cgs.Msun, a = 170.0*cgs.kpc):
    """
        Truncated flat dark matter potential. Default units are for 
        MW fit used in Gato et. al. 2013 from Lux et. al. 2010 and
        originally Wilkinson & Evans 1999
    """

    return (cgs.G*M/a) * np.log( ((r**2 + a**2)**0.5 + a) / r)


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
    
    f = interp.interp1d(data['r'], data['n'], kind='cubic')
    
    rmin = np.min(data['r'])
    rmax = np.max(data['r'])
    
    n = np.zeros(np.size(r))
    
    n[(r>rmin)*(r < rmax)] = f(r[(r>rmin)*(r<rmax)])
    n[r >=  rmax] = np.ones(np.size(n[r>=rmax]))*n[r<rmax][-1]
    n[r <=  rmin] = np.ones(np.size(n[r<=rmin]))*n[r>rmin][0]
    
    return n
    
    
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
       
def halo_gas_temperature(r, M = None, n=None, gamma = 1.6667, mu = cgs.mu,
                            **kwargs):
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
    
    if M == None:
        M = NFW_mass(r, **kwargs)
    
    if not n == None:
        Mgas = halo_gas_mass(r, n)
        if not (np.size(Mgas) == np.size(M)):
            Mgas = np.append(Mgas, Mgas[-1]) * cgs.Msun
        M = M + Mgas


    return gamma * cgs.G * mu * cgs.mp * M / (3.0 * r * cgs.kb)
        
def kaufmann(model, field):
    """
    Loads kaufman data.
    """
    file_dir = "/home/emerick/Research/dwarfs/analysis/halo/jana/"   
    data_file = file_dir + 'kaufmann' + model + ".dat"
    
    data = np.genfromtxt(data_file, names=True, skip_header=7)
    
    if field == 'density':
        field_data = data['3d_ghdens'] * 0.009420003
    
    elif field == 'r':
        field_data = data['radius']
    
    elif field == 'temperature':
        field_data = data['gh_temp']
    
    return field_data
    



