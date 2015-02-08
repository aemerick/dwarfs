from __future__ import division

import numpy as np
import cgs as cgs
import matplotlib.pyplot as plt

def NFW_DM(r, r_s=None, c=12., M200=1.0E12*cgs.Msun, rho_crit = 9.74E-30):
    """
        Compute the dark matter density profile for a 
        NFW halo. Taking the scale as R200/c

        Parameters
        -----------
        r : array
            Radius from center of halo in cm
        r_s : float, optional
            The scale paramter. If none is provided, it is computed 
            as R200 / c. If r_s is provided, c is ignored. Default = None
        c : float, optional
            Concentration parameter of the DM halo. This is 
            ignored if r_s is provided. Default c=12 (MW)
        M200 : float, optional
            Virial DM mass of halo. Default 1.0E12 Msun (MW)
        rho_crit : float, optional
            Critical density of the universe. This is a function of 
            redshift. Default 9.74E-30 (z = 0)
    """

    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    # scale radius is virial divided by c
    if r_s == None:
        r_s = R200 / c
    else:
        c = R200 / r_s

    # the scale density is given as
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    return rho_s / ( (r/r_s)*(1+r/r_s)**2)
    
    
    
def NFW_isothermal_gas(r, r_s=800.0*cgs.pc, c=None, M200=2.0E7*cgs.Msun,
                          T = 1.0E4, n_o = 0.27,  mu = 1.3, rho_crit=9.74E-30):
    """
        The isothermal density profile for a gas in HSE with a 
        NFW DM halo. The default parameters are Sextans from
        Gatto et. al. 2013
    """
    # calculate the virial radius
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    if r_s == None and c == None:
        print "must give r_s (scale radius) OR c (concentration parameter)"
        return 0.0
    elif r_s == None:
        r_s = R200 / c
    elif c == None:
        c = R200 / r_s
    
    # scale density for the DM halo
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    # constant in the exponential
    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu *cgs.mp/(cgs.kb*T)

    # central mass density 
    rho_o = n_o * cgs.mp * mu

    # gas profile
    rho = rho_o * np.exp(-C_NFW * (1.0 - np.log(1.0+r/r_s)/(r/r_s)))

    return rho


def plot_profile(r, profile, filename=None, **kwargs):

    function_dict = {'Isothermal NFW': NFW_isothermal_gas}
    
    rho = NFW_isothermal_gas(r, **kwargs)

    plt.plot(r/cgs.pc, rho, label=profile,lw=1.75,color='black')

    plt.semilogy()
    plt.xlabel(r'r (pc)')
    plt.ylabel(r'$\rho$ (g cm$^{-3}$)')
    
    if filename == None:
        filename = profile + "_gas_profile.png"

    plt.savefig(filename)
    plt.close()



r = np.logspace(-2,3.5,1000.0)*cgs.pc
plot_profile(r,'Isothermal NFW')
