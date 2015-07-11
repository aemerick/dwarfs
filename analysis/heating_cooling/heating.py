from yt import units as u
import cooling as cool
import numpy as np
import cgs as cgs

def LTE_heating_fit(r, density, T, order = 8,
                                outfile = 'HSE_heating_fit.dat'):
    """
    Supply density and temperature and radius. Fit profile using an
    nth order polynomial
    """
    
    # find the profile
    gamma, n_gamma, nn_delta = heating_balance(density, T, number_density=True)
    
    # find the polynomial fit
    pol = np.polyfit(r, gamma, order)
    p   = np.poly1d(pol)
    
    # flip around the polynomial order
    pol_flip = list(reversed(pol))
    
    # write out 
    f = open(outfile, 'w')
    for val in pol_flip:
        f.write("%16.16e\n"%(val))
        
    f.close()
    
    # return polynomials and profile
    return pol_flip, p(r)

def heating_balance(density, T, mu = 1.31, number_density = False):
    """
    calculates the heating balance needed for HSE against cooling
    for the density and temperature profiles given

    if density type == number, then density is taken as the number density
    """
    mh = 1.6733E-24 # mass of H in grams 

    if number_density:
        ndens = density
    else:
        ndens = density / (mh * mu)

    # n * Gamma - n*n*Lambda = 0 --- Heating balance
    Lambda = cool.radloss(ndens,T)
    Gamma = Lambda * ndens

    return Gamma, ndens*Gamma, ndens*ndens*Lambda


def metagalactic(n=0.,T=0.):
    """
    """

    return 0.889E-13 * 2.889E-13 # ev/s -> erg/s
    
def diffuse_UV(r, r_uv = 150.0*cgs.pc, pe_heat=2.0E-26):
    """
        Returns the diffuse UV heating rate.
        
        Parameters
        ----------
        r : array, ndarray
            Radial distance from dwarf center in cm
        r_uv : float, optional
            UV exponential scale in cm. Default: 150.0 pc
        pe_heat : float, optional
            UV heating rate. Defualt : 2.0E-26
    """
    
    return pe_heat * np.exp(- r / r_uv)

def IIK_2007(n, T, Gamma=2.0E-26):
   """
   
   """
   return Gamma
   
def lower_bound_heating(r, n, T, T_heat_min=10.0, T_heat_max=2.0E4, **kwargs):
    """
       The lower bound on the possible heating in the dwarf galaxy.
       This is taken as the HM12 metagalactic background + the diffuse_UV
       heating.
       
       Parameters
       ----------
       r : array, ndarra
           Radial distance from the center of the dwarf in cm
       T_heat_min : float, optional
           Heat for T > T_heat_min. Default: 0.0 K
       T_heat_max : float, optional
           Heat for T < T_heat_max. Default 2.0E4 K
           
    """
    if np.size(r) == 0:
        if r == 0.0:
            r = np.zeros(np.shape(T))
    if np.size(T) == 0 and np.size(r) > 0:
        T = np.ones(np.shape(r)) * T

    total_rate = np.zeros(np.shape(T))
    if np.size(T) > 1:
        select = (T > T_heat_min) * (T < T_heat_max)    
        total_rate[select] = diffuse_UV(r[select], **kwargs) + metagalactic()
    elif((T<T_heat_max) and (T>T_heat_min)):
        total_rate = diffuse_UV(r, **kwargs) + metagalactic()
    else:
        total_rate = 0.0


    return np.array(total_rate)
