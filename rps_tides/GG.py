import numpy as np
import cgs
import profiles as prof

from scipy.optimize import brentq

def stripping_radius(star, gas, Pram = None, rho = None, v = None,
                     rmax  = 40.0 * cgs.kpc,
                     alpha = 1.0, lu = 1.0): 
    """
    Given a stellar and gas surface density profile, which
    both must be functions of radius alone, evaluate the stripping
    radius for a given ram pressure stripping force, provided 
    either as P_ram or rho and v kwargs. All units are assumed
    to be in cgs. Stripping radius can optionally be converted to
    one's units of choice provided a conversion factor.
    """

    if Pram is None and ( (rho is None) or (v is None)):
        print "Error, must provide a density and velocity, OR Pram"
        raise ValueError
    elif Pram is None:
        Pram = rho * v * v

    # create the function to root solve
    func = lambda x : 2.0*np.pi*cgs.G*( star(x) + alpha*gas(x))*gas(x) - Pram

    # completely stripped
    if func(0) <= 0.0:
        return 0.0

    Rstrip = brentq(func, 1.0E-10*rmax, rmax)

    if Rstrip <= 5.0E-10*rmax:
        Rstrip = 0.0

    return Rstrip / lu
    

if __name__ == "__main__":
    # try with all of the defaults


    sigma_star = lambda x : prof.stellar_surface_density(x)
    sigma_gas  = lambda x : prof.gas_surface_density(x)

    Pram = 1.4E-13 # max of Salem et. al. 2015 LMC RPS test

    Rstrip = stripping_radius(sigma_star, sigma_gas, Pram = Pram, lu = cgs.kpc)

    print Rstrip, ' kpc'




