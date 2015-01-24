from yt import units as u
import cooling as cool
import numpy as np

def heating_balance(density, T, mu = 1.0, number_density = False):
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
    Lambda = cool.radloss(T)
    Gamma = Lambda * ndens

    return Gamma, ndens*Gamma, ndens*ndens*Lambda


def metagalactic():
    """
    """

    return 0.889E-13 * 2.889E-13 # ev/s -> erg/s
