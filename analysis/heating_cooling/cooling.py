import numpy as np


def radloss(temperatures):
    """
    Dalgarno and McCray cooling curve with ionization fraction 
    of 0.01
    """

    i = 0
    radloss = np.zeros(np.size(temperatures))

    for T in temperatures:
    

        if (T >= 1.70e4):
            if (T >= 5.62E5):
                if (T >= 2.76E6):
                    if (T >= 3.16E7):
                        loss = 3.090E-27 * T**0.5
                    else:
                        loss = 5.188E-21 / (T**0.33)
                else:
                    if ( T >= 1.78E6 ):
                        loss = 3.89E-4 / (T**(2.95))
                    else:
                        loss = 1.3E-22 * (T/1.5E5)**(0.01)
            else:
                if (T >= 7.94e4): 
                    if (T >= 2.51E5):
                        loss = 3.98E-11 / (T**2)
                    else:
                        loss = 6.31E-22*(T/1.0E-6)**(0.01)
                else:
                    if (T >= 3.98E4):
                        loss = 1.0E-31 * T **2
                    else:
                        loss = 1.479E-21/(T**0.216)  
        else:
            if (T >= 1.0E3):
                if (T >= 6.31E3):
                    if (T >= 1.0E4): 
                        loss = 7.63E-81*(T**13.8)
                    else:
                        loss = 3.13E-27*(T**0.396)
                else:
                    if ( T >= 3.16E3):
                        loss = 2.64E-26*(T**0.152)
                    else:
                        loss = 5.28E-27*(T**0.352)
            else:
                if (T >= 3.98E1):
                    if (T >= 2.00E2):
                        loss = 3.06E-27*(T**0.431)
                    else:
                        loss = 1.52E-28*(T**0.997)
                else:
                    if (T >= 2.51E1):
                        loss = 2.39E-29*(T**1.5)
                    else:
                        loss = 1.095E-32*(T**3.885)


        if (T <= 1.e1):
            loss = 0.0E-28


        radloss[i] = loss
        i = i + 1



    return radloss
            
def timescale(T, n, gamma = 1.66666667):
    """
    returns the cooling timescale using radloss function
    """

    k = 1.380658E-16 # erg / K

    return (gamma - 1.0)**(-1.0) * k*T/(n * radloss(T))
