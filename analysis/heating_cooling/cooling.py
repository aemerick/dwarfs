import numpy as np

def sw_dm(n,T, Tmin=6.0E4, Tmax=1.0E5):
    """
        Smoothly combining sarazin and white and dalgarno and mccray.
        Above 10E5 -> Sarazin and white. Below, Dalgarno.
    """

    SW =  sarazin_white(n,T)
    DM =  radloss(n,T)

    

    if np.size(SW) == 1:
        if (T < Tmax) and (T > Tmin):
            f = (T - Tmin)/(Tmax - Tmin)
            cooling = f * SW + (1-f)*DM
        elif (T >= Tmax):
            cooling = SW
        else:
            cooling = DM
    
        return cooling

    else:
        cooling = np.zeros(np.shape(T))

        f = (T - Tmin) / (Tmax - Tmin)

        f[f>1.] = 1.0
        f[f<0.] = 0.0
        cooling = f * SW + (1.0-f)*DM
        
        return cooling        

def sarazin_white(n, T):
    """
        Returns equation 27 in Sarazin and White 1987. This 
        is really only valid at 10**5 to 10**8 K

    Parameters:
    n : array
        Array of number densities (not used here... kept only for
        formatting)
    T : array
        Temperatures        
    """

    return 1.0E-22 *\
             4.700             * np.exp(-(T/3.5E5)**(4.5)) +\
             0.313   * T**(0.08) * np.exp(-(T/3.0E6)**(4.4)) +\
             6.420   * T**(-0.2) * np.exp(-(T/2.1E6)**(4.0)) +\
             4.39E-2 * T**(0.35)

def radloss(n,temperatures):
    """
    Dalgarno and McCray cooling curve with ionization fraction 
    of 0.01
    """

    def _DM(T):
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

        return loss



    if np.size(temperatures) > 1:
        i = 0
        radloss = np.zeros(np.shape(temperatures))
        
        imax, jmax = np.shape(temperatures)
        
        for i in np.arange(imax):
            for j in np.arange(jmax):
                loss = _DM(temperatures[i][j])
                radloss[i][j] = loss
    else:
        radloss = _DM(temperatures)


    return radloss
            
def timescale(n, T, gamma = 1.66666667):
    """
    returns the cooling timescale using radloss function
    """

    k = 1.380658E-16 # erg / K

    return (gamma - 1.0)**(-1.0) * k*T/(n * radloss(T))
    
def IIK_2007(n, T, Gamma=2.0E-26, mode='classic'):
    """
        Cooling curve used Inoue, Inutsuka & Koyama 2007, ApJL, 658, 99.
        This is what J. Grcevich used to establish equillibrium in Proeteus
        simulations of Leo T.
         
    """

    # L1 is curve from paper... not sure what L2 is... from J. Grcevich
    L1 =  Gamma*\
        (1.0E7*np.exp(-1.184E5/(T + 1000.)) + 1.4E-2*T**(0.5)*np.exp(-92./T))
    L2 = Gamma*\
        (1.0E4*np.exp(-1.2E7/(T+1.0E5)) + 1.75E-4*T**(0.5)*np.exp(-9.2E4/T))

    if mode == 'classic':
        total = L1

    elif mode == 'modified':
        total = L1 + L2

    elif mode == 'jana':
        # getcool function from J. Grcevich's code         
        mix1 = 0.5 * (1.0 + np.tanh((n-4.5-3.0)/5.2-4.0))
        mix2 = 0.5 * (1.0 - np.tanh((n-4.5-3.0)/5.2-4.0))

        total =  mix1*L1 + mix2*L2

    return total

