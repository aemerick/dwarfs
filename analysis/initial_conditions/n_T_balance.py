import numpy as np
from heating_cooling import cooling as cool
from heating_cooling import heating as heat
from scipy import optimize as opt

def find_equilibrium(n, Tmin=1.E2, Tmax=1.E7, tolerance=1.0E-8,
                      cooling_func=None, heating_func=None):
    """
       Find the equillibrium temperatures for a given n using 
       bisection search.
       
       Parameters
       ----------
       n  :  array, ndarray
             Array of number density values to find equil temperature
       Tmin : float, optional
             Minimum allowed temperature to search over. Default: 100 K
       Tmax : float, optional
             Maximum allowed temperature to search over. Default: 1.0E7K
       tolerance: float, optional
             Convergence tolerance. Default : 1.0E-8
       cooling_func: function, optional
           Supply the desired cooling function here. The default value is
           "None". If None is supplied, the cooling function from
           cool.radloss is used, the same as implemented in the FLASH code
       heating_func: function, optional
           Supply the desired heating function here. The default is "None".
           If None, minimum heating is used
       
    """
    
    if cooling_func == None:
        cooling_func = cool.radloss
    if heating_func == None:
        heating_func = heat.lower_bound_heating
    
    T_equil = np.zeros(np.size(n))
    
    for i in np.arange(np.size(n)):
        function = lambda T: heating_func(0.0,T) - n[i]*cooling_func(T)
        root = opt.brentq(function, Tmin, Tmax, xtol=tolerance)
        
        T_equil[i] = root
        
    return T_equil
    


# number densities
nmin, nmax = 1.0E-5, 1.0E3
Tmin, Tmax = 10.0  , 1.0E7

n = np.linspace(nmin, nmax, 1.0E3)

T_equil = find_equilibrium(n)

fig,ax1 = plt.subplots()

ax1.plot(n,T_equil,label="Temperature",lw=lw,color='black')
ax1.loglog()
ax1.set_xlabel(r'n (cm$^{-3}$)'); ax1.set_ylabel(r'T (K)')

ax2 = ax1.twinx()
ax2.plot(n, n*T_equil, label='Pressure',lw=lw,color='red')
ax2.set_ylabel(r'Pressure (cm$^{-3}$ K)')
ax2.semilogy()

plt.legend(loc='best',fancybox=True)
plt.savefig('nT_balance.png')
plt.close()


