import numpy as np
from heating_cooling import cooling as cool
from heating_cooling import heating as heat
from scipy import optimize as opt
from plotting import plotTools as myplot
import matplotlib.pyplot as plt
lw = 1.75


def find_equilibrium(n, Tmin=1.E2, Tmax=1.E7, tolerance=1.0E-12,
                      cooling_func=None, heating_func=None, cf_kwargs={}):
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
        heating_func = lambda n,T: heat.lower_bound_heating(0.0,n,T)
    
    T_equil = np.zeros(np.size(n))
    
    for i in np.arange(np.size(n)):
        function = lambda T: -heating_func(n[i],T) + n[i]*cooling_func(n[i],T,**cf_kwargs)

        if(np.sign(function(Tmin)) == np.sign(function(Tmax))):
            root = None
        else:
            root = opt.bisect(function, Tmin, Tmax, xtol=tolerance)
        
        T_equil[i] = root
        
    return T_equil
    

def plot_balance_map(nmin=1.0E-5,nmax=1.0E3,Tmin=10.0,Tmax=1.0E7,npoints=1.0E3):

    n = np.linspace(np.log10(nmin),np.log10(nmax),npoints)
    T = np.linspace(np.log10(Tmin),np.log10(Tmax),npoints)



    func = lambda n,T: heat.lower_bound_heating(0.0,n,T) - n*cool.radloss(T)
    
    nmesh,Tmesh,cmesh = myplot.color_mesh(func, n, T, xlog=True, ylog=True)

    cmesh = cmesh / np.abs(cmesh)
    cmin,cmax = np.min(cmesh),np.max(cmesh)

   

    plt.pcolormesh(nmesh,Tmesh, cmesh, cmap='RdBu',vmin=cmin,vmax=cmax)
    plt.axis([np.log10(nmin),np.log10(nmax),np.log10(Tmin),np.log10(Tmax)])
    plt.xlabel('log n')
    plt.ylabel('log T')
    plt.colorbar(label=r'$\Gamma$ - n$\Lambda$')
    plt.savefig('nT_cmap.png')
    plt.close()


#plot_balance_map()


