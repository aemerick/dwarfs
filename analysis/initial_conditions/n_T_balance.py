import numpy as np
from heating_cooling import cooling as cool
from heating_cooling import heating as heat
from scipy import optimize as opt
from plotting import plotTools as myplot
import matplotlib.pyplot as plt
lw = 1.75
def find_equilibrium(n, Tmin=1.E2, Tmax=1.E7, tolerance=1.0E-12,
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
        heating_func = lambda T: heat.lower_bound_heating(0.0,T)
    
    T_equil = np.zeros(np.size(n))
    
    for i in np.arange(np.size(n)):
        function = lambda T: -heating_func(T) + n[i]*cooling_func(T)

        if(np.sign(function(Tmin)) == np.sign(function(Tmax))):
            root = None
        else:
            root = opt.bisect(function, Tmin, Tmax, xtol=tolerance)
        
        T_equil[i] = root
        
    return T_equil
    

def plot_balance_map(nmin=1.0E-5,nmax=1.0E3,Tmin=10.0,Tmax=1.0E7,npoints=1.0E3):

    n = np.linspace(np.log10(nmin),np.log10(nmax),npoints)
    T = np.linspace(np.log10(Tmin),np.log10(Tmax),npoints)



    func = lambda n,T: heat.lower_bound_heating(0.0,T) - n*cool.radloss(T)
    
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


plot_balance_map()

# number densities
nmin, nmax = 1.0E-5, 1.0E3
Tmin, Tmax = 10.0  , 1.0E7

n = np.logspace(np.log10(nmin), np.log10(nmax), 1.0E3)

T_equil = find_equilibrium(n)
T_equil2 = find_equilibrium(n,cooling_func=cool.IIK_2007,heating_func=heat.IIK_2007)
fig,ax1 = plt.subplots()

ax1.plot(n,T_equil,label="FLASH",lw=lw,color='black', ls='-')
ax1.plot(n,T_equil2,label="IIK_2007",lw=lw,color='black', ls='--')
ax1.loglog()
ax1.set_xlabel(r'n (cm$^{-3}$)'); ax1.set_ylabel(r'T (K)')
ax1.set_ylim(Tmin,Tmax)
ax1.set_xlim(nmin,nmax)

ax2 = ax1.twinx()
ax2.plot(n, n*T_equil,lw=lw, color='blue',ls='-')
ax2.plot(n, n*T_equil2,lw=lw, color='blue',ls='--')
ax2.set_ylabel(r'Pressure (cm$^{-3}$ K)')

for t2 in ax2.get_yticklabels():
    t2.set_color('blue')

ax2.semilogy()

ax1.legend(loc='upper left',fancybox=True)
plt.savefig('nT_balance.png')
plt.close()


