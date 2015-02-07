import numpy as np
from heating_cooling import cooling as cool
from heating_cooling import heating as heat
from scipy import optimize as opt
from scipy import interpolate as interp
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
def pressure_eq(n,T,P, Tcorona = [5.0E5,5.0E6], Tdwarf=[1.0E4,3.0E4],
                       filename='equillibrium_nT.dat'):
    """
       Given n, T, and P, compute some equillibrium points. This operates
       by searching over each corona temperature value for the corresponding
       equilibrium value in the dwarf range of temperatues by using a 
       cubic spline fit the P(T) and n(T) for the dwarf values.

        Parameters
        ----------
        n : array
            Array of number densities for the cooling curve equilibrium.
        T : array
            Array of temperatures for the cooling curve equilibrium.
        P : array
            Array of pressures. Usually taken as n * T
        Tcorona: list, optional
            2d list of possible corona tempertures [min,max]. Default
            [5.0E5, 6.0E6]
        Tdwarf : list, optional 
            List of possible dwarf temperatures to search for
            pressure eq. with the corona [min,max]. Default [1.0E3,8.0E4]
        filename : string, optional
            File to write out parameters to. Default 'equilibrium_NT.dat'

        Returns
        -------
        dwarf_eq : dict
            Dictionary of dwarf equilibrium values, with n, T, P arrays
            given with the keys 'n', 'T', 'P'.
        corona_eq : dict
            Same as above, but for the corona values
    """

    # Prune out any nans if they do exist
    T = T[np.logical_not(np.isnan(T))]
    P = P[np.logical_not(np.isnan(T))]
    n = n[np.logical_not(np.isnan(T))]

    # since things will be functions of T, sort
    # them by T:
    sort = np.argsort(T)
    T = T[sort] ; P = P[sort] ; n = n[sort]

    print np.size(T), np.size(n), np.size(P)
    print np.min(T), np.max(T)

    T_cmin, T_cmax = Tcorona[0],Tcorona[1]
    Td_min, Td_max = Tdwarf[0],Tdwarf[1]
    print np.size(T[T>=T_cmin]), T_cmin, T[T>=T_cmin]
    print np.size(T[T<=T_cmax]), T_cmax, T[T<=T_cmax]
    # choose some possible corona temperatures from the evaluated T's
    T_corona = T[(T>=T_cmin)*(T<=T_cmax)]        
    P_corona = P[(T>=T_cmin)*(T<=T_cmax)]
    n_corona = n[(T>=T_cmin)*(T<=T_cmax)]   

    print 'num corona', np.size(T_corona)
    # Get the range of dwarf paramers for the temperature range.
    # use this to get them into functions of temperature
    dwarf_temperatures = T[(T>=Td_min)*(T<=Td_max)]
    dwarf_pressures =    P[(T>=Td_min)*(T<=Td_max)]
    dwarf_densities =    n[(T>=Td_min)*(T<=Td_max)]
    dwarf_Pmin, dwarf_Pmax = np.min(dwarf_pressures), np.max(dwarf_pressures)

    print 'num dwarf', np.size(dwarf_temperatures)
    # interpolate the P(T) and n(T) curves for the dwarf temperatures
    P_spline = interp.UnivariateSpline(dwarf_temperatures, dwarf_pressures)
    n_spline = interp.UnivariateSpline(dwarf_temperatures, dwarf_densities)    

    plt.plot(dwarf_temperatures,P_spline(dwarf_temperatures))
    plt.plot(dwarf_temperatures,n_spline(dwarf_temperatures))
    plt.loglog()
    plt.savefig('spline.png');    plt.close()

    # dwarf_equillibrium arrays
    dwarf_eq = {'n': np.zeros(np.size(P_corona)),
                'T': np.zeros(np.size(P_corona)),
                'P': np.zeros(np.size(P_corona))}

    # find all of the matched temperatues for 
    for i in np.arange(np.size(P_corona)):
        func = lambda T: P_spline(T) - P_corona[i]

        # find the actual root and save it
        if np.sign(func(Td_min)) == np.sign(func(Td_max)):
            #dwarf_eq['T'][i] = None
            #dwarf_eq['n'][i] = None
            #dwarf_eq['P'][i] = None
            print "Tdmin and Tdmax are both positive"
            print Td_min,Td_max,func(Td_min),func(Td_max)  
        else:
            root = opt.bisect(func, Td_min, Td_max)
     
        
            # save the equillibrium values
            dwarf_eq['T'][i] = root
            dwarf_eq['n'][i] = n_spline(root)
            dwarf_eq['P'][i] = P_spline(root)
        

    # now output everything to a file
    file = open(filename, 'w')
    
    file.write("# dwarf_n dwarf_T dwarf_P halo_n halo_T halo_P\n")
    fmt = "%6.6e %6.6e %6.6e %6.6e %6.6e %6.6e\n"
    for i in np.arange(np.size(P_corona)):
        
        file.write(fmt%(dwarf_eq['n'][i],dwarf_eq['T'][i],dwarf_eq['P'][i],
                        n_corona[i], T_corona[i], P_corona[i]))
 
    file.close()

    corona_eq = {'n': n_corona, 'T': T_corona, 'P':P_corona}

    return dwarf_eq, corona_eq
