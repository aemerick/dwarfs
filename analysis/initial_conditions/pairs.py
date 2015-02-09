import numpy as np
from n_T_balance import *
from scipy import interpolate as interp
from heating_cooling import cooling as cool
from heating_cooling import heating as heat

nmin = 1.0E-6
nmax = 100.0
n = np.logspace(np.log10(nmin), np.log10(nmax), 5.0E4)
#T = find_equilibrium(n,cooling_func=cool.IIK_2007,
#                       heating_func=heat.IIK_2007)
T= find_equilibrium(n,cooling_func=cool.sw_dm, heating_func=heat.metagalactic)

T = T[np.logical_not(np.isnan(T))]
n = n[np.logical_not(np.isnan(T))]
n = n[np.argsort(T)]
T = T[np.argsort(T)]
#print "--"
#print np.min(T), np.max(T)
#print T[0], T[1], T[-2], T[-1]
#function = interp.UnivariateSpline(T,n*T)

# find the pressure equillibrium for a few values of the corona
P_spline = interp.UnivariateSpline(T,n*T)
n_spline = interp.UnivariateSpline(T,n)

Tsample = np.logspace(np.log10(np.min(T)),np.log10(np.max(T)),1.0E4)

pressure_eq(n_spline(Tsample),Tsample,P_spline(Tsample),filename='sw_dm_equilibrium_vals.dat')

plt.scatter(T,n*T)
#plt.plot(n_spline(Tsample),n_spline(Tsample)*Tsample)
#Tfun = np.logspace(np.log10(np.min(T)),np.log10(np.max(T)),1000)
#vals = function(Tfun)

#plt.plot(Tfun,vals,
#           label='spline',lw=1.0,color='red')
plt.loglog()
#plt.xlim(1.0E2,1.0E8)
plt.savefig('t.png')
plt.close()
