import numpy as np
from n_T_balance import *

nmin = 1.0E-4
nmax = 1.0
n = np.logspace(np.log10(nmin), np.log10(nmax), 1.0E5)
T= find_equilibrium(n,cooling_func=cool.sw_dm,
                            heating_func=heat.metagalactic)

pressure_eq(n,T,n*T,filename='sw_dm_equilibrium_vals.dat')
