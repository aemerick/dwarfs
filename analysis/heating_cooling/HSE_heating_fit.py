import numpy as np
from scipy import stats
from yt import units as u

import matplotlib.pyplot as plt
# heating and cooling functions
import cooling as cool
import heating as heat
import dwarf   as dw

# read in some flash.par file
filename = 'flash.par'
order    = 6
out_file = "HSE_heating_fit.dat"

# 
sim = dw.simulation("dwarf_")

r   = np.logspace(14,21,1.0E6) # radii in cm
RM, r, rho, P, T = sim.get_initial_conditions(r)

r = r.value
# compute the heating balance
# gives gamma in erg/s
gamma, n_gamma, nn_delta = heat.heating_balance(rho.value,T.value)
gamma_z                  = heat.metagalactic()
                                                
gamma = gamma + gamma_z

# select over ONLY the dwarf
select = r<0.98*RM
gamma = gamma[select]
r     =     r[select]

# now fit the polynomial !!
pol = np.polyfit(r, gamma, order)
p   = np.poly1d(pol)


sum = 0
i = 1
pol_flip = reversed(pol)
for val in pol_flip:
    sum = sum + val * r**(i-1)
    i = i + 1

# doublecheck with a plot
plt.plot(r, sum  , label='checking implementation', color ='purple')
plt.plot(r, p(r) , label=(str(order) + ' order poly'), color='red')
plt.plot(r, gamma, label='Heating Curve', color='black')
plt.semilogy()
plt.legend(loc='best')
plt.savefig('HSE_fit.png')
plt.close()

#print pol

# now write the output file
f = open(out_file, 'w')
for val in pol_flip:
    print val
    f.write("%14.14e\n"%(val))

f.close()
