import numpy as np
from scipy import stats
from yt import units as u
import matplotlib.pyplot as plt
import cooling as cool
import heating as heat
import dwarf   as dw

#
lw = 1.75

#

# T array
T = np.logspace(1.0,9.0,1.0E4)
temperature_label = 'log T (k)'


# cooling
delta = cool.radloss(T)

plt.plot(T, delta, lw = lw, label = "Dalgarno and McCray (1972)",
              color = "black")
plt.loglog()
plt.ylim(1.0E-29, 1.0E-21)

plt.ylabel(r'log $\Lambda$(T) (erg cm$^{3}$ s$^{-1}$)')
plt.xlabel(temperature_label)
plt.legend(loc='best')
plt.savefig('cooling_curve.png')
plt.close()
#########################################

# load initial conditions for heating balance
sim = dw.simulation('dwarf_')

r = np.logspace(14,21,1.0E5) # radii in cm from dwarf center
RM, r, rho, P, T = sim.get_initial_conditions(r)

r = r.value

gamma, n_gamma, nn_delta = heat.heating_balance(rho.value, T.value)
gamma_z                  = heat.metagalactic()

plt.plot(T.value, gamma, label = "HSE Heating rate", lw = lw, color = 'black')
plt.plot(plt.xlim(),[gamma_z,gamma_z], label = "HM12, z=0", lw = lw, color = 'black', ls='--')
plt.xlabel(temperature_label)
plt.ylabel(r'$\Gamma$ (erg s$^{-1}$)')
plt.loglog()
plt.legend(loc='best')
plt.savefig('heating_curve.png')
plt.close()


r = (r*u.cm).convert_to_units('pc').value


plt.plot(r, n_gamma, label = 'HSE Heating - %5.4e'%(np.min(n_gamma)), lw = lw, color = 'black')
plt.plot(r, nn_delta, label = 'Cooling', lw = lw, color = 'red')
plt.xlabel('r (pc)')
plt.ylabel(r'Heating - Cooling (erg cm$^{-3}$ s$^{-1}$)')
plt.semilogy()
plt.legend(loc='best')
plt.savefig('radial_curves.png')
plt.close()

plt.plot(r,gamma,label='HSE Heating',lw=lw,color='black')
plt.plot(plt.xlim(),[gamma_z,gamma_z],label="HM12, z=0",lw=lw,color='black',ls='--')
plt.xlabel('r (pc)')
plt.ylabel(r'$\Gamma$ (erg s$^{-1}$)')
plt.semilogy()
#plt.legend(loc='best')
#plt.savefig('radial_heating.png')
#plt.close()


#### fit line to radial heating curve ####
slope, intercept, rval,pval,stderr = stats.linregress(r[r<295],np.log10(gamma[r<295]))
print "slope and intercept of linear function"
print "%5.4e %5.4e"%(slope,intercept)
plt.plot(r,10**(slope*r + intercept), label='linear',color='red')
# fit exponential
r = r[r<295]
gamma = gamma[r<295]

slope, intercept, rval,pval,stderr = stats.linregress(r, np.log(gamma))
print "now for the exponential (y = Aexp(Bx) "
print "B, 1/B, A"
print "%5.4e %5.4e %5.4e"%(slope,1/slope,np.exp(intercept))
plt.plot(r,np.exp(intercept)*np.exp(slope*r),label='exponential',color='green')



pol = np.polyfit(r,gamma,8)
p = np.poly1d(pol)
#i = 4
#y = np.zeros(np.size(r))
#for p in pol:
#    y = y + r**(i) * p
#    i = i - 1

print "polynomial", pol
plt.plot(r,p(r),label='n=4 polynomial',color='purple')


plt.legend(loc='best')
plt.savefig('radial_heating.png')
plt.close()


## 
#print "trying n = 4 polynomial fit"
#$ol = np.polyfit(r, gamma)
#print pol
#print 10**pol[0], 10**(
