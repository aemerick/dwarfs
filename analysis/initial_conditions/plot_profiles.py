import ic_list as icl
import numpy as np
import cgs as cgs
import matplotlib.pyplot as plt

leo_T = icl.dwarf_dict['Leo_test']

r = np.linspace(0.0,2000.0,10000)*cgs.pc
leo_T.find_density_profile(r,type='NFW_isothermal_rmatch')

fig, ax1 = plt.subplots()

ax1.plot(leo_T.rvals/cgs.pc, leo_T.rho, lw=1.75,color='black')
ax1.plot([leo_T.radius/cgs.pc,leo_T.radius/cgs.pc],plt.ylim(),lw=1.75,color='black',ls='--')
ax1.semilogy()
ax1.set_xlabel(r'r (pc)')
ax1.set_ylabel(r'$\rho$ (g cm$^{-3}$)')

ax2 = ax1.twinx()

Pdwarf = leo_T.rho / (leo_T.ic['mu_dwarf']*cgs.mp) * cgs.kb * leo_T.ic['T_dwarf']
Pcorona = leo_T.ic['n_halo'] * cgs.kb * leo_T.ic['T_halo']
ax2.plot(leo_T.rvals/cgs.pc, Pdwarf,color='red')
ax2.plot(plt.xlim(), [Pcorona,Pcorona], color='red',ls='--')
ax2.semilogy()
ax2.set_ylabel(r'Pressure')
plt.savefig('leo_t_test.png')
plt.close()


# gas mass 
r = r[r<=leo_T.radius]
rho = leo_T.rho[r<=leo_T.radius]
dr = r[1:] - r[0:-1]
M = 4.0*np.pi*np.trapz(r*r*rho, dx = dr)/3.0
print "within radius M = %5.4e -- Msun = %5.4e"%(M, M/cgs.Msun)
print "average n = %5.4e"%(np.average(rho)/(cgs.mp*leo_T.ic['mu_dwarf']))

r = r[r<=300.0*cgs.pc]
rho = rho[r<=300.0*cgs.pc]
dr = r[1:] - r[0:-1]
M = 4.0*np.pi*np.trapz(r*r*rho, dx = dr)/3.0
print "within 300 pc M = %5.4e -- Msun = %5.4e"%(M, M/cgs.Msun)
print "average n = %5.4e"%(np.average(rho)/(cgs.mp*leo_T.ic['mu_dwarf']))

