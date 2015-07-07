import numpy as np
import matplotlib.pyplot as plt
import cgs as cgs
import ic_generator as ic


MDM = 2.0E7 * cgs.Msun
sim_rhoAmbient = 2.4E-4 *cgs.mp*cgs.mu
sim_TAmbient   = 1.8E6

sim_TCloud = 1.0E4
sim_rhoCenter = 1.0*cgs.Msun/cgs.pc**3
sim_crho1    = 0.09 * cgs.mp * cgs.mu * 2
sim_rho1rm   = 6.2213294E-27
sim_rho2rm   = sim_rhoAmbient
sim_rhoRL    = sim_rhoAmbient
sim_RL       = 300.0*cgs.pc
sim_bParam   = 600.0*cgs.pc


r = np.linspace(0.1*cgs.pc,1200.0*cgs.pc,1000.0)

RM, dens, P, T= ic.spherical_NFW(r, sim_TCloud, sim_TAmbient, sim_bParam,
                  sim_RL, sim_rhoRL, sim_rhoCenter,
                  sim_rho1rm, sim_rho2rm, sim_crho1,M_DM=MDM)

print RM/cgs.pc

fig, ax1 = plt.subplots()
ax1.plot(r/cgs.pc,dens,label='rho',lw=1.75)
ax1.set_xlabel('r (pc)')
ax1.set_ylabel('rho (cgs)')
ax1.semilogy()

ax2 = ax1.twinx()
ax2.plot(r/cgs.pc,P,label='P',lw=1.75,color='red')
for t1 in ax2.get_yticklabels():
    t1.set_color('red')
ax2.semilogy()

plt.savefig('rho_T.png')
plt.close()
