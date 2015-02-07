import nupmy as np
from n_T_balance import * #ew
import matplotlib.pyplot as plt
lw = 1.75


# number densities
nmin, nmax = 1.0E-5, 1.0E5
Tmin, Tmax = 1.0  , 1.0E9

n = np.logspace(np.log10(nmin), np.log10(nmax), 1.0E4)

T_equil = find_equilibrium(n,heating_func=heat.metagalactic)
T_equil2 = find_equilibrium(n,cooling_func=cool.IIK_2007,
                              heating_func=heat.IIK_2007, 
                              cf_kwargs={'mode':'jana'})

T_equil3 = find_equilibrium(n,cooling_func=cool.sw_dm,heating_func=heat.metagalactic)

fig = plt.figure(figsize=[18,6])
ax1 = plt.subplot(131)
ax1.plot(n,T_equil ,label="FLASH",lw=lw,color='black', ls='-')
ax1.plot(n,T_equil2,label="IIK_2007",lw=lw,color='black', ls='--')
ax1.plot(n,T_equil3,label="SW_DM", lw=lw, color='black', ls='-.')

ax1.loglog()
ax1.set_xlabel(r'n (cm$^{-3}$)'); ax1.set_ylabel(r'T (K)')
ax1.set_ylim(Tmin,Tmax)
ax1.set_xlim(nmin,nmax)

ax2 = ax1.twinx()
ax2.plot(n, n*T_equil,lw=lw, color='blue',ls='-')
ax2.plot(n, n*T_equil2,lw=lw, color='blue',ls='--')
ax2.plot(n, n*T_equil3, lw=lw, color='blue', ls='-.')

ax2.set_ylabel(r'Pressure (cm$^{-3}$ K)')

for t2 in ax2.get_yticklabels():
    t2.set_color('blue')

ax2.semilogy()

ax1.legend(loc='lower left',fancybox=True)

ax1 = plt.subplot(132)
ax1.plot(T_equil, n*T_equil, label='FLASH',color='black')
ax1.plot(T_equil3,n*T_equil3,label='SWDM',color='red')
ax1.plot(T_equil2,n*T_equil2,label='IIK',color='purple')
ax1.set_xlabel('T (K)')
ax1.set_ylabel('Pressure')
ax1.loglog()
ax1.legend(loc='best')

ax1 = plt.subplot(133)
#ax1.plot(T_equil,  n,label='FLASH',color='black')
ax1.plot(T_equil3, n, label='SWDM',color='red')
ax1.plot(T_equil2, n,label='IIK',color='purple')
ax1.plot(T_equil,n,label='FLASH',color='black')
ax1.set_xlabel('T (K)')
ax1.set_ylabel('n')
ax1.loglog()
ax1.legend(loc='best')

plt.tight_layout()
plt.savefig('nT_balance.png')


plt.close()

