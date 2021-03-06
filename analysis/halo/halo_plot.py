import numpy as np
import matplotlib.pyplot as plt
import galaxy as gal # galaxy functions
import cgs    as cgs # cgs unit conversions

lw = 1.7
# plot density and temperature profiles for 
# things
r = np.linspace(1.0, 450.0, 1000.0) * cgs.kpc

n = {'miller': gal.halo_gas_ndensity(r),
     'miller_adj': gal.halo_gas_ndensity(r,n0=(0.46+0.74),rc=(0.35+0.29)*cgs.kpc,
                                           beta=(0.71-0.14),ambient=0.0),
     'isothermal': gal.gato_ndensity(r, 'isothermal'),
     'adiabatic' : gal.gato_ndensity(r, 'adiabatic'),
     'cooling'   : gal.gato_ndensity(r, 'cooling'),
     'kaufmann01' : gal.kaufmann('01', 'density'),
     'kaufmann02' : gal.kaufmann('02', 'density')}

T = {'NFW' : gal.halo_gas_temperature(r),
     'miller': gal.halo_gas_temperature(r,n=n['miller']),
     'isothermal': gal.halo_gas_temperature(r,n=n['isothermal']),
     'adiabatic': gal.halo_gas_temperature(r,n=n['adiabatic']),
     'cooling': gal.halo_gas_temperature(r,n=n['cooling']),
     'kaufmann01' : gal.kaufmann('01','temperature'),
     'kaufmann02' : gal.kaufmann('02','temperature')}

r = r / cgs.kpc
# number density
fig = plt.figure(figsize=[12,6])

ax1 = fig.add_subplot(121)
ax1.plot(r, n['miller'], label = 'Miller and Bregman 2013', color ='black',lw=lw)
ax1.plot(r, n['miller_adj'], label = 'Adjusted MB2013', color ='black',lw=lw,ls='--')
ax1.plot(r, n['isothermal'], label = "Gato Isothermal",
                                              lw = lw, color = 'blue',
                                              ls = '-')
ax1.plot(r, n['adiabatic'], label = "Gato Adiabatic",
                                              lw = lw, color = 'blue',
                                              ls = '--')
ax1.plot(r, n['cooling'], label = "Gato Cooling",
                                              lw = lw,color='blue',
                                              ls = '-.')

ax1.plot(gal.kaufmann('01','r'), n['kaufmann01'], label = 'Kaufmann NFW',
                                           color = 'purple', lw = lw)
ax1.plot(gal.kaufmann('02','r'), n['kaufmann02'], label = 'Kaufmann Realistic',
                                           color = 'purple', lw=lw, ls='--')


# minimum n out to 70 kpc from Grcevich & Putman 2009
ax1.plot([1.0,70.0],[2.0E-4,2.0E-4], label='Grcevich & Putman 2009',
                                color='green', lw=lw, ls='-')
ax1.plot([1.0,70.0],[3.0E-4,3.0E-4],color='green', lw=lw, ls='-')

ax1.scatter([425,425,425],[4.6E-5,1.5E-4,4.5E-4])
                                           
ax1.semilogy()
ax1.legend(loc='best')
ax1.set_xlabel(r'r (kpc)')
ax1.set_ylabel(r'n (cm$^{-3}$)')
ax1.set_xlim(1.0,np.max(r))                                              
ax2 = fig.add_subplot(122)

ax2.scatter([425,425,425],[7.4E5,7.5E5,8.1E5])
#ax2.plot(r, T['miller'], label = 'Miller and Bregman 2013', color ='black',lw=lw)

#ax2.plot(r, T['isothermal'], label = "Gato Isothermal",
#                                              lw = lw, color = 'blue',
#                                              ls = '-')
#ax2.plot(r, T['adiabatic'], label = "Gato Adiabatic",
#                                              lw = lw, color = 'blue',
#                                              ls = '--')
#ax2.plot(r, T['cooling'], label = "Gato Cooling",
#                                              lw = lw,color='blue',
#                                              ls = '-.')
ax2.plot(r, T['NFW'], label = 'HSE with NFW DM', color = 'red', lw = lw)
ax2.plot(gal.kaufmann('01','r'), T['kaufmann01'], label = 'Kaufmann NFW',
                                           color = 'purple', lw = lw)
ax2.plot(gal.kaufmann('02','r'), T['kaufmann02'], label = 'Kaufmann Realistic',
                                           color = 'purple', lw=lw, ls='--')
ax2.semilogy()
ax2.set_xlabel(r'r (kpc)')
ax2.set_ylabel(r'T (K)')
#ax2.legend(loc='best')
ax2.set_xlim(1.0,np.max(r))                                              

ax2.plot([50.0,90.0],[1.8E6,1.8E6],label=r'Gato T$_{min}$ - T$_{max}$',color='blue',
                                  lw=lw, ls='-')
ax2.plot([50.0,90.0],[3.0E6,3.0E6],color='blue',
                                  lw=lw, ls='-')
ax2.legend(loc='best')
plt.savefig('n_T_profile.png')
plt.close()                                              
                                             
