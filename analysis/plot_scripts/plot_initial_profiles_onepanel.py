import yt
import dwarf as dw
import numpy as np
import matplotlib.pyplot as plt
import cgs
from initial_conditions import ic_list as icl ; 
from analytic_model import dwarf_model as dw_model
from   matplotlib    import rc

from scipy import integrate

fsize = 18
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 3


simulation_dir = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'

model_names = ['LT_n075_v2_nh3', 'LT_n150_v2_nh3']
models = {}

for name in model_names:
    initial_conditions = icl.ic_object_dict[name]
    models[name] = dw_model.analytical_dwarf(name, initial_conditions.ic)
    
    

r = np.linspace(0.0, initial_conditions.ic['r_HI'], 100)


colors = {'n020' : 'black', 'n075': 'orange', 'n150':'blue'}
labels = {'n020' : r'n$_{\rm{o}}$ = 0.02 cm$^{-3}$',
          'n075' : r'n$_{\rm{o}}$ = 0.75 cm$^{-3}$',
          'n150' : r'n$_{\rm{o}}$ = 1.50 cm$^{-3}$'}

fig, ax1 = plt.subplots()

# now plot the gas profile
for name in model_names:
    
    if 'n020' in name:
        no = 'n020'
    elif 'n075' in name:
        no = 'n075'
    else:
        no = 'n150'
    
    ax1.plot(r/cgs.pc, models[name].gas_profile(r), label = labels[no],
                color = colors[no],
                ls = '-', lw = line_width)
                
                
    Mgas = models[name].ic['M_HI']/cgs.Msun
    power = np.floor(np.log10(Mgas)) 
    Mnumber = Mgas / 10.0**(power)
    
    # would be cool to add in the total gas mass for each of these
    if 'n020' in name:
        x = 50.0 ; y = models[name].gas_profile(x*cgs.pc)
        ax1.annotate(r'M$_{\rm{gas}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(Mnumber, power), 
                     xy=(x, y),
                     xytext=(x, 0.2*y), color=colors[no],
                     arrowprops=dict(facecolor=colors[no],shrink=0.05))
    elif 'n075' in name:
        x = 100.0 ; y = models[name].gas_profile(x*cgs.pc)
        ax1.annotate(r'M$_{\rm{gas}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(Mnumber, power), 
                     xy=(x, y),
                     xytext=(x, 0.25*y), color=colors[no],
                     arrowprops=dict(facecolor=colors[no],shrink=0.05) ) 
    elif 'n150' in name:
        x = 40.0 ; y = models[name].gas_profile(x*cgs.pc)
        ax1.annotate(r'M$_{\rm{gas}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(Mnumber, power), 
                     xy=(x, y),
                     xytext=(x, 2.2*y), color=colors[no],
                     arrowprops=dict(facecolor=colors[no],shrink=0.05) )                                         

        
        
# add burkert profile
#bname = 'Leo_T_Burkert'
#ic = icl.ic_object_dict[bname]
#models[bname] = dw_model.analytical_dwarf(bname, ic.ic)

#plt.plot(r/cgs.pc, models[bname].gas_profile(r), color = 'green', ls = '--', lw=line_width,label='Burkert')

#print 'central density', models[bname].gas_profile(0.01*cgs.pc) / (cgs.mp * ic.ic['mu_dwarf'])

#plt.semilogy()                
#plt.xlim(np.min(r)/cgs.pc, np.max(r)/cgs.pc)
#plt.xlabel(r'Radius (pc)')
#plt.ylabel(r'Gas Density (g cm$^{-3}$)')
#plt.legend(loc='best',fancybox=True)
#plt.savefig('LT_initial_gas_density_burkert.png')
#plt.close()

ax1.set_xlabel(r'Radius (pc)')
ax1.set_ylabel(r'Gas Density (g cm$^{-3}$)')
ax1.semilogy()
ax1.legend(loc='lower left')
r = np.linspace(0.0, initial_conditions.ic['r_DM'], 100)

ax2 = ax1.twinx()

# now plot the dark matter profile
ax2.plot(r/cgs.pc, models[models.keys()[0]].DM_profile(r),
                color = 'black', label= r'$\rho_{\rm{DM,NFW}}$',
                ls = '--', lw = line_width)

#ax2.plot(r/cgs.pc, models[bname].DM_profile(r), color = 'green', ls = '--', lw=line_width,label='Burkert')

ax2.legend(loc='upper right')

ax2.semilogy()                
ax2.set_xlim(np.min(r)/cgs.pc, np.max(r)/cgs.pc)
ax2.set_xlabel(r'Radius (pc)')
ax2.set_ylabel(r'Dark Matter Density (g cm$^{-3}$)')


x = 150.0 ; y = models[name].DM_profile(x*cgs.pc)
ax2.annotate(r'M$_{\rm{DM}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(7.3, 6), 
                     xy=(x, y),
                     xytext=(x, 0.25*y), color='black',
                     arrowprops=dict(facecolor='black',shrink=0.05) )    


ax1.minorticks_on()

fig.set_size_inches(9,8.0)
plt.tight_layout()
fig.savefig('LT_initial_profiles.png')

#.close()

#print "Leo T Burkert Parameters"
#ic.FLASH_readable_ic()



#fx = lambda x : models[bname].DM_profile(x) * 4.0 * np.pi * x * x
#print integrate.quad(fx, 0.0, 300.0 * cgs.pc)[0] / cgs.Msun
#fx = lambda x : models[bname].gas_profile(x) * 4.0 * np.pi * x * x

#print integrate.quad(fx, 0.0, 300.0*cgs.pc)[0] / cgs.Msun
