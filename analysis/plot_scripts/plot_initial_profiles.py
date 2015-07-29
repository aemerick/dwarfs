import yt
import dwarf as dw
import numpy as np
import matplotlib.pyplot as plt
import cgs
from initial_conditions import ic_list as icl ; 
from analytic_model import dwarf_model as dw_model
from   matplotlib    import rc

fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.5


simulation_dir = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'

model_names = ['LT_n020_v2_nh3', 'LT_n075_v2_nh3', 'LT_n150_v2_nh3']
models = {}

for name in model_names:
    initial_conditions = icl.ic_object_dict[name]
    models[name] = dw_model.analytical_dwarf(name, initial_conditions.ic)
    
    

r = np.linspace(0.0, initial_conditions.ic['r_HI'], 100)


colors = {'n020' : 'black', 'n075': 'red', 'n150':'blue'}
labels = {'n020' : r'n$_{\rm{o}}$ = 0.02 cm$^{-3}$',
          'n075' : r'n$_{\rm{o}}$ = 0.75 cm$^{-3}$',
          'n150' : r'n$_{\rm{o}}$ = 1.50 cm$^{-3}$'}

# now plot the gas profile
for name in model_names:
    
    if 'n020' in name:
        no = 'n020'
    elif 'n075' in name:
        no = 'n075'
    else:
        no = 'n150'
    
    plt.plot(r/cgs.pc, models[name].gas_profile(r), label = labels[no],
                color = colors[no],
                ls = '-', lw = line_width)
                
                
    Mgas = models[name].ic['M_HI']/cgs.Msun
    power = np.floor(np.log10(Mgas)) 
    Mnumber = Mgas / 10.0**(power)
    
    # would be cool to add in the total gas mass for each of these
    if 'n020' in name:
        x = 50.0 ; y = models[name].gas_profile(x*cgs.pc)
        plt.annotate(r'M$_{\rm{gas}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(Mnumber, power), 
                     xy=(x, y),
                     xytext=(x, 0.2*y), color=colors[no],
                     arrowprops=dict(facecolor=colors[no],shrink=0.05))
    elif 'n075' in name:
        x = 150.0 ; y = models[name].gas_profile(x*cgs.pc)
        plt.annotate(r'M$_{\rm{gas}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(Mnumber, power), 
                     xy=(x, y),
                     xytext=(x, 0.275*y), color=colors[no],
                     arrowprops=dict(facecolor=colors[no],shrink=0.05) ) 
    elif 'n150' in name:
        x = 175.0 ; y = models[name].gas_profile(x*cgs.pc)
        plt.annotate(r'M$_{\rm{gas}}$ = %.2f $\times$ 10$^{%1i}$ M$_{\odot}$'%(Mnumber, power), 
                     xy=(x, y),
                     xytext=(x, 2.2*y), color=colors[no],
                     arrowprops=dict(facecolor=colors[no],shrink=0.05) )                                         

plt.semilogy()                
plt.xlim(np.min(r)/cgs.pc, np.max(r)/cgs.pc)
plt.xlabel(r'Radius (pc)')
plt.ylabel(r'Gas Density (g cm$^{-3}$)')
plt.legend(loc='best',fancybox=True)
plt.savefig('LT_initial_gas_density.png')
plt.close()

r = np.linspace(0.0, initial_conditions.ic['r_DM'], 100)

# now plot the dark matter profile
plt.plot(r/cgs.pc, models[models.keys()[0]].DM_profile(r),
                color = 'black',
                ls = '-', lw = line_width)

plt.semilogy()                
plt.xlim(np.min(r)/cgs.pc, np.max(r)/cgs.pc)
plt.xlabel(r'Radius (pc)')
plt.ylabel(r'Dark Matter Density (g cm$^{-3}$)')

plt.savefig('LT_initial_DM_density.png')
plt.close()

