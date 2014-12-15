import dwarf as danal
import matplotlib.pyplot as plt
import yt
import glob
import numpy as np

# operate on all sorted checkpoint files
ds_names = glob.glob('*_chk_*')
ds_names.sort() 
#ds_names = [ds_names[0]] + [ds_names[10]] + [ds_names[100]] + [ds_names[1000]]
m = np.zeros(np.size(ds_names))
t = np.zeros(np.size(ds_names))


# need a better way to do this..
RM = None; i = 0
for ds_name in ds_names:
    print ds_name
    ds = yt.load(ds_name)

    dw = danal.dwarf(ds, 'flash.par', rm=RM)
    if RM == None:
        RM = dw.rmatch(1000)
    



    m[i] = dw.param_contained('total_mass', r = RM).convert_to_units('Msun')
    t[i] = (dw.time).convert_to_units('Myr')
    
    i = i + 1
plt.plot(t,m/m[0],color='black',lw=1.5, label=r"M$_o$=%3.2E M$_{\odot}$"%(m[0]))
plt.xlabel('t (Myr)')
plt.ylabel(r'M(t)/M$_o$')
plt.legend(loc='best')
#plt.semilogy()
plt.ylim(0.0,1.0)
plt.savefig('plots/mass_t.png')
plt.close()    
    

    
