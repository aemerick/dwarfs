import dwarf as danal
import matplotlib.pyplot as plt
import yt
import glob
import numpy as np
import cgs as cgs

# operate on all sorted checkpoint files
ds_names = glob.glob('*_plt_cnt_*')
ds_names.sort() 

m = np.zeros(np.size(ds_names))
t = np.zeros(np.size(ds_names))


# need a better way to do this..
RM = None; i = 0
f = open('mass_t.dat', 'w')
f.write("#t m\n")

sim = danal.simulation('dwarf_fullp_hdf5_plt_cnt_')
RM = 360.3*cgs.pc

for ds_name in ds_names:
    print ds_name
    ds = yt.load(ds_name)

    dw = danal.dwarf(ds, 'flash.par', rm=RM)
  #  if RM == None:
   #     RM = dw.rmatch(1000)
    



    m[i] = dw.param_contained('total_mass', r = RM).convert_to_units('Msun')
    t[i] = (dw.time).convert_to_units('Myr')
    f.write("%8.8e %8.8e\n"%(t[i],m[i]))
    
    i = i + 1
f.close()
plt.plot(t,m/m[0],color='black',lw=1.5, label=r"M$_o$=%3.2E M$_{\odot}$"%(m[0]))
plt.xlabel('t (Myr)')
plt.ylabel(r'M(t)/M$_o$')
plt.legend(loc='best')
#plt.semilogy()
plt.ylim(0.0,1.0)
plt.xlim(t[0],t[-1])
plt.savefig('plots/mass_t.png')
plt.close()    
    
#np.savetxt('mass_t.dat', (t,m), fmt='%6.5e',header="#t m")
    
