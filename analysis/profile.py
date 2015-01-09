import dwarf as danal
import matplotlib.pyplot as plt
import yt
import glob
import numpy as np

# get all checkpoint files as a list
ds_names = glob.glob('*_chk_*')
ds_names.sort()

RM = None; i = 0
for ds_name in ds_names:

    print ds_name
    ds = yt.load(ds_name)
    dw = danal.dwarf(ds,'flash.par',rm=RM)
    if RM == None:
        RM = dw.rmatch(1000)

    sphere = dw.ds.sphere(dw.center, (2.0*RM.convert_to_units('cm'),'cm'))

    bins, field = dw.profile('density', nbin = 20, weight = 'cell_mass',
                  xfield='radius',xmin=0.0*yt.units.cm,xmax=2.0*RM,data_source=sphere)

    field = np.append(field, field[-1])

    bins = bins.convert_to_units('pc').value
    

    plt.step(bins, field.value, label='Density', lw = 1.5, where='post')

    plt.semilogy()

    plt.xlabel('R (pc)')
    plt.ylabel('Density (cgs)')

    RMpc = RM.convert_to_units('pc').value
    plt.ylim(1.0E-28,5.0E-23)
    plt.xlim(0,1.5*RMpc)
    plt.plot([RMpc,RMpc],plt.ylim(), color='black', ls='-.', lw=1.5)

    plt.savefig('plots/density/%0004i'%(i))
    plt.close()
    i = i + 1
    
