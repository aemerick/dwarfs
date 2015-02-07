import numpy as np
import cooling as cool
import matplotlib.pyplot as plt
lw = 1.75
T = np.logspace(1.0,9.0,10000.0)

labels = {'DM': 'Dalgarno and McCray 1972',
          'SW': 'Sarazin and White 1987',
          'SWDM': 'SW + DM',
          'IIK' : 'IIK 2007'}
functions = {'DM': cool.radloss,
             'SW': cool.sarazin_white,
             'SWDM': cool.sw_dm,
             'IIK' : cool.IIK_2007}
colors = {'DM': 'black',
          'SW': 'blue',
          'SWDM': 'red',
          'IIK' : 'purple'}


for name in labels:
    L = functions[name](np.zeros(np.size(T)),T)
    
    plt.plot(T,L, label=labels[name], color=colors[name], lw=lw,ls='-')


plt.xlabel(r'T (K)')
plt.ylabel(r'$\Lambda$ (erg cm$^{3}$ s$^{-1}$)')
plt.loglog()
plt.legend(loc='best',fancybox=True)
plt.savefig('cooling_curves.png')
plt.close()
