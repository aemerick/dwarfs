import matplotlib.pyplot as plt
import numpy as np
import cgs as cgs
from matplotlib import rc
import dwarf as dw


from astropy import units as u

fsize = 17
rc('text', usetex=True)
rc('font', size=fsize)
lw = 2.9
point_size = 30

# 
lvl3_ls = ':'
lvl4_ls = '-.'
lvl5_ls = '-'
lvl6_ls = '--'


#####
fpath = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'
fpathv2 = fpath + 'SN_LT_n150_v2_nh4/resolution_study/'
fpathv4 = fpath + 'SN_LT_n150_v4_nh4/resolution_study/'
filename = 'grav_bound_mass.dat'
outname = 'LT_n150_nh4_resolution.png'

fig, ax = plt.subplots(1,2)
def print_results(name, data):
    strformat='%50s %8.8f %8.8f'
    
    output_string  = name.replace(fpath,'').replace(filename,'')
    stripping_time =  dw.predict_stripping_time(data['t'], data['m'])
    
    gas_fraction = data['m'][ np.argmin(np.abs(data['t'] - 2000.0)) ] /data['m'][0]
    
    
    print strformat%(output_string, stripping_time, gas_fraction)
    
name = 'lvl3/'

fname = fpathv2 + name + filename
data = np.genfromtxt(fname, names=True)
ax[0].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
          ls = lvl3_ls, label = '39.06 pc')
print_results(name, data)

fname = fpathv4 + name + filename
data = np.genfromtxt(fname, names=True)
ax[1].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
          ls = lvl3_ls, label = '39.06 pc')
print_results(name, data)

# ---------------------------------

name = 'lvl4/'
fname = fpathv2 + name + filename

data = np.genfromtxt(fname, names=True)
ax[0].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
            ls = lvl4_ls, label = '19.54 pc')
print_results(name, data)

fname = fpathv4 + name + filename

data = np.genfromtxt(fname, names=True)
ax[1].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
            ls = lvl4_ls, label = '19.54 pc')
print_results(name, data)



name = 'lvl5/'
fname = fpathv2 + name + filename

data = np.genfromtxt(fname, names=True)
ax[0].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
            ls = lvl5_ls, label = '9.77 pc')
print_results(name, data)

fname = fpathv4 + name + filename

data = np.genfromtxt(fname, names=True)
ax[1].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
            ls = lvl5_ls, label = '9.77 pc')
print_results(name,data)


name = 'lvl6/'
fname = fpathv2 + name + filename

data = np.genfromtxt(fname, names=True)
ax[0].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
            ls = lvl6_ls, label = '4.89 pc')
print 'lvl6 v2', data['t'][-1]
print_results(name, data)

fname = fpathv4 + name + filename
data = np.genfromtxt(fname, names=True)
ax[1].plot(data['t'], data['m']/data['m'][0], lw=lw, color = 'black',
            ls = lvl6_ls, label = '4.89 pc')
print 'lvl6 v4', data['t'][-1]
print_results(name, data)


for a in ax:
    a.set_xlabel(r'Time (Myr)')
    a.set_ylabel(r'M(t) / M$_{\rm{o}}$')
    a.set_xlim(0.0,2000.0)

    a.set_ylim(0.0,1.2)

    a.minorticks_on()

    a.legend(loc='best')

fig.set_size_inches(12,6)
plt.tight_layout()
fig.savefig(outname)

