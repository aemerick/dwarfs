import numpy as np
import matplotlib.pyplot as plt
import dwarf as dw

from matplotlib import rc

fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
lw     = 3
point_size = 30


fpath = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'
file_name = 'grav_bound_mass.dat'

fig, ax = plt.subplots(2,2)


####
nsn_ls = '--'
sn_ls  = '-'

###### 
def print_results(name, data):
    strformat='%50s %8.8f %8.8f'
    
    output_string  = name.replace(fpath,'').replace(file_name,'')
    stripping_time =  dw.predict_stripping_time(data['t'], data['m'])
    
    gas_fraction = data['m'][ np.argmin(np.abs(data['t'] - 2000.0)) ] /data['m'][0]
    
    
    print strformat%(output_string, stripping_time, gas_fraction)

    
def compute_dMdt(data):
    
    dMdt = np.diff(data['m']) / np.diff(data['t'])
    
    
    return dMdt

def bin_data(x, y, xmin, xmax, bins= None, dx = 500):
    
    if bins == None:
        nbins = (xmax - xmin) / dx + 1
    
        bins = np.linspace(xmin, xmax, nbins)
        
    

    binned_data = np.zeros(np.size(bins))
    
    for i in np.arange(np.size(bins) - 1):

        binned_data[i] = np.average(y[ (x>=bins[i])*(x<bins[i+1])])
    
    binned_data[-1] = binned_data[-2]
    return bins, binned_data

sign = -1

####### SN_LT_n075_v2_nh4 #########
bins = np.linspace(0.0, 2000.0, 11)

name  = 'SN_LT_n075_v2_nh4/'
fname = fpath + name + '0000_cfloor_nsn/' + file_name

data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)


ax[(0,0)].step(binned_t,sign* binned_dMdt, lw=lw, color='black', where='post',
                          ls = nsn_ls, label = 'No SN')

fname = fpath + name + '0000_cfloor_global/' + file_name
data = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)

ax[(0,0)].step(binned_t,sign* binned_dMdt, lw=lw, color='black', where = 'post',
                          ls = sn_ls, label = 'SN')




####### SN_LT_n075_v4_nh4 #########
name  = 'SN_LT_n075_v4_nh4/'

fname = fpath + name + '0000_cfloor_nsn/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(0,1)].step(binned_t,sign* binned_dMdt, lw=lw, color='black', where='post',
                          ls = nsn_ls, label = 'No SN')


fname = fpath + name + '0000_cfloor_global/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(0,1)].step(binned_t,sign* binned_dMdt, lw=lw, color ='black', where='post',
                          ls = sn_ls, label = 'SN')

fname = fpath + name + '0000_cfloor_global_x2/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(0,1)].step(binned_t,sign* binned_dMdt, lw=lw, color ='red', where='post',
                          ls = ':', label = 'SN x 2')

####### SN_LT_n150_v2_nh4 #########
name  = 'SN_LT_n150_v2_nh4/'

fname = fpath + name + '0000_cfloor_nsn/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(1,0)].step(binned_t,sign* binned_dMdt, lw=lw, color='black', where='post',
                          ls = nsn_ls, label = 'No SN')


fname = fpath + name + '0000_cfloor_global/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(1,0)].step(binned_t,sign* binned_dMdt, lw=lw, color ='black', where='post',
                          ls = sn_ls, label = 'SN')

####### SN_LT_n150_v4_nh4 #########
name  = 'SN_LT_n150_v4_nh4/'

fname = fpath + name + '0000_cfloor_nsn/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(1,1)].step(binned_t,sign* binned_dMdt, lw=lw, color='black', where='post',
                          ls = nsn_ls)#, label = 'No SN')


fname = fpath + name + '0000_cfloor_global/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(1,1)].step(binned_t,sign* binned_dMdt, lw=lw, color ='black', where='post',
                          ls = sn_ls)#, label = 'SN')

fname = fpath + name + '0000_cfloor_global_x2/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(1,1)].step(binned_t,sign* binned_dMdt, lw=lw, color ='red', where='post',
                          ls = ':', label = 'SN x 2')

fname = fpath + name + '0000_cfloor_global_x5/' + file_name
data  = np.genfromtxt(fname, names=True)
dMdt = compute_dMdt(data)

binned_t, binned_dMdt = bin_data(data['t'][:-1], dMdt, 0.0, 2000.0, bins = bins)
ax[(1,1)].step(binned_t,sign* binned_dMdt, lw=lw, color ='green', where='post',
                          ls = ':', label = 'SN x 5')


################################################################
# annotations #
xpos = 1120.0
ypos = 0.95

xy = (xpos, ypos)
ax[(0,0)].annotate(r'v = 200.0 km/s',xy=xy,xytext=xy,textcoords='data')
xy = (xpos, ypos-0.085)
ax[(0,0)].annotate(r'n$_{\rm{o}}$ = 0.75 cm$^{-3}$',xy=xy,xytext=xy,textcoords='data')


xy = (xpos, ypos)
ax[(0,1)].annotate(r'v = 400.0 km/s',xy=xy,xytext=xy,textcoords='data')
xy = (xpos, ypos-0.085)
ax[(0,1)].annotate(r'n$_{\rm{o}}$ = 0.75 cm$^{-3}$',xy=xy,xytext=xy,textcoords='data')

xy = (xpos, ypos)
ax[(1,0)].annotate(r'v = 200.0 km/s',xy=xy,xytext=xy,textcoords='data')
xy = (xpos, ypos-0.085)
ax[(1,0)].annotate(r'n$_{\rm{o}}$ = 1.50 cm$^{-3}$',xy=xy,xytext=xy,textcoords='data')

#ax[(1,1)].legend(loc='upper right')
#xy = (50,0.1)
xy = (xpos, ypos)
ax[(1,1)].annotate(r'v = 400.0 km/s',xy=xy,xytext=xy,textcoords='data')
xy = (xpos, ypos-0.085)
ax[(1,1)].annotate(r'n$_{\rm{o}}$ = 1.50 cm$^{-3}$',xy=xy,xytext=xy,textcoords='data')

########## Legend and Labels ################
for a in [(0,0), (0,1), (1,0), (1,1)]:
    ax[a].set_xlabel(r'Time (Myr)')
    ax[a].set_ylabel(r'Mass Loss Rate (M$_{\odot}$ Myr$^{-1}$)')

    #ax[a].set_ylim(0.0,1.05)
    ax[a].set_xlim(0.0,2000.0)

    #ax[a].legend(loc='best')
for i in [0,1]:
    ax[(i,0)].set_ylim(0.0, 120.0)
    ax[(i,1)].set_ylim(0.0, 500.0)

ax[(1,0)].legend(loc='lower left')
ax[(1,1)].legend(loc='lower left')

fig.set_size_inches(9,8)
plt.tight_layout()
fig.savefig('SN_dMdt.png')
