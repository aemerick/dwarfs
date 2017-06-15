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
    
####### SN_LT_n075_v2_nh4 #########
name  = 'SN_LT_n075_v2_nh4/'
fname = fpath + name + '0000_cfloor_nsn_redo/' + file_name

data  = np.genfromtxt(fname, names=True)
ax[(0,0)].plot(data['t'], data['m']/data['m'][0], lw=lw, color='black',
                          ls = nsn_ls, label = 'No SN')
print_results(fname, data)




#fname = fpath + name + '0000_cfloor_nsn/' + file_name
#data  = np.genfromtxt(fname, names=True)
#ax[(0,0)].plot(data['t'], data['m']/data['m'][0], lw=lw, color='red',
                          #ls = nsn_ls, label = 'No SN')
#print_results(fname, data)


#fname = fpath + name + '0000_cfloor_global/' + file_name
#data = np.genfromtxt(fname, names=True)
#ax[(0,0)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='black',
#                          ls = sn_ls, label = 'SN')
#print_results(fname, data)

fname = fpath + name + '0000_cfloor_global_redo/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(0,0)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='black',
                          ls = sn_ls, label = 'SN')
print_results(fname, data)


####### SN_LT_n075_v4_nh4 #########
name  = 'SN_LT_n075_v4_nh4/'

fname = fpath + name + '0000_cfloor_nsn/' + file_name
data  = np.genfromtxt(fname, names=True)
ax[(0,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color='black',
                          ls = nsn_ls, label = 'No SN')
print_results(fname, data)


fname = fpath + name + '0000_cfloor_global/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(0,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='black',
                          ls = sn_ls, label = 'SN')
print_results(fname, data)

fname = fpath + name + '0000_cfloor_global_x2/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(0,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='blue',
                          ls = ':', label = 'SN x 2')
print_results(fname, data)

####### SN_LT_n150_v2_nh4 #########
name  = 'SN_LT_n150_v2_nh4/'

fname = fpath + name + '0000_cfloor_nsn/' + file_name
data  = np.genfromtxt(fname, names=True)
ax[(1,0)].plot(data['t'], data['m']/data['m'][0], lw=lw, color='black',
                          ls = nsn_ls, label = 'No SN')
print_results(fname, data)


fname = fpath + name + '0000_cfloor_global/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(1,0)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='black',
                          ls = sn_ls, label = 'SN')
print_results(fname, data)

####### SN_LT_n150_v4_nh4 #########
name  = 'SN_LT_n150_v4_nh4/'

fname = fpath + name + '0000_cfloor_nsn/' + file_name
data  = np.genfromtxt(fname, names=True)
ax[(1,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color='black',
                          ls = nsn_ls)#, label = 'No SN')
print_results(fname, data)


fname = fpath + name + '0000_cfloor_global/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(1,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='black',
                          ls = sn_ls)#, label = 'SN')
print_results(fname, data)

fname = fpath + name + '0000_cfloor_global_x2/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(1,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='blue',
                          ls = ':', label = 'SN x 2')
print_results(fname, data)

fname = fpath + name + '0000_cfloor_global_x5/' + file_name
data = np.genfromtxt(fname, names=True)
ax[(1,1)].plot(data['t'], data['m']/data['m'][0], lw=lw, color ='orange',
                          ls = '-.', label = 'SN x 5')
print_results(fname, data)


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
    ax[a].set_ylabel(r'M(t) / M$_{\rm{o}}$')

    ax[a].set_ylim(0.0,1.05)
    ax[a].set_xlim(0.0,2000.0)

    #ax[a].legend(loc='best')

    ax[a].minorticks_on()
    
ax[(1,0)].legend(loc='lower left')
ax[(1,1)].legend(loc='lower left')

fig.set_size_inches(9,8)
plt.tight_layout()
fig.savefig('SN_mass_loss.png')
