import matplotlib.pyplot as plt
import numpy as np
import cgs as cgs
#import dwarf_model as dw_model
#from initial_conditions import ic_list as icl ;
from matplotlib import rc

from astropy import units as u

fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.5
point_size = 30


# load all of the simulations from the initial conditions
#LT_names = []
#LT_dict  = {}
#for ic in icl.ic_object_dict.keys():
#    if ic.startswith('LT_n'):
#        LT_names.append(ic)
#        LT_dict[ic] = icl.ic_object_dict[ic]
#

# plot the points
nh3_color = 'blue' ; nh4_color = 'black'; nh5_color = 'red'

cool_ps = 'D' ; SN_ps = 's'; nosn_ps = 'o'


#tstrip_data = np.genfromtxt('tstrip_times.dat', names=True)

data = np.genfromtxt('./../sim_table.dat',names=True,dtype=None)


dwarf_select = 'LT_n150'

for i in np.arange(np.size(data['n'])):

    name = data['name'][i]
    
    if dwarf_select in name:
        if 'nh3' in name:
            color = nh3_color
        elif 'nh4' in name:
            color = nh4_color
        else:
            color = nh5_color

       
        if 'SN' in name:
            fs = color 
            if data['SN'][i]:
                marker = SN_ps
            else:
                marker = cool_ps
        else:
            fs = 'none'
            marker = nosn_ps    

    
        P_RPS = data['nh'][i] * cgs.mp * 0.62 * (data['v'][i]*1.0E5)**2
      
    
        plt.scatter( data['tstrip'][i] , P_RPS, marker = marker, s = point_size , color= color, facecolor=fs)

    
if '150' in dwarf_select:
    dwarf_label = r'n$_{\rm{o}}$ = 1.50 cm$^{-3}$'
elif '75' in dwarf_select:
    dwarf_label = r'n$_{\rm{o}}$ = 0.75 cm$^{-3}$'
else:
    dwarf_label = r'n$_{\rm{o}}$ = 0.20 cm$^{-3}$'

plt.ylim(1.0E-18,5.0E-12)
plt.xlim(0.0,20.0)    


# plot the alphas
alpha = [5.0, 3.0, np.pi/2.0, 1.0]
linestyle = ['-', '--', '-.', ':']
color = ['black','gray','purple']
#n_sample = np.linspace(0.0,2.0,50)
#for i in np.arange(len(alpha)):
# .
#    j = 0
#    for name in ['LT_n075_v2_nh5']:
#        M_o = LT_dict[name].ic['M_HI'] + LT_dict[name].ic['M_DM']
#        r_o = LT_dict[name].ic['r_HI']
#
#        P_anal   = alpha[i] * cgs.G * n_sample * M_o / (r_o) * 1.31 * cgs.mp
#        plt.plot(n_sample, P_anal, ls = linestyle[i], lw = 2.5, color = 'black')
#
#        j = j + 1

#p_stable = np.genfromtxt('./../analytic_model/p_ram_stable.dat',names=True)



#plt.plot(p_stable['n'], p_stable['10'], ls = '-', lw = 2.5, color = 'black')
#plt.plot(p_stable['n'], p_stable['25'], ls = '--', lw = 2.5, color = 'black')
#plt.plot(p_stable['n'], p_stable['50'], ls = '-.', lw = 2.5, color = 'black')
#plt.plot(p_stable['n'], p_stable['100'], ls = ':', lw = 2.5,  color = 'black')

#plot Slater Bell 2014 lines
SB_50 = 10.0**(-12.8)
SB_90 = 10.0**(-14.8)

plt.plot(plt.xlim(), 2*[SB_50], ls='--', lw=2.5, color='green')
plt.plot(plt.xlim(), 2*[SB_90], ls='-', lw=2.5, color='green')

SB_90_t = 2.0
plt.plot(2*[SB_90_t], plt.ylim(), ls='-', lw=2.5, color = 'green')

plt.semilogy()
plt.ylabel(r'RPS Pressure (dyne cm$^{-2}$)')
plt.xlabel(r'Stripping Time (Gyr)')

# do the legend
plt.scatter(-1,-1, color = nh3_color, marker='o', s=point_size, 
                  label = r'n$_{\rm{halo}}$ = 10$^{-3}$ cm$^{-3}$')
plt.scatter(-1,-1, color = nh4_color, marker='o', s=point_size, 
                  label = r'n$_{\rm{halo}}$ = 10$^{-4}$ cm$^{-3}$')
plt.scatter(-1,-1, color = nh5_color, marker='o', s=point_size, 
                  label = r'n$_{\rm{halo}}$ = 10$^{-5}$ cm$^{-3}$')

plt.scatter(-1,-1, color ='black', marker = cool_ps, s= point_size,
                  label = r'Cooling')
plt.scatter(-1,-1, color ='black', marker = SN_ps, s= point_size,
                  label = r'SN')

#
# Put 50 and 90% labels from Slater + Bell lines
#
# vertical:
lx = np.max(plt.xlim()) -2.0
plt.plot(-1,-1,color='green',ls='-',lw=2.5,label = r'Slater \& Bell 2014')
xy = (lx, 1.1*SB_90)
plt.annotate(r'90\%', color = 'green', xy=xy, xytext=xy, textcoords='data')
xy = (lx, 1.1*SB_50)
plt.annotate(r'50\%', color = 'green', xy=xy, xytext=xy, textcoords='data')

# horizontal:
xy = (SB_90_t + 0.3, 2.0E-12)
plt.annotate(r'90\%', color = 'green', xy=xy, xytext=xy, textcoords='data')


xmin,xmax=plt.xlim()
plt.axes().xaxis.set_ticks(np.arange(xmin,xmax+1,1.0),minor=True)

plt.legend(loc='lower left',fancybox=True,ncol=2)

plt.annotate(dwarf_label, xy=(xmax-5.,2.0E-12),xytext=(xmax-5.,2.0E-12),textcoords='data')

plt.savefig(dwarf_select + '_RPS_tstrip.png')

