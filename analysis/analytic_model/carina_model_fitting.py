from __future__ import print_function

# coding: utf-8

# In[2]:

#get_ipython().magic(u'matplotlib inline')

from scipy import interpolate
import numpy as np
import cgs as cgs
import dwarf_model as dw_model
from initial_conditions import ic_list as icl ;
from matplotlib import rc
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import sys
#from __future__ import print_function


# In[3]:

carina = icl.ic_object_dict['CarinaMidMed']

carina_orbit = np.genfromtxt("./../orbits/carina_orbit_tab.dat")
t_orbit = carina_orbit[:,0] * cgs.Myr
v_orbit = carina_orbit[:,2] * 1.0E5
rho_halo = carina.ic['n_halo'] * cgs.mp * carina.ic['mu_halo']

anal_carina = dw_model.analytical_dwarf('CarinaMidMed',carina.ic)
anal_carina.setup_orbit(t_orbit,rho_halo,v_orbit)


# In[4]:

#anal_carina.ic


# In[5]:

fpath="/home/emerick/Research/dwarfs/flash_runs/"
adiabatic = np.genfromtxt(fpath + "carina_adiabatic/carina_adiabatic_mass_dt5Myr.dat",names=True)
adiabatic['t'] = adiabatic['t'] * cgs.Myr
adiabatic['m'] = adiabatic['m'] * cgs.Msun


# In[6]:


def function_to_fit(x, p0, p1, p2, dwarf, method):

    def _mass_function(a,rpsb,khb,t,d,m):
        M_fit, R_fit = d.evolve(t, ['RPS','KH'],
                                   physics_kwargs={'RPS':{'alpha':a,'beta':rpsb, 'method':m,
                                                                          'T_galaxy':d.ic['T_dwarf'],
                                                                          'mu_galaxy':d.ic['mu_dwarf']},
                                                   'KH':{'beta':khb}})
        return M_fit
    # parameters for fit
    alpha     = p0
    rps_beta  = p1
    kh_beta   = p2
    
    t = np.linspace(0.0,np.max(t_orbit)/cgs.Myr,2000.0)*cgs.Myr
    # evolve to get m as function of time
    m = _mass_function(alpha,rps_beta,kh_beta,t,dwarf,method)
    
    # recast as callable function
    function = interpolate.UnivariateSpline(t, m,k=3)
    
    return function(x)


# In[7]:

#func = lambda x,alpha,beta: function_to_fit(x,alpha,beta,anal_carina,method='shock')

#popt, pconv = curve_fit(func, adiabatic['m'], adiabatic['t'], p0=[0.75,25.0])


# In[8]:

# time to do some manual fitting
alphas     = np.linspace(0.1,    0.7,  20)
rps_betas  = np.linspace(0.01,  0.25,  20)
kh_betas   = np.linspace(0.95,  1.25,  20)


# In[ ]:

nparams = np.size(alphas)*np.size(rps_betas)*np.size(kh_betas)
shock_error  = np.zeros(nparams)
shock_params = np.array([[0.,0.,0.]]*nparams)*1.0
i=0
f = open('shock_parameter_study.dat','w')
f.write("#alpha rps_beta kh_beta error\n")
out_format = "%8.8E %8.8E %8.8E %8.8E\n"
for a in alphas:
    for rps_b in rps_betas:
        for kh_b in kh_betas:
            m = function_to_fit(adiabatic['t'],a,rps_b,kh_b,anal_carina,method='shock');
            shock_error[i] = np.sum((m - adiabatic['m'])**2)
            shock_params[i][0] = a * 1.0
            shock_params[i][1] = rps_b * 1.0
            shock_params[i][2] = kh_b*1.0
            print('%i of %i - alpha = %.3f , rps_beta = %.3f, kh_beta %.3f, error = %.5e'%(i+1,nparams,a,rps_b,kh_b,shock_error[i]),file=sys.stderr)
            f.write(out_format%(a,rps_b,kh_b,shock_error[i]))
            i = i + 1
f.close()
shock_params = np.array(shock_params)


# In[40]:

print(np.argmin(shock_error), np.min(shock_error) , np.median(shock_error))
print(shock_params[np.argmin(shock_error)])

min_err_index = np.argmin(shock_error)
best_shock_alpha = shock_params[min_err_index][0]
best_shock_beta_rps  = shock_params[min_err_index][1]
best_shock_beta_kh = shock_params[min_err_index][2]
string = "best alpha  = %.4f - best beta = %.4f - best kh = %.4f\n"%(best_shock_alpha,best_shock_beta_rps, best_shock_beta_kh)
print(string)
outfile = open('outfile.dat','w')
outfile.write('shock\n')
outfile.write(string)
outfile.close()

"""

nnparams = np.size(alphas)*np.size(rps_betas)*np.size(kh_betas)
sound_error  = np.zeros(nparams)
sound_params = np.array([[0.,0.,0.]]*nparams)*1.0
i=0
f = open('sound_parameter_study.dat','w')
f.write("#alpha rps_beta kh_beta error\n")
out_format = "%8.8E %8.8E %8.8E %8.8E\n"
for a in alphas:
    for rps_b in rps_betas:
        for kh_b in kh_betas:
            m = function_to_fit(adiabatic['t'],a,rps_b,kh_b,anal_carina,method='sound');
            sound_error[i] = np.sum((m - adiabatic['m'])**2)
            sound_params[i][0] = a * 1.0
            sound_params[i][1] = rps_b * 1.0
            sound_params[i][2] = kh_b*1.0
            print('%i of %i - alpha = %.3f , rps_beta = %.3f, kh_beta = %.3f, error = %.5e'%(i+1,nparams,a,rps_b,kh_b,sound_error[i]),file=sys.stderr)
            f.write(out_format%(a,rps_b,kh_b,sound_error[i]))
            i = i + 1
f.close()
sound_params = np.array(sound_params)


# In[ ]:

print(np.argmin(sound_error), np.min(sound_error) , np.median(sound_error))
print(sound_params[np.argmin(sound_error)])

min_err_index = np.argmin(sound_error)
best_sound_alpha = sound_params[min_err_index][0]
best_sound_beta_rps  = sound_params[min_err_index][1]
best_sound_beta_kh = sound_params[min_err_index][2]

string = "best alpha  = %.4f - best beta = %.4f - beta kh - %.4f\n"%(best_sound_alpha,best_sound_beta_rps,best_sound_beta_kh)
print(string)
outfile.write('sound\n')
outfile.write(string)


# In[ ]:

alpha_shock = best_shock_alpha
rps_beta_shock  = best_shock_beta_rps#0.18181818 #best_beta
kh_beta_shock = best_shock_beta_kh# 0.1 #best_kh_beta
alpha_sound = best_sound_alpha #0.5         0.18181818  0.1  
rps_beta_sound = best_sound_beta_rps
kh_beta_sound = best_sound_beta_kh

m_shock = function_to_fit(adiabatic['t'],alpha_shock,rps_beta_shock,kh_beta_shock,anal_carina,method='shock')
m_sound = function_to_fit(adiabatic['t'],alpha_sound,rps_beta_sound,kh_beta_sound,anal_carina,method='sound')

error_shock = np.sum(((m_shock-adiabatic['m'])/adiabatic['m'][0])**2)
error_sound = np.sum(((m_sound-adiabatic['m'])/adiabatic['m'][0])**2)

plt.plot(adiabatic['t']/cgs.Myr, adiabatic['m']/adiabatic['m'][0],color='red',lw=2,label='sim')
plt.plot(adiabatic['t']/cgs.Myr, m_shock/m_shock[0] ,color='black',lw=2,label='model - shock')
plt.plot(adiabatic['t']/cgs.Myr, m_sound/m_sound[0] ,color='blue',lw=2,label='model - sound')

plt.ylim(0,1)
plt.legend(loc='best')


# In[ ]:

outfile.close()


# In[ ]:

f = open('alpha_beta_rps.dat','w')
string = "%.6f %.6f %.6f %8.8e %.6f %.6f %.6f %8.8e\n"

f.write("alpha_shock rps_beta_shock kh_beta_shock error_shock alpha_sound rps_beta_sound kh_beta_sound error_sound\n")
for i in np.arange(np.size(shock_error)):
    f.write(string%(shock_params[i][0], shock_params[i][1], shock_params[i][2], shock_error[i], sound_params[i][0], sound_params[i][1], sound_params[i][2], sound_error[i]))
    
f.close()


# In[ ]:




# In[ ]:


"""
