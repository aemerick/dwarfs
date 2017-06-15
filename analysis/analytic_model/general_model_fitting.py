from __future__ import print_function

from scipy import interpolate
import numpy as np
import cgs as cgs
import dwarf_model as dw_model
from initial_conditions import ic_list as icl ;
from matplotlib import rc
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import sys


simulation_dir = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'
model_name = "LT_n075_v2_nh4"
sim_mass   = model_name + '_cont_mass_1-0.dat'

outfolder = '/home/emerick/Research/dwarfs/flash_runs/leo_T/analytical_model_fits/'


initial_conditions = icl.ic_object_dict[model_name.replace('v4','v2')]
anal_model = dw_model.analytical_dwarf(model_name, initial_conditions.ic)
anal_model.load_simulation_mass(simulation_dir +model_name +'/' + sim_mass)

if 'v4' in model_name:
    velocity = 400.0E5
elif 'v2' in model_name:
    velocity = 200.0E5

anal_model.setup_orbit(0.0, initial_conditions.ic['n_halo'], velocity)

t_vals = anal_model.simulation_data['time'] * cgs.Myr
m_vals = anal_model.simulation_data['mass']

zeros = np.where(m_vals <= 100.0)

if np.size(zeros) >= 1:
   
    max_fit_time = np.float(t_vals[zeros[0][0]])
else:
    max_fit_time = np.float(t_vals[-1])

# In[6]:


def function_to_fit(x, p0, p1, p2, dwarf, method, t_max = None, exclusive=False):

    def _mass_function(a,rpsb,khb,t,d,m,rpskhexclusive):
        M_fit, R_fit = d.evolve(t, ['RPS','KH'],
                                   physics_kwargs={'RPS':{'alpha':a,'beta':rpsb, 'method':m,
                                                                          'T_galaxy':d.ic['T_dwarf'],
                                                                          'mu_galaxy':d.ic['mu_dwarf']},
                                                   'KH':{'beta':khb}},
                                                   RPS_KH_exclusive = rpskhexclusive)
        return M_fit
    # parameters for fit
    alpha     = p0
    rps_beta  = p1
    kh_beta   = p2
    
    
    npoints = int(t_max/cgs.Myr)
    
    t = np.linspace(0.0, t_max, npoints)
    
    # evolve to get m as function of time
    m = _mass_function(alpha,rps_beta,kh_beta,t,dwarf,method,exclusive)
    
    # recast as callable function
    function = interpolate.UnivariateSpline(t, m, k=3)
    
    return function(x)


# time to do some manual fitting
alphas     = np.linspace(0.25,  4.0,   100)
rps_betas  = np.linspace(0.5,   2.0,  100)
kh_betas   = np.linspace(1.0,   6.0 , 100)




nparams = np.size(alphas)*np.size(rps_betas)*np.size(kh_betas)
shock_error  = np.zeros(nparams)
shock_params = np.array([[0.,0.,0.]]*nparams)*1.0
i=0
f = open(outfolder + model_name + '_shock_param.dat','w')
f.write("#alpha rps_beta kh_beta error\n")
out_format = "%.4f %.4f %.4f %.14E\n"
for a in alphas:
    for rps_b in rps_betas:
        for kh_b in kh_betas:
            m = function_to_fit(anal_model.simulation_data['time'] * cgs.Myr ,a,rps_b,kh_b, anal_model, method='shock', t_max = max_fit_time, exclusive = True);
            shock_error[i] = np.sum((m/cgs.Msun - anal_model.simulation_data['mass'] )**2)
            shock_params[i][0] = a * 1.0
            shock_params[i][1] = rps_b * 1.0
            shock_params[i][2] = kh_b*1.0
            #print('%i of %i - alpha = %.3f , rps_beta = %.3f, kh_beta %.3f, error = %.10e'%(i+1,nparams,a,rps_b,kh_b,shock_error[i]),file=sys.stderr)
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
outfile = open(outfolder + model_name + '_best_shock_param.dat','w')
outfile.write('shock\n')
outfile.write(string)
outfile.close()
