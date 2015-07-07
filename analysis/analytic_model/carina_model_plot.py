import numpy as np
import matplotlib.pyplot as plt
import dwarf as dw
import dwarf_model as dw_model
import cgs as cgs
from initial_conditions import ic_list as icl ;
from scipy import interpolate

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
#best alpha  = 0.2500 - best beta = 0.1000 - best kh = 1.1455
#sound
#best alpha  = 1.0455 - best beta = 0.8000 - beta kh - 1.7500


# parameters
shock_par = [0.2500,0.1000,1.1455]
sound_par = [1.0455,0.8000,1.7500]

# load IC's for carina
carina = icl.ic_object_dict['CarinaMidMed']

# load orbit (not needed in this script yet)
carina_orbit = np.genfromtxt("./../orbits/carina_orbit_tab.dat")
t_orbit = carina_orbit[:,0] * cgs.Myr
v_orbit = carina_orbit[:,2] * 1.0E5
rho_halo = carina.ic['n_halo'] * cgs.mp * carina.ic['mu_halo']

anal_carina = dw_model.analytical_dwarf('CarinaMidMed',carina.ic)
anal_carina.setup_orbit(t_orbit,rho_halo,v_orbit)

# load adiabatic simulation data
fpath="/home/emerick/Research/dwarfs/flash_runs/"
adiabatic = np.genfromtxt(fpath + "carina_adiabatic/carina_adiabatic_mass_dt5Myr.dat",names=True)
adiabatic['t'] = adiabatic['t'] * cgs.Myr
adiabatic['m'] = adiabatic['m'] * cgs.Msun

# plot
alpha_shock = shock_par[0]
rps_beta_shock  = shock_par[1]
kh_beta_shock = shock_par[2]
alpha_sound = sound_par[0]
rps_beta_sound = sound_par[1]
kh_beta_sound = sound_par[2]

m_shock = function_to_fit(adiabatic['t'],alpha_shock,rps_beta_shock,kh_beta_shock,anal_carina,method='shock')
m_sound = function_to_fit(adiabatic['t'],alpha_sound,rps_beta_sound,kh_beta_sound,anal_carina,method='sound')

error_shock = np.sum(((m_shock-adiabatic['m'])/adiabatic['m'][0])**2)
error_sound = np.sum(((m_sound-adiabatic['m'])/adiabatic['m'][0])**2)

plt.plot(adiabatic['t']/cgs.Myr, adiabatic['m']/adiabatic['m'][0],color='red',lw=2,label='sim')
plt.plot(adiabatic['t']/cgs.Myr, m_shock/m_shock[0] ,color='black',lw=2,label='model - shock')
plt.plot(adiabatic['t']/cgs.Myr, m_sound/m_sound[0] ,color='blue',lw=2,label='model - sound')

plt.ylim(0,1)
plt.legend(loc='best')
plt.savefig('carina_model_fitting.png')
