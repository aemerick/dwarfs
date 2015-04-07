# this will contain a bunch of stuff to
# solve the anaylytical model we are developing

import numpy as np
import cgs as cgs
from initial_conditions import profiles as prof

# given an orbit of a dwarf galaxy as a function of time
#    - Given: r (position), v_dwarf, n_halo, T_halo 
#    - return r_dwarf and M_gas as a function of that
#             time including: RPS and tidal stripping
#
#
# zeroth order approximation is to calculate the RPS time
# scale at each timestep, and do a linear interpolation of 
# the mass loss i.e.
#  M_i+1 = M_i * (t_i+1 - t_i)/timescale
#
# with a little work, the orbital input can time evolve with
# the RPS and tidal stripping mass loss....
#
#
# 

class analytical_dwarf:

    def __init__(self, name, ic, orbit, potential):
        """ 
            Supply ic's and orbital parameters
        """

        self.ic    = ic
        self.orbit = orbit

    def initialize_gas(self, R, potential):

        if potential == 'NFW':
            self.rho_o = prof.NFW_isothermal_gas(R, r_s=self.ic['b'],
                                                 M200=self.ic['M200'],
                                                 T = self.ic['T_dwarf'], 
                                                 n_o = self.ic['n_o'],
                                                 mu=self.ic['mu_dwarf'],
                                                 rho_crit=self.ic['rho_crit'],
                                      Pcorona = self.ic['n_halo']*cgs.kb*self.ic['T_halo'])

    def gas_profile(self, R, potential='NFW'):
        
        if 'potential' in self.ic.keys():
            potential = self.ic['potential']
            
        if potential == 'NFW':
            return prof.NFW_isothermal_gas(R,    r_s=self.ic['b'],
                                                 M200=self.ic['M200'],
                                                 T = self.ic['T_dwarf'], 
                                                 n_o = self.ic['n_o'],
                                                 mu=self.ic['mu_dwarf'],
                                                 rho_crit=self.ic['rho_crit'],
                                      Pcorona = self.ic['n_halo']*cgs.kb*self.ic['T_halo'])
        
        
    def DM_profile(self, R, potential='NFW'):
        
        if 'potential' in self.ic.keys():
            potential = self.ic['potential']
            
        if potential == 'NFW':
            return prof.NFW_DM(R, r_s=self.ic['b'], M200=self.ic['M200'], rho_crit=self.ic['rho_crit'],
                               decay=False, r_decay=None)
        

def evolve_satellite(t, included_physics, rho_halo, v, rho_gas, M_o, R_o, physics_kwargs={}):
    """
    This function evolves the satellite given a list of desired physical processes
    to model for mass loss in the dwarf galaxy. The included_physics list is currently
    limited to 'RPS' and 'KH'. Physics kwargs should be a dictionary of kwargs with the
    kwargs for a given physics module assigned to the name of the physics, i.e.
    physics_kwargs = { 'RPS': kwargs_dict_for_RPS, 'KH': kwargs_dict_for_KH}.
      
    """
    
    # included physics is going to be a list of the physics "modules" to evovel
    # right now, options should be 'KH' and 'RPS'
    
    # do a check of input parameters. Are they functions or constants?
    
    if not hasattr(rho_halo, '__call__'):
        halo_gas_density = lambda x : rho_halo
        
    if not hasattr(v, '__call__'):
        galaxy_velocity = lambda x : v
        
    # galaxy gas density should be function of radius in galaxy!!!
    if not hasattr(v, '__call__'):
        galaxy_gas_density = lambda x : rho_gas # constant!
    

    KH_const = 0.0; RPS_const = 0.0
    if 'KH' in included_physics:
        KH_const = 1.0
    if not 'KH' in physics_kwargs.keys():
        physics_kwargs['KH'] = {}
    
    
    if 'RPS' in included_physics:
        RPS_const = 1.0
    if not 'RPS' in physics_kwargs.keys():
        physics_kwargs['KH'] = {}
    
    # need to come up with some way to make a function on the fly... constants is fine
    # but if this gets complicated then.... yaa.....
    
    ode_function = lambda y, t :\
                 KH_const *  _KH_evolution(y, t, halo_gas_density, galaxy_velocity,
                                                 galaxy_gas_density,physics_kwargs['KH'])+\
                RPS_const * _RPS_evolution(y, t, halo_gas_density, galaxy_velocity,
                                           galaxy_gas_density, galaxy_gas_density(0.0),physics_kwargs['RPS'])
    # now solve!
    soln = integrate.odeint(ode_function, [M_o, R_o], t)
    
    return soln


def _KH_evolution(y, t, halo_gas_density, galaxy_velocity, galaxy_gas_density):

    Mi = y[0]
    Ri = y[1]
    
    num = 1.0
    if Mi <= 0 or Ri <=0:
        Ri = 0.0
        num = 0.0
        
    Mdot = -np.pi * Ri**2 * halo_gas_denisty(t) * galaxy_velocity(t)
    
    Rdot = Mdot / (4.0*np.pi*Ri**2 * galaxy_gas_density(Ri))
        
    return [Mdot, Rdot]

def RPS_evolution(y, t, halo_gas_density, galaxy_velocity, galaxy_gas_density, rho_gas_o,
                  method='shock', T_galaxy = None, mu_galaxy = None):
    """
        TO DO: make a kwarg to switch between mdot methods...
               this one, nichols, and then using the sound speed scaling 
               
        if method == 'sound' then sound speed determines stripping rate and 
                   temperature and mean molecular weight is needed
    """
    Mi = y[0]
    Ri = y[1]
       
    num = 1.0
    if Mi <= 0 or Ri <=0:
        Ri = 0.0
        num = 0.0
    
    if method == 'shock':
        # using the shock speed to give the RPS timescale
        # taking mdot as current mass divided by timescale
        vshock = (4./3.)*(halo_gas_density(t)/rho_gas_o)**0.5 * galaxy_velocity(t)
        tau_rps = 2.0*Ri/ vshock
        mdot = -Mi / tau_rps      
        rdot = mdot / (4.0*np.pi*Ri**2 * galaxy_gas_density(Ri))
    elif method == 'sound':
        c_sound = ((5.0/3.0)*cgs.kb*T_galaxy/(cgs.mp*mu_galaxy))**0.5
        rdot = -c_sound
        mdot = 4.0*np.pi*Ri**2 * galaxy_gas_density(Ri) * rdot
   
    return [mdot, rdot]



def evolve_1D_satellite(t, r, v, rho_halo, rho_o, r_o, M_gas, T_halo, T_dwarf,
                     rps=True, KH=True, tidal=False):

    nsteps = np.size(t)
    M = np.zeros(nsteps)
    r_dwarf = np.zeros(nsteps)

    # all possible timescales
    timescales_dict = {'RPS'   : (rps, rps_timescale),
                       'KH'    : (KH , KH_timescale),
                       'tidal' : (tidal, tidal_timescale)}


    if np.size(r) == 1:
        r = r * np.ones(nsteps)
    if np.size(rho_o) == 1:
        rho_o = rho_o * np.ones(nsteps)
    if np.size(rho_halo) == 1:
        rho_halo = rho_halo * np.ones(nsteps)
    if np.size(v) == 1:
        v = v* np.ones(nsteps)

    M[0] = M_gas

    for i in np.arange(1,nsteps):
        
        tau_rps   = rps_timescale(rho_halo[i], rho_o[i], v[i], r_o)
        tau_KH    = KH_timescale(M[i-1], rho_halo[i], v[i], T_halo,T_dwarf)
        tau_tidal = tidal_timescale()

        inv_tau = int(rps)/tau_rps + int(KH)/tau_KH + int(tidal)/tau_tidal

        M[i] = M[i-1] * (1.0 - dt * inv_tau)



    return M

def tidal_timescale():

    return 0

def rps_timescale(rho_halo, rho_o, v_dwarf, r_o,
                   n_halo=None, n_dwarf=None,
                   mu_halo=0.61, mu_dwarf=1.30):
    """
        Given the central gas mass density, halo gas mass density,
        dwarf velocity, computes timescale of RPS
    """

    if n_halo is not None:
        rho_halo = n_halo * cgs.mp * mu_halo
    if n_dwarf is not None:
        rho_o    = n_dwarf * cgs.mp * mu_halo

    v_shock = 4.0/3.0 * np.sqrt(rho_halo / rho_o) * v_dwarf

    tau = 2.0 * r_o / v_shock

    return tau

def KH_timescale(M_gas, r_o, rho_halo, rho_dwarf, v_dwarf,
                 T_halo, T_dwarf, n_halo=None, mu_halo=0.61,mu_dwarf=1.31,
                 gamma=5.0/3.0):
    """
    KH timescale for stripping
    """
    if n_halo is not None:
        rho_halo = n_halo *cgs.mp * mu_halo

#    M_rate = np.pi * r_o**2 * rho_halo * v_dwarf
    cs_dwarf = np.sqrt(gamma * cgs.kb * T_dwarf / (cgs.mp*mu_dwarf)) 
    cs_halo  = np.sqrt(gamma * cgs.kb * T_halo / (cgs.mp*mu_halo ))

    M_rate = np.pi * r_o**2 * rho_halo * v_dwarf * (cs_dwarf/cs_halo)

    return M_gas / M_rate
