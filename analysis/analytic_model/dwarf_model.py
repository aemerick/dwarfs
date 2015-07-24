# this will contain a bunch of stuff to
# solve the anaylytical model we are developing

import numpy as np
import cgs as cgs
from initial_conditions import profiles as prof

import copy

from scipy import integrate
from scipy import interpolate


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

    def __init__(self, name, ic):
        """ 
            Initialize analytical dwarf model given a name and a dictionary of initial conditions
        """

        self.ic    = ic
        
        self.R_o   = self.ic['RM']
        integrand = lambda x: x*x*self.gas_profile(x)
        self.M_o   = 4.0*np.pi* integrate.quad(integrand, 0.0, self.R_o)[0]

        
    def setup_orbit(self, t, halo_gas_density, galaxy_velocity):
        """
        Setup orbit parameters of dwarf galaxy... namely the halo gas density dwarf is
        moving through over time, and the galaxy's velocity. Sets these up as 
        functions of time using cubic spline interpolation.
        
        Parameters:
        t  : array
            Array of time values for orbit
        halo_gas_density : float or array
            Array of halo gas densities at the associated times. Requires mass density,
            but makes assumptions to convert to mass density if number density is provided.
            If a single float is provided, a constant function is made
        galaxy_velocity : float or array
            Array of galaxy velocities at the associated times. Makes function using
            cubic spline interpolation. If float, a constant function is made
            
        Returns
        -------
        Nothing. Sets up 'halo_density' and 'galaxy_velocity' parameters.
        """
    
        if any( [halo_gas_density > 1.0E-10] ) : # convert to mass density
            halo_gas_density = halo_gas_density * self.ic['mu_halo'] * cgs.mp
        
        # if t is an array, then use a cubic spline to make a function from the orbital
        # data. If t is a single value, then halo gas dnesity and velocity are constants..
        # make them into functions anyway to make rest of everything work...
        if np.size(halo_gas_density) > 1 : 
            self.halo_density    = interpolate.UnivariateSpline(t, halo_gas_density,k=3)
        else:
            self.halo_density = lambda x: halo_gas_density
        
        if np.size(galaxy_velocity) > 1:
            self.galaxy_velocity = interpolate.UnivariateSpline(t, galaxy_velocity ,k=3)
        else:
            self.galaxy_velocity = lambda x: galaxy_velocity
   

    def gas_profile(self, R, potential='NFW'):
        """
        Evauluate the gas profile at given radii. Gas profile picked from profiles module
        """
        
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
        """
        Evaluate dark matter profile at given radii. DM profile picked from profiles module        
        """
        
        if 'potential' in self.ic.keys():
            potential = self.ic['potential']
            
        if potential == 'NFW':
            return prof.NFW_DM(R, r_s=self.ic['b'], M200=self.ic['M200'], rho_crit=self.ic['rho_crit'],
                               decay=False, r_decay=None)
        
        
    def evolve(self, t, included_physics, physics_kwargs={}, **kwargs):
        """
        Wrapper around evolve_satellite function. Automatically supplies orbital params
        and galaxy profile functions to evolve_satellite based upon the initial conditions 
        already set up in the model.
        
        Parameters
        -----------
        t : array
            Array of time values to evolve over
        included_physics : list
            See 'evolve_satellite' for more details. List of physics to include. Currently
            'RPS' or 'KH'
        physics_kwargs : dict
            Dict of dicts to evolve satellite over.
            
        Returns
        -------
        M, R : array, array
            Arrays of mass and radius evolution over time.
        """
        M, R = \
                evolve_satellite(t, included_physics, self.halo_density, self.galaxy_velocity,
                                              self.gas_profile, self.DM_profile, self.M_o, self.R_o,
                                              physics_kwargs=physics_kwargs, **kwargs)
                
        self.M = M
        self.R = R
        self.t = t
        
        return M, R
        
        
    def load_simulation_mass(self, file_name):
        """
        Loads simulated mass evolution
        """
                
        data = np.genfromtxt(file_name, names=True)
        
        if not hasattr(self, 'simulation_data'):
            self.simulation_data = {}
        
        
        self.simulation_data['mass'] = data['m'] # assuming in solar masses
        self.simulation_data['time'] = data['t'] # assuming in millions of years
        
        
    def evaluate_KH_timescale(self):
        
        self.KH_timescale = np.zeros(np.size(self.R))

        i = 0
        for R,t in zip(self.R,self.t):
            self.KH_timescale[i] = KH_condition_timescale(R, self.DM_profile, self.gas_profile,
                                   self.halo_density(t)/(self.ic['mu_halo']*cgs.mp), 
                                   self.galaxy_velocity(t))
            i = i + 1
            
        return self.KH_timescale

                                   
                                   
def evolve_satellite(t, included_physics, halo_gas_density, galaxy_velocity, galaxy_gas_density, rho_DM, M_o, R_o, physics_kwargs={}, RPS_KH_exclusive = False):
    """
    This function evolves the satellite given a list of desired physical processes
    to model for mass loss in the dwarf galaxy. The included_physics list is currently
    limited to 'RPS' and 'KH'. Physics kwargs should be a dictionary of kwargs with the
    kwargs for a given physics module assigned to the name of the physics, i.e.
    physics_kwargs = { 'RPS': kwargs_dict_for_RPS, 'KH': kwargs_dict_for_KH}.
    
    Parameters
    t : array
        Array of time values to calculate the ODE's over. 
    included_physics: list
        List of physics to include. Current options are:
            'KH' - Kelvin Helmholtz instability.
            'RPS' - Ram Pressure stripping. RPS mass loss rate calculated from
                    rate at which its radius decreases over time, a function
                    of the forward shock.
    halo_gas_density : const or function
        A function of time only that gives the halo gas density as the galaxy
        moves through in its orbit. If a single value is given, assumed constant over
        time
    galaxy_velocity : float or function
        A function of time only giving the galaxy velocity as it moves through its orbit.
        If a single value is given, assumed constant over time
    galaxy_gas_density : float or function
        A function of radius from galaxy center only. Gives galaxy gas density as a 
        function of r from center for a spherical dwarf galaxy model. If constant, assumed
        to be a constant density.
    rho_DM : float or function
        A function of radius from galaxy center only. Gives galaxy dark matter density as a 
        function of r from the center for a spherical dwarf galaxy model. If constnat, assumed
        to be a constant density
    M_o : float
        Initial total gas mass of galaxy
    R_o : float
        Initial (gas) radius of galaxy
    RPS_KH_exclusive: logical, optional
        If set to True, makes RPS and KH instability stripping mutually exclusive. If 
        RPS condition is met, RPS is assumed to dominate mass stripping and KH is turned off.
        KH only occurs, therefore, when RPS condition is not met. Default False.
    
    
        
    Returns
    -------
    M, R: array
        Arrays of gas mass and radius evolution of dwarf galaxy over time.
      
    """
    # included physics is going to be a list of the physics "modules" to evovel
    # right now, options should be 'KH' and 'RPS'

    physics_kwargs_copy = copy.deepcopy(physics_kwargs)
    # do a check of input parameters. Are they functions or constants?
    
    if not hasattr(halo_gas_density, '__call__'):
        halo_gas_density = lambda x : halo_gas_density
        
    if not hasattr(galaxy_velocity, '__call__'):
        galaxy_velocity = lambda x : galaxy_velocity
        
    # galaxy gas density should be function of radius in galaxy!!!
    if not hasattr(galaxy_gas_density, '__call__'):
        galaxy_gas_density = lambda x : galaxy_gas_density # constant!
    
    # assume KH and RPS are off unless in list of included physics
    KH_const = 0.0; RPS_const = 0.0
    
    if 'KH' in included_physics:
        KH_const = 1.0

    if not 'KH' in physics_kwargs_copy.keys(): # bookkeeping if off
        physics_kwargs_copy['KH'] = {}
    
    if 'RPS' in included_physics:
        RPS_const = 1.0
        
    if not 'RPS' in physics_kwargs_copy.keys(): # bookkeeping if off
        physics_kwargs_copy['RPS'] = {}
    
    # if alpha is contained in physcis kwargs... strip it to be 
    # used in the RPS condition function call, as it is not used in the
    # RPS mass loss rate calculation
    if 'alpha' in physics_kwargs_copy['RPS'].keys():
        alpha = physics_kwargs_copy['RPS']['alpha']
        physics_kwargs_copy['RPS'].pop('alpha',None)
    else:
        alpha = 1.0
        
    # need to come up with some way to make a function on the fly... constants is fine
    # but if this gets complicated then.... yaa.....
    
    ode_function = lambda y, t, A, B:\
                 A *  _KH_evolution(y, t, halo_gas_density, galaxy_velocity,
                                                 galaxy_gas_density, **physics_kwargs_copy['KH'])+\
                 B * _RPS_evolution(y, t, halo_gas_density, galaxy_velocity,
                                           galaxy_gas_density,
                                           galaxy_gas_density(0.0), **physics_kwargs_copy['RPS'])
    
    # write a loop here... solve step by step
    M = np.zeros(np.size(t)); R = np.zeros(np.size(t))
    M[0] = M_o; R[0] = R_o
    keep_looping = True; i = 0; turn_KH_off = 1.0 
    while (i < np.size(t) - 1) and keep_looping:
        
        # check if ram pressure stripping occurs
        if 'RPS' in included_physics:
            # integrate and test around the current radius
#            rps_cond = _RPS_condition(np.linspace(0.9999*R[i],1.0001*R[i],5), rho_DM, galaxy_gas_density, 
#                                               halo_gas_density(t[i]), galaxy_velocity(t[i]), alpha=alpha)

            rps_cond = _RPS_condition(R[i], rho_DM, galaxy_gas_density, halo_gas_density(t[i]),
                                                    galaxy_velocity(t[i]), alpha=alpha)
            
            # if RPS is valid at current radius, use it... otherwise set to zero
            if rps_cond  > 0:
                RPS_const = 1.0
            else:
                RPS_const = 0.0 
        
            if RPS_KH_exclusive and RPS_const == 1.0:   # turn KH off
                turn_KH_off = 0.0
            elif RPS_KH_exclusive and RPS_const == 0.0: # turn KH on
                turn_KH_off = 1.0
            else:                          # else just keep it the same
                turn_KH_off = KH_const 
               
        ode_function_args = (KH_const * turn_KH_off, RPS_const,)
        
       
        
        soln = integrate.odeint(ode_function, [M[i],R[i]], t[i:i+2], 
                                    args=ode_function_args,
                                    mxhnil=0, ixpr=False)
        M[i+1] = soln[1,0]; R[i+1] = soln[1,1]
        
        i = i + 1
        
        simple_check = M[i] + ode_function([M[i],R[i]], t[i], *ode_function_args)[0] * (t[i] - t[i-1])
        
        if M[i] <= 0.0 or R[i] <= 0.0 or simple_check <= 0.0:
            M[i] = 0.0; R[i] = 0.0
            keep_looping = False
            
    
    return M, R

def _RPS_condition(r, DM_density, gas_density, halo_density, galaxy_velocity, alpha=1.0):
    # need DM density profile, gas density profile, and 
    
    # find the mass profile from density profile
    DM_integrand = lambda x: x*x*DM_density(x)
    gas_integrand = lambda x: x*x*gas_density(x)

    return_float = False
    if np.size(r) == 1:
        r = np.array([0.0,r])
        return_float = True

    M_DM      = np.zeros(np.size(r)-1)
    M_gas     = np.zeros(np.size(r)-1)
    
    for i in np.arange(1,np.size(r)):
        M_DM[i-1] = 4.0*np.pi* integrate.quad(DM_integrand, 0.0, r[i])[0]
        M_gas[i-1] = 4.0*np.pi* integrate.quad(gas_integrand, 0.0,r[i])[0]
    
    # Now, calculate the RHS and LHS of the stripping condition
    RHS = (M_DM + M_gas) * gas_density(r[1:]) / r[1:] * cgs.G * alpha
    LHS = halo_density * galaxy_velocity * galaxy_velocity
    
    if return_float:
        return np.float((LHS-RHS)[0])
    else:
        return LHS - RHS
    

def _KH_evolution(y, t, halo_gas_density, galaxy_velocity, galaxy_gas_density, beta=1.0):

    Mi = y[0]
    Ri = y[1]
    
    num = 1.0
    if Mi <= 0 or Ri <=0:
        Mdot = 0.;Rdot=0.
    else:
        Mdot = -np.pi * Ri**2 * halo_gas_density(t) * galaxy_velocity(t)
        Rdot = Mdot / (4.0*np.pi*Ri**2 * galaxy_gas_density(Ri))
        
    return np.array([Mdot*beta, Rdot*beta])

def _RPS_evolution(y, t, halo_gas_density, galaxy_velocity, galaxy_gas_density, rho_gas_o,
                  method='shock', T_galaxy = None, mu_galaxy = None, beta = 1.0):
    """
        This function can be supplied to an ode solver (taken care of by the 
        evolve_satellite function) to model mass loss and radius change of a 
        spherical galaxy undergoing ram pressure stripping. RPS mass loss
        rate is calculated in two methods, set by the 'method' kwarg.
        
        Parameters
        -----------
        y   :  2 element array or list
            Two element array or list with the first element containing mass 
            and the second radius
        t   : array
            Dependent variable for ODE. Time in cgs
        halo_gas_density : function
            Function of time only that gives the gas density of the halo
            the galaxy is moving through.
        galaxy_velocity: function
            Function of time only that gives the velocity of the galaxy
        galaxy_gas_density: function
            Function of R only. The spherically symettric gas (mass) density
            profile of the dwarf galaxy. 
        rho_gas_o : real
            Value of the above at R = 0 (central gas density).
        
        Returns
        -------
        [Mdot, Rdot] : List 
            
        
        
    
        TO DO: make a kwarg to switch between mdot methods...
               this one, nichols, and then using the sound speed scaling 
               
        if method == 'sound' then sound speed determines stripping rate and 
                   temperature and mean molecular weight is needed
    """
    Mi = y[0]
    Ri = y[1]
       
    num = 1.0
    if Mi <= 0 or Ri <=0:
        mdot=0.; rdot=0.
    else:
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

    return np.array([mdot*beta, rdot*beta])



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

def KH_timescale(M_gas, r_o, rho_halo, v_dwarf,
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

def KH_condition_timescale(r, DM_density, gas_density, n_halo, v_dwarf,
                                  mu_halo=0.61, mu_dwarf=1.31):
   
    # find the mass profile from density profile
    DM_integrand = lambda x: x*x*DM_density(x)
    gas_integrand = lambda x: x*x*gas_density(x)
    
    M_DM  = 4.0*np.pi* integrate.quad(DM_integrand,  0.0, r)[0]
    M_gas = 4.0*np.pi* integrate.quad(gas_integrand, 0.0, r)[0]       
        
    F = M_gas / (M_DM + M_gas)    
    
    M_DM = M_DM / cgs.Msun
    v_dwarf = v_dwarf / cgs.km
        
        
    # From McCarthy et. al. 2008 ... MDM is in Msun below, rho_halo in cc .. v_dwarf in km/s  
    # gives timescale in years
    tau_kh = 2.19E9*(F/0.1) * (M_DM / 1.0E9)**(1.0/7.0) * (n_halo / 1.0E-4)**(-1.0) * (v_dwarf/1.0E3)**(-1)
  
    return tau_kh * cgs.yr
