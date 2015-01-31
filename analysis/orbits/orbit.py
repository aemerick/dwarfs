from __future__ import division

import numpy as np
import cgs   as cgs


class dwarf_orbit:
    """ 
        dwarf orbit object for integrating dwarf orbit and calculating
        density, temperature, and velocity as a function of orbit
        time / distance to GC
    """

    def __init__(self, name):
        self.name = name

    def set_x0(self,x0):
        """ 
            sets initial position in cartesian coordinates
            as a numpy array (x0)
        """
            self.x0 = np.array(x0)

    def set_v0(self, v0):
        """
            initial cartesian velocity as a 3D numpy array
        """
            self.v0 = np.array(v0)

    def calculate_initials(self):
        print "does nothing... calculate x0 and v0 from obs"

    def integrate_orbit(self, t_end = 500.0*cgs.Myr, dt=1.0E9, **kwargs):
        """
            integrate the orbit given some dt and nsteps
            in order for this to work, must have x0, v0, and 
            accel_function set

            Parameters
            ----------
            t_end : float, optional
                End time of integration in seconds. Default 500 Myr
            dt    : float, optional
                Timestep size in seconds. Default around 0.1 kyr
        """
        nsteps = np.ceil(t_end / dt)

        print "integrating orbit for " + self.name
        print "for %5.4e Myr"%(t_end/cgs.Myr)
        print "Using %2.2e timesteps at dt = %5.4e"%(nsteps,dt)

        t,x,v = leapfrog_integrate(self.acceleration_function, self.x0,
                           self.v0, dt, nsteps, kwargs)


        self.t = t
        self.x = x
        self.v = v

        self.r  = np.sqrt(np.sum(x**2, axis=-1))
        self.vr = np.sqrt(np.sum(v**2, axis=-1))

    def calculate_density(self,density_type="MB13_adj", mu=cgs.mu **kwargs):
        """
            returns the number density and density assuming
            ionized primordial gas for the provided density type.
            Options are MB13, MB13_adj, Gato_isothermal, 
            Gato_adiabatic, Gato_cooling, Kaufmann_realistic

            ONLY "MB13", "MB13_adj" are function profiles, the rest are
            interpolated from data. This means there is a rmin and rmax
            allowed for the other profiles.
        """

        if density_type == 'MB13':
            n = gal.halo_gas_ndensity(r,kwargs)

        elif density_type == 'MB13_adj':
            n = gal.halo_gas_ndensity(r, n0=.46+.74, rc=(.35+.29)*cgs.kpc,
                                        beta=.74-.14, ambient=0.0)

        elif density_type == 'Gato_isothermal':
            n = gal.gato_ndensity(r,'isothermal')
        elif density_type == 'Gato_adiabatic':
            n = gal.gato_ndensity(r,'adiabatic')
        elif density_type == 'Gato_cooling':
            n = gal.gato_ndensity(r,'cooling')
        elif density_type == 'Kaufmann_realistic':
            n = gal.kaufmann('02', 'density')

        self.n = n
        self.rho = n * cgs.mp * mu
                         
    def calculate_temperature(self, T_type='NFW', **kwargs):
        """
            Calculates the gas temperature at the intergated
            positions for the given model.
        Parameters:
        -----------
            T_type : string, optional
            Temperature profile type. Defuault is the profile from virial
            equillibrium with a spherical NFW DM halo
        """
   
        if T_type == 'NFW':
            T = gal.halo_gas_temperature(r, kwargs)
        elif T_type == 'Kaufmann_realistic':
            T = gal.kaufmann('02', 'temperature')

        self.T = T
 
    def set_acceleration(self, function):
        """
            sets the desired acceleration function where function 
            must be a function with variables t and x. This will
            be supplied to leapfrog integrator.
        """
        self.acceleration = function
    
    def write_out_wind(self, filename=self.name + "_wind.dat"):
        """
        Write out the integrated orbital time, radius, total velocity,
        number density, density, and temperature at all timesteps for 
        the integrated orbit. This file is to be read into FLASH
        for establishing the wind flow conditions
   
        Parameters
        ----------
        filename: string, optional
            Name of output file. Default to object name with "_wind.dat" ext
        """
        ok_to_write = True
        
        if not hasattr(self, 'n'):
           ok_to_write = False
           print "Need to calculate density before writing out"
        if not hasattr(self, 'T'):
           ok_to_write = False
           print "Need to calculate temperature before writing out"


        if ok_to_write:
            f = f.open(filename)
        
            fmt = "%8.8E %8.8E %8.8E %8.8E %8.8E %8.8E\n"
            for i in np.arange(np.size(r)):
                f.write(fmt%( self.t[i],   self.r[i], self.vr[i],
                              self.n[i], self.rho[i], self.T[i])
            f.close()
        
def spherical_NFW(t, xyz, c=12, M200=1.0E12*cgs.Msun , rho_crit = 9.74E-30):
    """
    acceleration function for the spherical NFW profile
    """


    R200 = (3.0 * M200 / (4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
#    fc   = np.log(1.0 + c) - c/(1.0+c)    
    rscale = c / R200

    r = np.sqrt(np.sum(xyz**2, axis=-1))
    val = np.log( r/rscale + 1.0)/r**1.5 - 1.0/(r**2.5 + rscale)

    return val[:,np.newaxis] * xyz * cgs.G * M200

def leapfrog_integrate(acceleration_func, x0, v0, dt, nsteps, t1=0., args=()):
    """ A simple implementation of Leapfrog integration 
    
    Parameters
    ----------
    acceleration_func : callable
        A function that computes the acceleration at a given position.
    x0 : numpy.ndarray
        Initial position(s).
    v0 : numpy.ndarray
        Initial velocity(s).
    dt : numeric
        Timestep.
    nsteps : int
        Number of steps to run for.
    t1 : numeric (optional)
        Initial time.
    args : iterable (optional)
        Any other arguments the acceleration function takes.
        
    Returns
    -------
    t : numpy.ndarray
        Array of times.
    x : numpy.ndarray
        Array of positions.
    v : numpy.ndarray
        Array of velocities.
    """
    
    # ensure that the initial conditions are arrays and at least 2D
    x0 = np.atleast_2d(x0).copy()
    v0 = np.atleast_2d(v0).copy()
    norbits,ndim = x0.shape
    
    # wrapper around the acceleration function so we can call it with just the position, x
    acc = lambda t,x: acceleration_func(t,x,*args)
    
    all_x = np.zeros((nsteps,norbits,ndim))
    all_v = np.zeros((nsteps,norbits,ndim))
    t = np.zeros(nsteps)
    
    all_x[0] = x0
    all_v[0] = v0
    
    # velocity at 1/2 step 
    v_iminus1_2 = v0 + acc(t1, x0)*dt/2.
    x_iminus1 = x0.copy()
    for i in range(1,nsteps):
        t[i] = t[i-1] + dt
        x_i = x_iminus1 + v_iminus1_2*dt # full step
        a_i = acc(t[i], x_i)
        v_i = v_iminus1_2 + a_i*dt/2. # half step
        v_iplus1_2 = v_i + a_i*dt/2. # half step
        
        all_x[i] = x_i
        all_v[i] = v_i
        
        x_iminus1 = x_i
        v_iminus1_2 = v_iplus1_2
    
    return t, all_x, all_v

