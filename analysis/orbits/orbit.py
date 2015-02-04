from __future__ import division

import numpy as np
import cgs   as cgs
from halo import galaxy as gal

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

    def integrate_orbit(self, t_end = 1.0E17, dt=1.0E11,
                        verbose=True, **kwargs):
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
        nsteps = int(np.ceil(t_end / dt))

        print "integrating orbit for " + self.name
        print "for %5.4e Myr"%(t_end/cgs.Myr)
        print "Using %2.2e timesteps at dt = %5.4e"%(nsteps,dt)

    

        t,x,v = leapfrog_integrate(self.acceleration_function, self.x0,
                           self.v0, dt, nsteps, verbose, kwargs)


        self.t = t
        self.x = x
        self.v = v

        self.r  = np.sqrt(np.sum(x**2, axis=-1)).flatten()
        self.vr = np.sqrt(np.sum(v**2, axis=-1)).flatten()

    def calculate_density(self,density_type="MB13_adj", mu=cgs.mu, **kwargs):
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
            n = gal.halo_gas_ndensity(self.r,kwargs)

        elif density_type == 'MB13_adj':
            n = gal.halo_gas_ndensity(self.r, n0=.46+.74, rc=(.35+.29)*cgs.kpc,
                                        beta=.74-.14, ambient=0.0)

        elif density_type == 'Gato_isothermal':
            n = gal.gato_ndensity(self.r,'isothermal')
        elif density_type == 'Gato_adiabatic':
            n = gal.gato_ndensity(self.r,'adiabatic')
        elif density_type == 'Gato_cooling':
            n = gal.gato_ndensity(self.r,'cooling')
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
            T = gal.halo_gas_temperature(self.r, **kwargs)
        elif T_type == 'Kaufmann_realistic':
            T = gal.kaufmann('02', 'temperature')

        self.T = T
 
    def set_acceleration(self, function):
        """
            sets the desired acceleration function where function 
            must be a function with variables t and x. This will
            be supplied to leapfrog integrator.
        """
        self.acceleration_function = function
    
    def write_out_wind(self, filename = None):
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
        if filename == None:
            filename = self.name + "_wind.dat"

        ok_to_write = True
        
        if not hasattr(self, 'n'):
           ok_to_write = False
           print "Need to calculate density before writing out"
        if not hasattr(self, 'T'):
           ok_to_write = False
           print "Need to calculate temperature before writing out"


        if ok_to_write:
            header = "#t       r       v        n        rho     T\n"


            f = open(filename,'w')
            f.write(header)
            fmt = "%8.8E %8.8E %8.8E %8.8E %8.8E %8.8E\n"
            for i in np.arange(np.size(self.r)):
                f.write(fmt%( self.t[i],   self.r[i], self.vr[i],
                              self.n[i], self.rho[i], self.T[i]))
            f.close()
        
def spherical_NFW(t, xyz, c=12, M200=1.0E12*cgs.Msun , rho_crit = 9.74E-30):
    """
    acceleration function for the spherical NFW profile
    """


    R200 = (3.0 * M200 / (4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
#    fc   = np.log(1.0 + c) - c/(1.0+c)    
    rscale = R200 / c

#    print R200 / cgs.kpc
#    print rscale / cgs.kpc

    r = np.sqrt(np.sum(xyz**2, axis=-1))
    val = np.log( r/rscale + 1.0)/r**3.0 - 1.0/((rscale*r**2)*(r/rscale + 1.0))

    return -1.0* val[:,np.newaxis] * xyz * cgs.G * M200

def leapfrog_integrate(acceleration_func, x0, v0, dt, nsteps, 
                       verbose=True, t1=0.0,args=()):
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



    if verbose:
        print "Entering integration loop"
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

        if i % 1E5 == 0.0: # print out progress
            print "t = %4.2e, i = %4.4i"%(t,i)

        print "Exiting integration loop. Finishing integration"
    else: # to avoid looping over if statements for no reason
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

def gal_to_cartesian(l, b, d, xyz_sun=np.array([8.0,0.0,0.0])*cgs.kpc):
    """
        Converts galactic longitude and latitude to 
        galactocentric cartesian coordinates

    Parameters
    ----------
    l : float
        Galactic longitude in degrees
    b : float
        Galactic latitude in degrees
    d : float
        Distance to object from sun in cm
    xyz_sun : numpy array
        Coordinates of sun in galactocentric cartesian
    
    Returns
    -------
    xyz : numpy array
        3D numpy array containing x,y,z galactocentric cartisain coords
    """

    l *= np.pi / 180.0
    b = (b)*(np.pi / 180.0)
    
    # make life easier by taking sins and cosines and saving
    cosl = np.cos(l)# * np.sign(l)
    cosb = np.cos(b)# * np.sign(b)
    sinb = np.sin(b)
    sinl = np.sin(l)

    # convert to heliocentric cartesian coordinates
    x = (d * cosb * cosl)
    y = (d * cosb * sinl)
    z = (d * sinb       )
 
    xyz = np.array([x,y,z])
    # convert to galactocentric
    xyz += xyz_sun



    return xyz


def convert_proper_motion(l, b, mu_l, mu_b, d, rv,
                          lsr_vel = np.array([-10.0,5.25,7.17])*cgs.km,
                          vc      = 237.0 * cgs.km):
    """
        Function converts proper motion measurements to 
        galactocentric cartesian velocities.

        Parameters
        ----------
        l : float
            galactic longitude in degrees
        b : float
            galactic latitude in degrees
        mu_l : float
            proper motion in l in degrees per second
        mu_b : float
            proper motion in b in degrees per second
        d    : float
            heliocentric distance to object in cm
        rv   : float
            heliocentric radial velocity in cm/s
        lsr_vel : numpy array, optional
            LSR velocity of sun in cm/s. Default 11.1,12.24,7.25 km/s
        vc : float, optional
            circular velocity at position of sun in cm/s. Default 212 km/s
    
        Returns
        --------
        v_xyz : numpy array
            3D array containing galactocentric cartesian velocities in 
            cm/s       
    """

    l *= np.pi / 180.0
    b = (b)*np.pi/180.0
#    b = (90.0-b)*np.pi/180.0
    mu_l = mu_l * np.pi / 180.0
    mu_b = mu_b * np.pi / 180.0

    # save sines and cosines for convenience
    cosl = np.cos(l)# * np.sign(l)
    cosb = np.cos(b)# * np.sign(b)
    sinl = np.sin(l)
    sinb = np.sin(b)

    # find the heliocentric cartesian velocities
    vx = cosb*cosl*rv + d*cosb*sinl*mu_l + d*cosl*sinb*mu_b
    vy = cosb*sinl*rv - d*cosb*cosl*mu_l + d*sinl*cosb*mu_b
    vz = sinb*rv      - d*cosb*mu_b


    #vx = cosl * sinb * rv - (d*sinl*sinb*mu_l) + (d*cosl*cosb*mu_b)
    #vy = sinl * sinb * rv + (d*sinb*cosl*mu_l) + (d*sinl*cosb*mu_b)
    #vz =        cosb * rv + (d*sinb*mu_b)

    

    # now convert from heliocentric to galactocentric
    v_xyz = np.array([vx,vy,vz])

    print 'bfore change', v_xyz /cgs.km
    v_xyz = v_xyz + lsr_vel
    v_xyz[1] = v_xyz[1] + vc # add circular velocity in y

    
    return v_xyz


class orbit_params:

    known_orbits = {'Draco':{
                          'l':86.3730,         # deg
                          'b':34.7088,         # deg
                          'd':82.4*cgs.kpc,    # +/- 5.8 kpc
                          'e': 0.29,           # ellipticity
                          'rv': -293.3*cgs.km, # +/- 1.0
                          'references': ['Pryor et. al. 2015'],
                          'mu_l' : -23.1 * cgs.mas / (100.0*cgs.yr), # +/-6.3
                          'mu_b' : -16.3 * cgs.mas / (100.0*cgs.yr),  # +/-6.3
                          'v_cyl' : np.array([27.0,89.0,-212.0])*cgs.km,
                          'v_r'   : -98.5 *cgs.km, # +/- 2.6,
                          'v_t'   : 210.0 *cgs.km, # +/- 25 
                          }
                    }

    def __init__(self, name = None):
        if name == None:
            print "please choose a orbit from the following list"
            print "use known_orbits.set_name function to choose"
            self.list_known_orbits()

        else:
            
            self.set_name(name)
    
    def set_name(self, name):
        self.name = name

        # save orbital parameter dict
        self.params = orbit_params.known_orbits[name]

        # give all of the important parameters their
        # own variable
        self.l    = self.params['l']
        self.b    = self.params['b']
        self.d    = self.params['d']
        self.rv   = self.params['rv']
        self.refs = self.params['references']
        self.mu_l = self.params['mu_l']
        self.mu_b = self.params['mu_b']

        # get the galactocentric cartesian coordinates
        self.xyz = gal_to_cartesian(self.l, self.b, self.d)

        # get the galactocentric cartesian velocities
        self.v_xyz = convert_proper_motion(self.l, self.b, 
                                           self.mu_l, self.mu_b, 
                                           self.d, self.rv)


        

    def list_known_orbits(self):
        print orbit_params.known_orbits.keys()

    
    



