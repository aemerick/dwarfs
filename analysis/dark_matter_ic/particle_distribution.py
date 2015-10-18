"""
particle_distribution

                description: Particle distribution class. Generates a distribution of dark matter
                             particles following a provided distribution function generated from a desired
                             density profile. Following BT and Kazantizidis et. al. 2005
                created:                 10.14.15
                author:                 Andrew J. Emerick, Columbia University ; AMNH
                contact:                 emerick@astro.columbia.edu
"""

from __future__ import division

import numpy as np
from scipy      import optimize     # for root finder
from scipy      import interpolate  # for optimizing generating the PD
from scipy      import integrate

import cgs as cgs

import multiprocessing
from multiprocessing import Pool
import os

OUTPUT = multiprocessing.Queue()

def _while_loop(pd, nmax, max_loop):
    """
    Due to the way the parallel multiprocessing works, this must be defined at the top level
    ... there may be a better way to do this....
    """
    
    n_particles = 0
    loop_counter = 0
    F_max = np.max(pd.DF.f)

    if pd.optimize:
        relative_potential = pd._interpolate_relative_potential
    else:
        relative_potential = pd.DF.relative_potential
              
            
    pos = np.zeros((nmax, 3))
    vel  = np.zeros((nmax, 3))   
        
    while (( n_particles < nmax) and (loop_counter < max_loop)):
            
        r = pd._choose_position()
        Psi = relative_potential(r)    
   
        v = pd._choose_velocity(r, Psi)
     
        E = Psi + 0.5 * v * v
        
        f_E = pd.DF.interpolate_f(E)
            
        F = np.random.rand() * F_max
            
        if F <= f_E: # accept particle
            index = n_particles 
                
            # convert position to cartesian using random theta and phi
            theta = np.random.rand() * np.pi
            phi   = np.random.rand() * 2.0 * np.pi
                
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
                
            pos[index] = r * np.array([x,y,z])
                
            # repeat for velocity using new random numbersw
            theta = np.random.rand() * np.pi
            phi   = np.random.rand() * 2.0 * np.pi
               
            vx = np.sin(theta) * np.cos(phi)
            vy = np.sin(theta) * np.sin(phi)
            vz = np.cos(theta)
                
            vel[index] = v * np.array([vx,vy,vz])
                
             
            n_particles = n_particles + 1
                
        else:
              
            continue

            loop_counter = loop_counter + 1
    
    OUTPUT.put([pos,vel])
    return pos, vel
                



class particle_distribution:
    
    def __init__(self, DF = None, N = None, M = None, optimize = False):
        """
        Initialize the particle_distribution class by providing a distribution function object,
        and desired number of particles (or particle mass).
        
        Parameters:
        DF  : distribution_function object
            Distribution function object as defined in df.py. Must already have the DF evaluated
            and interpolateable, or loaded from a file and interpolateable
        N   : integer
            Desired number of particles. Particle mass is computed accordingly from the total
            system mass         
        M   : float, optional
            Desired particle mass (in cgs) can be provided instead of particle number. Default None
        optimize : logical, optional
            Turn on optimization? Default is False. Optimze on operates substantially faster for all
            particle counts. Gains improve for larger particle count. See explanation below:
            
            By default, the cumulative mass 
            function and potential of the halo are computed exactly, which requires quite a bit 
            of numerical integration. This means things run pretty slow, even if using the
            parallelized method. Initial tests give roughly 3-5 hours for 10k particles on 4 cores. 
            This isn't bad, but 1E6 particles would be unfeasable. Optimize True intead tabulates 
            cumulative mass and potential functions at the start in log R and interpolates
            along them when they need to be evaluated. This should work fine, but does require
            care. 
        """
        
        self.DF = DF
        
        
        if (N == None and M == None and DF == None):
            self.M_part = None
            self.N_part = None
            self.small_r = None
            return

        
        if (not M == None):
            N = int(self.DF.dprof.M_sys) / M
            self.M_part = M
        elif (not N == None):
            self.M_part = self.DF.dprof.M_sys / (1.0 * N)
        
            
        self.N_part   = N
        self.optimize = optimize
        self.small_r  = self.DF.dprof.small_r
        
        if self.optimize:
            _my_print("Optimization turned on. Sanity check your results")
            self._optimize_npoints = 1.0E3
            
            # tabulate profiles for interpolation
            self._tabulate_cumulative_mass()
            self._tabulate_relative_potential()
            
            # set small r as smallest r in the interpolated profiles
            self.small_r = np.min(np.array([self._cumulative_mass_r,self._relative_potential_r]))
            self.small_r = 10.0**self.small_r
            
        else:
            _my_print("Optimization turned off. Will run slow, even if parallelized")
        
      
        
    def generate_particle_distribution(self, max_loop = np.inf, outfile=None):
        """
        Given a distribution function (df) object, uses the accept-reject method to populate the halo
        with N particles. Requires that object is initialized with either desired number of 
        particles or particle mass.
        
        Paramters
        ----------
        max_loop : integer, optional
            Maximum number of loops over accept-reject to make before quitting. This is mainly
            for testing, to prevent an infinte loop. Default is np.inf
            
        outfile  : string, optional
            Filepath to write out particle distribution. If nothing is supplied, the result is 
            not written to file. Write out can be done at any time by calling the 
            particle_distribution.write_pd function. Default None
            
        Returns
        -------
            Nothing. Sets IC's as particle_distribution.pos and particle_distribution.vel 
        """
    
        self.pos   = np.zeros((self.N_part, 3))
        self.vel   = np.zeros((self.N_part, 3))
         
        

        F_max = np.max(self.DF.f)
        

        n_particles = 0
        loop_counter = 0
        
        if self.optimize:
            relative_potential = self._interpolate_relative_potential
        else:
            relative_potential = self.DF.relative_potential
            
        
        
        # Continue until max number of particles chosen, or until max loop counter
        while ((n_particles < self.N_part) and (loop_counter < max_loop)):
            
            # choose random position, eval potential, choose velocity
            r     = self._choose_position()
            
            Psi   = relative_potential(r)    
            v     = self._choose_velocity(r, Psi)
        
            E     = Psi + 0.5 * v * v
        
            # interpolate along DF to find f(E) of chosen particle
            f_E = self.DF.interpolate_f(E)
            
            # random number from 0 to F_max for accept reject
            F = np.random.rand() * F_max
            
            if F <= f_E: # accept particle

                
                # convert position to cartesian using random theta and phi
                theta = np.random.rand() * np.pi
                phi   = np.random.rand() * 2.0 * np.pi
                
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                
                # save particle position
                self.pos[n_particles] = r * np.array([x,y,z])
                
                # repeat for velocity using new random numbers
                theta = np.random.rand() * np.pi
                phi   = np.random.rand() * 2.0 * np.pi
                
                vx = np.sin(theta) * np.cos(phi)
                vy = np.sin(theta) * np.sin(phi)
                vz = np.cos(theta)
                
                # save particle velocity
                self.vel[n_particles] = v * np.array([vx,vy,vz])
                
            
                n_particles = n_particles + 1
                
                            
            if (loop_counter % 50) == 0:
                _my_print("Have %4i particles. On loop %6i"%(n_particles, loop_counter))
            loop_counter = loop_counter + 1
                
                
        if (not outfile == None):
            self.write_pd(outfile)
                
        return self.pos, self.vel

    def parallel_generate_particle_distribution(self, max_loop = np.inf, Ncore = 1, outfile=None):
        """
        Given a distribution function (df) object, uses the accept-reject method to populate the halo
        with N particles. Alternatively, one can specify the desired particle mass instead.
        """
    
        self.pos   = np.zeros((self.N_part, 3))
        self.vel   = np.zeros((self.N_part, 3))
                         
              
        # start running
        nmax = self.N_part / Ncore
        #pool = Pool(processes = Ncore)
        #pool.apply_async(_while_loop,)
        #result = pool.map(_while_loop, args=(self, nmax, max_loop,))
        #print result.get(timeout = 100)
        #p = Process(target=_while_loop, args=(nmax, max_loop,))
        jobs = []
        for i in np.arange(Ncore):
            p = multiprocessing.Process(target=_while_loop, args=(self, nmax, max_loop,))
            jobs.append(p)
            p.start()
        
        for p in jobs:
            p.join()
           
        results = [OUTPUT.get() for p in jobs]
        
        results = np.array(results)
        
        pos = results[:,0]
        pos = pos.reshape(self.N_part,3)
        self.pos = pos
        
        vel = results[:,1]
        vel = vel.reshape(self.N_part,3)
        self.vel = vel
        
        
        if (not outfile == None):
            self.write_pd(outfile)
                
        return self.pos, self.vel
    
    
    def write_pd(self, outfile):
        
        f = open(outfile, 'w')
        
        header = "#m x y z vx vy vz\n"
        fmt    = "%12.12E %12.12E %12.12E %12.12E %12.12E %12.12E %12.12E\n"
        
        f.write(header)
        for i in np.arange(self.N_part):
            f.write(fmt%(self.M_part, self.pos[i,0], self.pos[i,1], self.pos[i,2], 
                                      self.vel[i,0], self.vel[i,1], self.vel[i,2]))
        
        
        f.close()
        
        return
    

    def _choose_velocity(self, r, Psi):
        
        u = np.random.rand()
        
        v_max = np.sqrt(2.0 * Psi)
        
        return u * v_max
    
    def _choose_position(self):
        """
        Chooses a random position for the particle weighted by the 
        cumulative mass function. This is necessary as the particle
        distribution is not uniform in radius. 
        """
        
        
        # function used to find R given M
        def _root_function(r, func, uval, m_tot):
            
            return uval * m_tot - func(r)
        
        # optimization switches
        if self.optimize:
            mass_func = self._interpolate_cumulative_mass
        else:
            # use exact profile
            mass_func = self.DF.dprof.cumulative_mass
        
        
        umin = mass_func(1.00001*self.DF.dprof.small_r) / self.DF.dprof.M_sys
        umax = mass_func(0.99999*self.DF.dprof.large_r) / self.DF.dprof.M_sys
        
        failed = True
        # root finder may fail occasionally if r is too close to zero
        # keep drawing random number until it works.
        # alternate soln would be to draw from M(small_r)/M_tot to
        # M(large_r) / M_tot instead of 0 to 1... 
        #while failed:
        u = np.random.rand()*(umax - umin) + umin
            
        #    try:
        r = optimize.brentq(_root_function, self.small_r, self.DF.dprof.large_r, 
                                    args=(mass_func ,u,self.DF.dprof.M_sys,))
        #        failed = False
                
        #    except: 
        #        failed = True

        
        
        return r

    def _interpolate_cumulative_mass(self, r):
        """
        At runtime, tabulates cumulative mass function over log r and interpolates along it
        """
        
        # interpolate
        #spline = interpolate.UnivariateSpline(self._cumulative_mass_r,
        #                                      self._cumulative_mass_m)
        
        # linear interpolation is more reliable, assuming number of points
        # is large enough
        spline = interpolate.interp1d(self._cumulative_mass_r, self._cumulative_mass_m)
        
        return 10.0**spline(np.log10(r))
    
    def _interpolate_relative_potential(self, r):
        """
        At runtime, tabulates cumulative mass function over log r and interpolates along it
        """
        
        # interpolate
        #spline = interpolate.UnivariateSpline(self._relative_potential_r,
        #                                      self._relative_potential_psi, k = 1)
        
        # linear interpolation is more reliable assuming number of points
        # is large enough
        spline = interpolate.interp1d(self._relative_potential_r, self._relative_potential_psi)
        
        return 10.0**spline(np.log10(r))

        
    def _tabulate_cumulative_mass(self):
        
        rmax  = self.DF.dprof.large_r
        rmin  = self.DF.dprof.small_r*1.000001

        
        r   = np.logspace(np.log10(rmin), np.log10(rmax), self._optimize_npoints)
       
        
        # compute cumulative mass
        cumulative_mass = self.DF.dprof.cumulative_mass(r)
       
        # save logged values 
        self._cumulative_mass_r     = np.log10(r)
        self._cumulative_mass_m     = np.log10(cumulative_mass)
                
        _my_print("Completed mass tabulation")
        return
        
    def _tabulate_relative_potential(self):
        
        rmax  = self.DF.dprof.large_r    
        rmin  = self.DF.dprof.small_r*1.000001

                
        r     = np.logspace(np.log10(rmin), np.log10(rmax), self._optimize_npoints)

        # compute relative potential
        relative_potential = self.DF.relative_potential(r)

        # save logged values        
        self._relative_potential_r     = np.log10(r)
        self._relative_potential_psi   = np.log10(relative_potential)
        
        _my_print("Completed potential tabulation")
        return
        

    def load_particle_ic(self, file_name):
        """
        Given a file name, loads particle initial conditions.
        Particles are assumed to have the same mass.
        File must have headers labeled:
          m x y z vx vy vz
        """
   

        data = np.genfromtxt(file_name, names = True)

        self.N_part = np.size(data['x'])

        self.pos = np.array([data['x'], data['y'], data['z']])
        self.pos = self.pos.T.reshape(self.N_part,3)
        self.vel = np.array([data['vx'], data['vy'], data['vz']])
        self.vel = self.vel.T.reshape(self.N_part,3)
        
        self.M_part = data['m'][0] # assuming all particles have same mass

        _my_print('loaded %6i particles from '%(self.N_part) + file_name)
        return

    def r(self):
        """
        Returns radial position from cartesian coordinates
        """

        r = np.sqrt(self.pos[:,0]**2 + self.pos[:,1]**2 + self.pos[:,2]**2)

        return r


    def density_profile(self, nbins):
        """
        Using a binning procedure, computes the density profile from the particle r)
        positions and velocities
        """


        r = self.r()

        r_bins = np.linspace(0.0, np.max(r), nbins + 1)

        # now bin with np hist
        r_hist, r_bins = np.histogram(r, bins = r_bins)

        # calculate the volume in each bin
        volume = 4.0 * np.pi * (r_bins[1:]**3 - r_bins[:-1]**3) / 3.0

        # now calculate the bin centers
        r_cent = 0.5*(r_bins[1:] + r_bins[:-1])

        # number density
        density = self.M_part * r_hist / volume

        return r_cent, density

    def potential_from_particles(self, nbins):
        """
        Uses the binned density profile from the particles to calculate 
        (approximately) what the total gravitational potential should be.
        This should be compared to the exact potential generated by the 
        N-body calculation. With enough particles, they should be the same.
        """

        r, dens = self.density_profile(nbins)

        dens_function = interpolate.UnivariateSpline(r, dens)


        integrand_1 = lambda x : x * x * dens_function(x)
        integrand_2 = lambda x :     x * dens_function(x)
        
        rmin, rmax = np.min(r), np.max(r)
        
        pot = np.zeros(nbins)       
        for i in np.arange(nbins):
            

            A = integrate.quad(integrand_1, rmin,  r[i])[0]
            B = integrate.quad(integrand_2,  r[i], rmax)[0]

            pot[i] = A/r[i] + B

        pot = - 4.0 * np.pi * cgs.G * pot

        return r, pot
    
def _my_print(string):   
    print "[Particle Distribution]: ", string
