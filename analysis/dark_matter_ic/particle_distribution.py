from __future__ import division

import numpy as np
from scipy      import optimize     # for root finder

import cgs as cgs

import multiprocessing
from multiprocessing import Pool
import os

output = multiprocessing.Queue()
def _while_loop(pd, nmax, max_loop):
    n_particles = 0
    loop_counter = 0
    F_max = np.max(pd.DF.f)

    
            
    pos = np.zeros((nmax, 3))
    vel  = np.zeros((nmax, 3))   
        
    while (( n_particles < nmax) and (loop_counter < max_loop)):
            
        r = pd._choose_position()
        Psi = pd.DF.relative_potential(r)    
   
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
                
            vel[index] = np.array([vx,vy,vz])
                
             
            n_particles = n_particles + 1
                
        else:
              
            continue

            loop_counter = loop_counter + 1
    
    output.put([pos,vel])
    return pos, vel
                



class particle_distribution:
    
    def __init__(self, DF, N, M = None):
        
        self.DF = DF
        
        if (not M == None):
            N = int(self.DF.dprof.M_sys) / M
            self.M_part = M
        else:
            self.M_part = self.DF.dprof.M_sys

            
        self.N_part = N
        
        # calculate some parameters important for DF
        #_initialize_parameters()
        
        
    def generate_particle_distribution(self, max_loop = np.inf, outfile=None):
        """
        Given a distribution function (df) object, uses the accept-reject method to populate the halo
        with N particles. Alternatively, one can specify the desired particle mass instead.
        """
    
        self.pos   = np.zeros((self.N_part, 3))
        self.vel   = np.zeros((self.N_part, 3))
         
        # need to program the accept reject technique

        F_max = np.max(self.DF.f)
        print "%.4E"%(F_max)

        n_particles = 0
        loop_counter = 0
        while ((n_particles < self.N_part) and (loop_counter < max_loop)):
            
            r = self._choose_position()
            Psi = self.DF.relative_potential(r)    
   
            v = self._choose_velocity(r, Psi)
        
            E = Psi + 0.5 * v * v
        
            f_E = self.DF.interpolate_f(E)
            
            F = np.random.rand() * F_max
            
            if F <= f_E: # accept particle

                
                # convert position to cartesian using random theta and phi
                theta = np.random.rand() * np.pi
                phi   = np.random.rand() * 2.0 * np.pi
                
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                
                self.pos[n_particles] = r * np.array([x,y,z])
                
                # repeat for velocity using new random numbersw
                theta = np.random.rand() * np.pi
                phi   = np.random.rand() * 2.0 * np.pi
                
                vx = np.sin(theta) * np.cos(phi)
                vy = np.sin(theta) * np.sin(phi)
                vz = np.cos(theta)
                
                self.vel[n_particles] = v * np.array([vx,vy,vz])
                
             #   print 'succesfully made particle number ', n_particles
                n_particles = n_particles + 1
                
            else:
               # print 'failed with F = %0.4E and F_E = %0.4E'%(F, f_E)
	       # print 'r = %.4E, PSI = %.4E'%(r/cgs.kpc,Psi)
                #print 'v %0.3E'%(v / cgs.km)
                continue

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
           
        results = [output.get() for p in jobs]
        
        results = np.array(results)
        print np.shape(results)
        
        print self.N_part
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
    
#    def load_pd(self, infile):
#        data = np.genfromtxt(infline, names=True)
#        
#        self.pos
    
    def _choose_velocity(self, r, Psi):
        
        u = np.random.rand()
        
        v_max = np.sqrt(2.0 * Psi)
        
        return u * v_max
    
    def _choose_position(self):
    
        # pick a random number 0 to 1
        
        u = np.random.rand()
        
        def _root_function(r, func, uval, m_tot):
            
            return uval * m_tot - func(r)
        
        r = optimize.brentq(_root_function, self.DF.dprof.small_r, self.DF.dprof.large_r, 
                                 args=(self.DF.dprof.cumulative_mass,u,self.DF.dprof.M_sys,))
        
        return r

    def load_particle_ic(self, file_name):

        data = np.genfromtxt(file_name, names = True)

        self.N_part = np.size(data['x'])

        self.pos = np.array([data['x'], data['y'], data['z']])
        self.pos = self.pos.T.reshape(self.N,3)
        self.vel = np.array([data['vx'], data['vy'], data['vz']])
        self.vel = self.vel.T.reshape(self.N,3)

        print 'loaded %6i data points from '%(self.N) + file_name
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
        density = r_hist / volume

        return r_cent, density

