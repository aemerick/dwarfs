from __future__ import division

import numpy as np
from scipy      import optimize     # for root finder

import cgs as cgs


class particle_distribution:
    
    def __init__(self, DF, N, M = None):
        
        self.DF = DF
        
        if (not M == None):
            N = int(self.DF.dprof.M_sys) / M
            self.M_part = M
            
        self.N_part = N
        
        
        # calculate some parameters important for DF
        #_initialize_parameters()
        
        
    def generate_particle_distribution(self):
        """
        Given a distribution function (df) object, uses the accept-reject method to populate the halo
        with N particles. Alternatively, one can specify the desired particle mass instead.
        """
    
        self.pos   = np.zeros((self.N_part, 3))
        self.vel   = np.zeros((self.N_part, 3))
         
        # need to program the accept reject technique

        F_max = np.max(self.DF.f)
    
        n_particles = 0
        while n_particles < self.N_part:
            
            r = _choose_position()
            Psi = self.DF.relative_potential(r)    
   
            v = _choose_velocity(r, Psi)
        
            E = Psi - 0.5 * v * v
        
            f_E = self.DF.interpolate_f(E)
            
            F = np.random.rand() * F_max
            
            if F <= f_E: # accept particle
                n_particles = n_particles + 1
                
                # convert position to cartesian using random theta and phi
                theta = np.random.rand() * np.pi
                phi   = np.random.rand() * 2.0 * np.pi
                
                x = np.sin(theta) * cos(phi)
                y = np.sin(theta) * sin(phi)
                z = np.cos(theta)
                
                self.pos[n_particles] = r * np.array([x,y,z])
                
                # repeat for velocity using new random numbers
                theta = np.random.rand() * np.pi
                phi   = np.random.rand() * 2.0 * np.pi
                
                vx = np.sin(theta) * cos(phi)
                vy = np.sin(theta) * sin(phi)
                vz = np.cos(theta)
                
                self.vel[n_particles] = np.array([vx,vy,vz])
                
            else:
                continue
                
                
        return self.pos, self.vel
        
    
    def _choose_velocity(self, r, Psi):
        
        u = np.random.rand()
        
        v_max = np.sqrt(2.0 * Psi)
        
        return u * v_max
    
    def _choose_position(self):
    
        # pick a random number 0 to 1
        
        u = np.random.rand()
        
        def _root_function(r, func, uval, m_tot):
            
            return uval - m_tot - func(r)
        
        r = optimize.brentq(_root_function, self.DF.dprof.small_r, self.DF.dprof.large_r, 
                                 args=(self.DF.dprof.cumulative_mass,u,self.DF.dprof.M_sys,))
        
        return r