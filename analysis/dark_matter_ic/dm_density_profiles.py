"""
dm_density_profiles

                description: Density profile class for dark matter profiles
                created:                 09.21.15
                author:                  Andrew J. Emerick, Columbia University ; AMNH
                contact:                 emerick@astro.columbia.edu
"""


from __future__ import division

import numpy as np

from scipy.misc import derivative
from scipy import integrate
import cgs as cgs


# actually... just make generalized density function and derivatives using
# the forms from Kazantzidis paper... 
# that way many density profiles can be described just by using 3 functions
# and changing alpha, beta, and gamma

class general_dm_profile:

    def __init__(self, name):
 
        #density_profile.__init__(self, name, self.density)
        self.name    = name

        # set values for large and small r ... this is for integration purposes
        # mostly. Using np.inf for large r does not work always. Must be careful with this one
        # default values for large and small r are below
        # when r_vir is set, large r is automaticall scaled to 10,000 * virial radius
        self.small_r = 1.0E-100 * cgs.pc
        self.large_r = 1.0E4    * cgs.kpc

        self._set_params()
 
    def _set_params(self):
        self.profile_shape_params = []
        self.r_s      = None
        self.r_vir    = None
        self.r_decay  = None
        self.rho_s    = None
        self.M_vir    = None
        self.rho_crit = None

    def set_params(self, profile_shape_params = None, r_s = None, r_vir = None, r_decay = None,
                   rho_s = None, M_vir = None, rho_crit = None):

        if (not profile_shape_params == None):
            if (len(profile_shape_params) == 3):
            
                self.profile_shape_params = profile_shape_params
                self.calculate_epsilon()
            else:
                print 'Profile shape must have 3 parameters (alpha, beta, gamma)'

        if (not r_s      == None):
            self.r_s      = r_s      
            self.small_r  = 1.0E-6 * self.r_s
            self.large_r  = 100.0  * self.r_s;
            self.calculate_epsilon()

        if (not r_vir    == None): 
            self.r_vir = r_vir ; self.calculate_epsilon()

        if (not r_decay  == None): self.r_decay  = r_decay  ; self.calculate_epsilon()
        if (not rho_s    == None): self.rho_s    = rho_s
        if (not M_vir    == None): self.M_vir    = M_vir
        if (not rho_crit == None): self.rho_crit = rho_crit
            
        if ((not r_vir == None or not r_decay == None) and (not self.r_vir == None and not self.r_decay == None)):
            self.large_r = self.r_vir + 3.0 * self.r_decay
            
            
        
        if self._check_params():
            self._calculate_system_mass()
            
    def _check_params(self):
        """
        make sure everything is set
        """
        everything_OK = True
        
        params = [self.M_vir, self.r_vir, self.r_s]
        
        if any(p == None for p in params):
            everything_OK = False 
            
        elif (self.rho_s == None):
            self.calculate_rho_s() ; everything_OK = True
            
        if (not len(self.profile_shape_params) == 3):
            everything_OK = False 
            
        if ((self.profile_shape_params[1] <= 3) and (self.r_decay == None)):
            everything_OK = False
            
        return everything_OK
        
    def _calculate_system_mass(self):
        """
        Calculates total system mass, or mass out to large_r. This shouldn't be more than a few percent larger
        than M_vir (more for systems requireing the exponential density cutoff)
        """
        
        
        self.M_sys = self.cumulative_mass(self.large_r)
        
        return
        

    def calculate_rho_s(self):
        """
        Calculates and sets the central density if M_vir, r_vir, and r_s are known
        """
        alpha,beta,gamma = self.profile_shape_params 
       
        r_s = self.r_s
 
        # need to do some units magic here to make sure things work with the scipy integrate function

        integrand = lambda x : x*x / ((x/r_s)**(gamma) * (1.0 + (x/r_s)**(alpha))**((beta-gamma)/alpha))

        rmin = self.small_r ; rmax = self.r_vir

        
 
        self.rho_s = (self.M_vir/(4.0*np.pi)) * (integrate.quad(integrand, rmin, rmax)[0])**(-1.0)

    def calculate_epsilon(self):
        """
        Calculates the exponent in the exponential decay portion of the density profile
        such that the logarithmic slope at rvir is continious
        """

        if (hasattr(self, 'r_vir') and hasattr(self, 'r_s') and hasattr(self, 'profile_shape_params')\
                                   and hasattr(self, 'r_decay')):

            try:
                c = self.r_vir / self.r_s
                alpha,beta,gamma = self.profile_shape_params

                self.epsilon = (-gamma - beta * c**(alpha)) / (1.0 + c**(alpha))\
                                     + self.r_vir / self.r_decay     
            except TypeError:
                self.epsilon = None
                
             

    def _set_values_check(self):
        if ((self.rho_s == None) and (not self.M_vir == None and not self.r_vir == None)):
            self.calculate_rho_s()
        elif (self.rho_s == None):
            print "Error: Must set both M_vir and r_vir OR rho_s"

        if (not hasattr(self, 'epsilon')):
            self.calculate_epsilon()
        elif (self.epsilon == None):
            self.calculate_epsilon()

    def density(self, r):


        self._set_values_check()
        alpha, beta, gamma = self.profile_shape_params
        
        r = np.asarray(r)
        scalar_input = False
        if r.ndim == 0:
            r = r[None]
            scalar_input = True
        
        rho = np.zeros(np.shape(r))

        # compute the density
        c = (r / self.r_s)
        rho = self.rho_s / ( c**gamma * (1.0 + c**alpha)**((beta-gamma)/alpha) )
   

        # use exponential cutoff if beta <= 3
        if beta <= 3:
            # now calculate for values greater than the virial radius
            c = self.r_vir / self.r_s
  
            rho[r > self.r_vir ] = self.rho_s / ( c**gamma * (1.0 + c**alpha)**((beta-gamma)/alpha))
     
            rho[r > self.r_vir ] = rho[r>self.r_vir] * (r[r>self.r_vir]/self.r_vir)**(self.epsilon) *\
                                                        np.exp(-(r[r>self.r_vir]-self.r_vir)/self.r_decay)
       

        if scalar_input:
            return np.squeeze(rho)
        else:
            return rho
        


    def first_derivative(self, r):
        """
        Computes, analytically, the first derivative of the dark matter density profile 
        with respect to radius. Where available, the coded derivative makes use of the definitions
        for the density profile.
        """
        self._set_values_check()
        alpha, beta, gamma = self.profile_shape_params

        r = np.asarray(r)
        scalar_input = False
        if r.ndim == 0:
            r = r[None]
            scalar_input = True

        first_deriv = np.zeros(np.shape(r))


        
        first_deriv = -1.0 * self.density(r) / r * \
                                    (beta * (r/self.r_s)**alpha + gamma) / ((r/self.r_s)**alpha + 1.0)


        if beta <= 3:
            first_deriv[ r > self.r_vir ] = self.density(r[r>self.r_vir]) *\
                                      ((self.epsilon/r[r>self.r_vir]) - 1.0 / self.r_decay)
        



        if scalar_input:
            return np.squeeze(first_deriv)
        else:
            return first_deriv


    def second_derivative(self, r):
        """
        Computes, analytically, the second derivative of the dark matter density profile 
        with respect to radius. Where available, the coded derivative makes use of the definitions
        for the density profile and first derivative.
        """

        self._set_values_check()
        alpha, beta, gamma = self.profile_shape_params

        r = np.asarray(r)
        scalar_input = False
        if r.ndim == 0:
            r = r[None]
            scalar_input = True


        second_deriv = np.zeros(np.shape(r))


        c = r / self.r_s
        second_deriv = self.first_derivative(r) * self.first_derivative(r) / self.density(r)+\
                                      self.density(r) *\
                            ( (alpha*(c)**alpha * (beta*(c)**alpha + gamma)) +\
                              ((c)**alpha + 1.0)*(-alpha*beta*(c)**alpha + beta*c**alpha + gamma)) /\
                              (r * r * ((c)**alpha + 1.0)**2 )


        if beta <= 3:
            second_deriv[r>self.r_vir] = self.first_derivative(r[r>self.r_vir])*\
                                         (self.epsilon / r[r>self.r_vir] - 1.0/self.r_decay) -\
                                         self.density(r[r>self.r_vir])*self.epsilon / (r[r>self.r_vir])**2



        if scalar_input:
            return np.squeeze(second_deriv)
        else:
            return second_deriv

    def cumulative_mass(self, r):
        """
        Uses the defined density function to compute the cumulative mass interior to some radius r.
        """ 
        self._set_values_check()
        alpha, beta, gamma = self.profile_shape_params

        r = np.asarray(r)
        scalar_input = False
        if r.ndim == 0:
            r = r[None]
            scalar_input = True
   
        mass = np.zeros(np.shape(r))
        integrand = lambda x : x * x * self.density(x)

        prev_mass = 0.0; rlow = self.small_r
        for i in np.arange(np.size(r)):
            mass[i] = integrate.quad(integrand, rlow, r[i])[0] + prev_mass
            prev_mass = 1.*mass[i] ; rlow = 1.*r[i]

        mass = mass * 4.0 * np.pi

        if scalar_input:
            return np.squeeze(mass)
        else:
            return mass

    def potential(self, r):

        """
        Uses the defined density function to compute the cumulative mass interior to some radius r.

        potential defined as :
             pot = -4*pi*G [ M(r) / r + integral(r*rho(r) dr, r, infinity) ]

             in the below, A = M(r) / r ... B = integral(r*rho(r) * dr, r, infinity) 

             first portion caclulated using defined cumulative mass function
             second portion integrated here with scipy
        """ 
        self._set_values_check()
        alpha, beta, gamma = self.profile_shape_params

        tolerance = self.small_r
        
        r = np.asarray(r)
        scalar_input = False
        if r.ndim == 0:
            r = r[None]
            scalar_input = True
   
        pot = np.zeros(np.shape(r))

        # second integrand
        integrand = lambda x : x * self.density(x)


        # integral at zero can just be written as the below
        pot[r <= tolerance] = -4.0*np.pi*cgs.G * integrate.quad(integrand, self.small_r, self.large_r)[0]

        # this is first
        A = (self.cumulative_mass(r[r > tolerance]) / (4.0 * np.pi)) / r[r > tolerance]
 

        # compute B:
        # this integral behaves very badly if the upper bound is too large... not sure
        # what is happening (difference of small ##'s somewhere?)... but just need it large enough
        # esp as the density distribution eventually truncates... error becomes small after enough
        # virial radii

        B = np.zeros(np.shape(r[r > tolerance]))
        i = 0
        for rval in r[r > tolerance]:
            B[i] = integrate.quad(integrand, rval, self.large_r)[0]
            i = i + 1

        pot[ r > tolerance ] = -4.0 * np.pi * cgs.G * (A + B)

        if scalar_input:
            return np.squeeze(pot)
        else:
            return pot

    def dPhi_dr(self, r):
        """
        Computes derivative of potential with respect to radius
        """

        return cgs.G * self.cumulative_mass(r) / r**2.0

    def d2Phi_dr2(self, r):
        """
        Second derivative of the potential with respect to radius
        """

        return 4.0 * np.pi * cgs.G * self.density(r) - 2.0 * self.dPhi_dr(r) / r
