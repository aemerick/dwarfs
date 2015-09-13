import numpy as np

from scipy.misc import derivative
from scipy import integrate
import cgs as cgs

class density_profile:

    def __init__(self, name, density_function, *density_args):

        self.name = name
        self._density_args = density_args

    def _define_density_function(self, density_function):

       self.density = lambda x : density_function(x, self._density_args)


    def first_derivative(self, r):

        return derivative(self.denisty, r, 1, args=self._density_args)

    def second_derivative(self, r):
        
        return derivative(self.density, r, 2, args=self._density_args)

    def cumulative_mass(self, r):
    
        integrand = lambda x : x * x * self.density(x)


# actually... just make generalized density function and derivatives using
# the forms from Kazantzidis paper... 
# that way many density profiles can be described just by using 3 functions
# and changing alpha, beta, and gamma

class general_dm_profile:

    def __init__(self, name):
 
        #density_profile.__init__(self, name, self.density)
        self.name    = name
        self.small_r = 1.0E-6 * cgs.pc

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

        if (not r_s      == None): self.r_s      = r_s      ; self.calculate_epsilon()
        if (not r_vir    == None): self.r_vir    = r_vir    ; self.calculate_epsilon()
        if (not r_decay  == None): self.r_decay  = r_decay  ; self.calculate_epsilon()
        if (not rho_s    == None): self.rho_s    = rho_s
        if (not M_vir    == None): self.M_vir    = M_vir
        if (not rho_crit == None): self.rho_crit = rho_crit

    def calculate_rho_s(self):
        """
        Calculates and sets the central density if M_vir, r_vir, and r_s are known
        """
        alpha,beta,gamma = self.profile_shape_params 
       
        r_s = self.r_s
 
        # need to do some units magic here to make sure things work with the scipy integrate function

        integrand = lambda x : x*x / ((x/r_s)**(gamma) * (1.0 + (x/r_s)**(alpha))**((beta-gamma)/alpha))

        rmin = self.small_r ; rmax = self.r_vir

        
 
        self.rho_s = (self.M_vir/4.0*np.pi) * (integrate.quad(integrand, rmin, rmax)[0])**(-1.0)
        
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

        # calculate for values less than the virial radius
        c = (r[r <= self.r_vir] / self.r_s)
        rho[r <= self.r_vir] = self.rho_s / ( c**gamma * (1.0 + c**alpha)**((beta-gamma)/alpha) )
 
        
      
        # now calculate it for values greater than the virial radius
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


        rtemp = r [ r <= self.r_vir]
        first_deriv[r <= self.r_vir] = -1.0 * self.density(rtemp)/rtemp * \
                                    (beta * (rtemp/self.r_s)**alpha + gamma) / ((rtemp/self.r_s)**alpha + 1.0)


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


        rtemp = r [ r <= self.r_vir]
        c = rtemp / self.r_s
        second_deriv[r<=self.r_vir] = self.first_derivative(rtemp) * self.first_derivative(rtemp) / self.density(rtemp)+\
                                      self.density(rtemp) *\
                            ( (alpha*(c)**alpha * (beta*(c)**alpha + gamma)) +\
                              ((c)**alpha + 1.0)*(-alpha*beta*(c)**alpha + beta*c**alpha + gamma)) /\
                              (rtemp*rtemp * ((c)**alpha + 1.0)**2 )


        second_deriv[r>self.r_vir] = self.first_derivative(r[r>self.r_vir])*\
                                     (self.epsilon / r[r>self.r_vir] - 1.0/self.r_decay) -\
                                     self.density(r[r>self.r_vir])*self.epsilon / (r[r>self.r_vir])**2



        if scalar_input:
            return np.squeeze(second_deriv)
        else:
            return second_deriv

    def cumulative_mass(self, r):
        self._set_values_check()
        alpha, beta, gamma = self.profile_shape_params

        r = np.asarray(r)
        scalar_input = False
        if r.ndim == 0:
            r = r[None]
            scalar_input = True
   



        mass = np.zeros(np.shape(r))
        integrand = lambda x : 4.0 * np.pi * x * x * self.density(x)

        prev_mass = 0.0; rlow = self.small_r
        for i in np.arange(np.size(r)):
            mass[i] = integrate.quad(integrand, rlow, r[i])[0] + prev_mass
            prev_mass = mass[i] ; rlow = r[i]

        if scalar_input:
            return np.squeeze(mass)
        else:
            return mass
