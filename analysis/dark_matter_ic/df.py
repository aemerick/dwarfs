import numpy as np

from scipy.misc import derivative
from scipy import integrate
import cgs as cgs

# likely not necessary to import
import dm_density_profiles as dm_prof


#
# Phi and Psi are not the same !!!!!!!!!!!!!!!!!
#

class DF:

    def __init__(self, dprof):

        self.dprof = dprof

        

    def dPhi_dr(self, r):
        """
        First derivative of the potential with respect to radius
        """

        return cgs.G * self.dprof.cumulative_mass(r) / r**2.0

    def d2Phi_dr2(self, r):
        """
        Second derivative of the potential with respect to radius
        """

        return 4.0 * np.pi * cgs.G * self.dprof.density(r) - 2.0 * self.dPhi_dr(r) / r

  
    def tabulate_df(self, E):
        """
        Tabulates the distribution function
        """

        
    def _d2rho_dPhi2(self,r):
        """
        Computes the second derivative of density with respect to potential using the chain rule
        """

        dr_dPhi = 1.0 / self.dPhi_dr(r)

        result = self.dprof.second_derivative(r) * (dr_dPhi)**2.0
       
        result = result - self.dprof.first_derivative(r) * self.d2Phi_dr2(r) * (dr_dPhi)**3.0

        return result

 
