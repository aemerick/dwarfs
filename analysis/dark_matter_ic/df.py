from __future__ import division

import numpy as np

# for numerical derivation and integration
from scipy.misc import derivative
from scipy import integrate
from scipy import optimize


# constants in cgs units and unit conversions
import cgs as cgs

# likely not necessary to import
import dm_density_profiles as dm_prof


#
# Phi and Psi are not the same !!!!!!!!!!!!!!!!!
#

class DF:

    def __init__(self, dprof):

        self.dprof = dprof

        self.max_E = self.relative_potential(0.0)    

    def compute(self, E, M_DM, N_DM):
        """
        Tabulates the distribution function
        """
        self.M_DM = M_DM ; self.N_DM = N_DM

        # integrand:
        def _integrand(x, E_val):

            r_val = self._r_root_finder(x)

            if x == E_val:
                integrand_value =  2.0 * np.sqrt(E_val) * self._d2rho_dPsi2(r_val)

            else:
                integrand_value =  (1.0 / np.sqrt(E_val - x)) * self._d2rho_dPsi2(r_val)

            return integrand_value


        E = np.asarray(E)
        scalar_input = False
        if E.ndim == 0:
            E = E[None]
            scalar_input = True


        self.f = np.zeros(np.shape(E))
        

        # be careful about the choice of the lower bound of the integral
        lower_bound = self.relative_potential(0.9999*self.dprof.large_r)

        # change this to go from eval to eval

        i = 0; f_prev = 0.0
        for E_value in E:
            print "%03i Computing value for E = %.3E"%(i,E_value),
            self.f[i] = integrate.quad(_integrand, lower_bound, E_value, args=(E_value,))[0] + f_prev

            print " - f = %0.3E"%(self.f[i])
            f_prev = self.f[i] ; lower_bound = E_value
            i = i + 1

        # multiply by the constants out front


        normalization = self.M_DM**(3.0/2.0) / self.N_DM

        self.f = self.f * normalization / (np.sqrt(8.0) * np.pi*np.pi)

        self.E = E

        if scalar_input:
            return np.squeeze(self.f)
        else:
            return self.f



        
    def _d2rho_dPsi2(self,r):
        """
        Computes the second derivative of density with respect to potential using the chain rule

        NEED TO FIX FOR THE WHOLE RELATIVE POTENTIAL CRAP... make sure the sign is correct

        """

        dr_dPsi = 1.0 / self._dPsi_dr(r)

        result = self.dprof.second_derivative(r) * (dr_dPsi)**2.0
       
        result = result - self.dprof.first_derivative(r) * self._d2Psi_dr2(r) * (dr_dPsi)**3.0

        return result

 
    def _r_root_finder(self, psi_value):

        """ 
        Solves the equation   phi_value - dprof.phi(r) = 0.0 for r
        """
        def _eq_to_solve(x, pval):
            return pval - self.relative_potential(x)


        psi_value = np.asarray(psi_value)
        scalar_input = False
        if psi_value.ndim == 0:
            psi_value = psi_value[None]
            scalar_input = True

        r = np.zeros(np.shape(psi_value))


        i = 0
        for p in psi_value:
            try:
                r[i] = optimize.brentq(_eq_to_solve, 0.0, self.dprof.large_r, args=(p,))
            except:
                print "Failing in brentq"
                print psi_value, p
            i = i + 1

        if scalar_input:
            return np.squeeze(r)
        else:
            return r


    def relative_potential(self, r):
        """
        Relative potential (greek letter capital Psi) is defined here as -1.0 times the potential (greek
        eltter lower case phi). Psi written here always refers to relative potential.
        """

        return -1.0 * self.dprof.potential(r)



    def _dPsi_dr(self, r):
        """
        Derivative of relative potential with respect to r... Psi = - phi
        """
        return -1.0 * self.dprof.dPhi_dr(r)


    def _d2Psi_dr2(self, r):
        """
        Second derivative of relative potential with respect to r... Psi = - phi
        """
        return -1.0 * self.dprof.d2Phi_dr2(r)







