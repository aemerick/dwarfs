from __future__ import division

import numpy as np

# for numerical derivation and integration
from scipy      import integrate    # for quad integration
from scipy      import optimize     # for root finder


# constants in cgs units and unit conversions
import cgs as cgs



class DF:

    def __init__(self, dprof):
        """
        Initialize the DF class by passing a dark matter density profile object
        
        Parameters
        ----------
        dprof : general_dm_profile object
            Dark matter profile object from dm_density_profiles.py
        """
        
        self.dprof = dprof
        
    def load_df(self, filename):
        """
        Loads f(E) from a file. File assumed to have header names of "E" and "f"
        
        Parameters
        ----------
        filename : string
            Path of file to read from
        """

        try:
            data = np.genfromtxt(filename, names=True)
        
            self.E = data['E']
            self.f = data['f']
            
            npoints = np.size(self.f)
            print "%4i DF points successfuly loaded from "%(npoints) + filename
            
        except:
            print "Error in loading from file. Check header names. Should be E and f"
        
        
        
        return
        
        
    def compute(self, n_points, verbose = True, E = None, filename = None):
        """
        Tabulates the distribution function. DF is computed as given in BT XXXX.
        Energy values and f(E) are stored as parameters of the object.
        
        Parameters
        ----------
        n_points : int
            Number of points to compute for the distribution function. The energy range
            is automatically computed over the minimum / maximum energy values possible as determined
            by the 'small_r' and 'large_r' parameters from the dm density class. This is a numerically
            tractable way of integrating from Psi = 0 to Psi = E, as the former occurs as r -> inf and
            latter only at r = 0
        verbose : boolean, optional
            Turns on/off verbose output. If off, nothing is printed. Default : True
        E : ndarray, optional
            Energy values to compute DF (i.e. f(E) ). This overrides the automatic calculation of E. 
            Be careful with this. Default : None
        filename : string, optional
            If filename is provided, writes out E and f(E) to a text file. Default : None
            
         Returns;
         f : ndarray
            Distribution function as computed at energy values E. 
        """

        # if E is not specified, choose it 
        #
        if E == None:
            E_max = self.relative_potential(self.dprof.small_r)
            E_min = self.relative_potential(0.9*self.dprof.large_r)

            E = np.logspace( np.log10(E_min), np.log10(E_max), n_points)

           
        # integrand in the DF integral
        def _integrand(x, E_val):

            r_val = self._r_root_finder(x)

            if (E_val - x) == 0.0:
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
            if verbose:
                print "%03i Computing value for E = %.3E"%(i,E_value),
            self.f[i] = integrate.quad(_integrand, lower_bound, E_value, args=(E_value,))[0] #+ f_prev


            self.f[i] = self.f[i] + 1.0/np.sqrt(E_value) * (self.dprof.first_derivative(0.99*self.dprof.large_r))*(1.0/(self._dPsi_dr(0.99*self.dprof.large_r)))

            if verbose:
                print " - f = %0.3E"%(self.f[i])

            i = i + 1

        # multiply by the constants and normalize
        normalization = 1.0/self.dprof.M_sys

        self.f = self.f * normalization / (np.sqrt(8.0) * np.pi*np.pi)

        self.E = E

        # write out if filename is given
        if (not filename == None):
            f = open(df_filename,'w')
            f.write("# E f\n")
            for i in np.arange(np.size(self.f)):
                f.write("%.8E %.8E\n"%(self.E[i],self.f[i]))
            f.close()
        
        
        if scalar_input:
            return np.squeeze(self.f)
        else:
            return self.f



        
    def _d2rho_dPsi2(self,r):
        """
        Computes the second derivative of density with respect to relative potential using the chain rule
        """

        dr_dPsi = 1.0 / self._dPsi_dr(r)

        result = self.dprof.second_derivative(r) * (dr_dPsi)**2.0
       
        result = result - self.dprof.first_derivative(r) * self._d2Psi_dr2(r) * (dr_dPsi)**3.0

        return result

 
    def _r_root_finder(self, psi_value):

        """ 
        Solves the equation phi_value - dprof.phi(r) = 0.0 for r 
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
        Derivative of relative potential with respect to r... Again, Psi = - phi
        """
        return -1.0 * self.dprof.dPhi_dr(r)


    def _d2Psi_dr2(self, r):
        """
        Second derivative of relative potential with respect to r... Again, Psi = - phi
        """
        return -1.0 * self.dprof.d2Phi_dr2(r)





def hernquist_df(E, M, a):
    """
    Computes the Hernquist distribution function analytically. This is used for testing purposes. Should 
    get errors of a few percent with at most 8-10% closer to the maximum energy value (i.e. r -> 0 in the 
    potential).
    
    Paramters
    ---------
    E : ndarray
        Energy values to evaluate f(E)
    M : float
        System mass
    a : float
        System scale radius
        
    Returns
    ---------
    f : ndarray
        f(E), as evaluated at E
    """

    # scale E
    E = E * a / (cgs.G * M)

    F = (np.sqrt(2.0) * (2.0*np.pi)**3 * (cgs.G * M * a)**(3.0/2.0))**(-1.0)
 
    F = F * np.sqrt(E) / (1.0 - E)**2
 
    F = F * ( (1.0 - 2.0*E)*(8.0*E*E - 8.0*E - 3.0) + (3.0 * np.arcsin(np.sqrt(E)))/(np.sqrt(E*(1.0-E))))

    return F




