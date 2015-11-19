"""
df

                description: Distribution function class to generate a DF given a DM
                             density profile. Following BT and Kazantizidis et. al. 2005
                created:                 09.21.15
                author:                 Andrew J. Emerick, Columbia University ; AMNH
                contact:                 emerick@astro.columbia.edu
"""


from __future__ import division

import numpy as np

# for numerical derivation and integration
from scipy      import integrate    # for quad integration
from scipy      import optimize as opt     # for root finder
from scipy      import interpolate

# constants in cgs units and unit conversions
import cgs as cgs

LR_FACTOR = 1.0

class DF:

    def __init__(self, dprof, optimize = False, optimize_npoints = 1.0E3):
        """
        Initialize the DF class by passing a dark matter density profile object
        
        Parameters
        ----------
        dprof : general_dm_profile object
            Dark matter profile object from dm_density_profiles.py
        """
        
        self.dprof = dprof
        self.optimize = optimize
        self.dprof.set_optimization( self.optimize )
        
        if self.optimize:
            self._optimize_npoints = optimize_npoints
            self._tabulate_relative_potential()
        
        
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
            _my_print("%4i DF points successfuly loaded from "%(npoints) + filename)
            
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
            E_min = self.relative_potential(LR_FACTOR*self.dprof.large_r)
            #E_min = 0.0
            
            E = np.zeros( n_points + np.int(0.1*n_points))
            E[:n_points] = np.logspace( np.log10(E_min), np.log10(E_max*0.9), n_points)

            E[n_points:] = np.logspace( np.log10(E_max*0.9), np.log10(E_max), np.int(n_points*0.1))
      #  print np.min(E), np.max(E), E[0]
        #   E = np.linspace(E_min, E_max, n_points)
        # integrand in the DF integral
        def _integrand(x, E_val):

            r_val = self._r_root_finder(x)

            if np.abs((E_val - x)) <= 1.1E-14:
                integrand_value =  2.0 * np.sqrt(E_val) * self._d2rho_dPsi2(r_val)

            else:
                integrand_value =  (1.0 / np.sqrt(E_val - x)) * self._d2rho_dPsi2(r_val)

            return integrand_value

        def _recast_integrand(u, E_val):
            
            x = E_val * np.sin(u)**2.0
            r_val = self._r_root_finder(x)
            
            integrand_value = 2.0 * np.sqrt(E_val) * np.sin(u) * self._d2rho_dPsi2(r_val)

            return integrand_value
        
        
        
        E = np.asarray(E)
        scalar_input = False
        if E.ndim == 0:
            E = E[None]
            scalar_input = True


        self.f = np.zeros(np.shape(E))
        

        # be careful about the choice of the lower bound of the integral
        #lower_bound = self.relative_potential(LR_FACTOR*self.dprof.large_r)
        lower_bound = 0.0
        #lower_bound = self.relative_potential(LR_FACTOR*self.dprof.large_r)
        
        # change this to go from eval to eval
        normalization = 1.0/self.dprof.M_sys
        normalization =  normalization / (np.sqrt(8.0) * np.pi*np.pi)
        
        # save energy values
        self.E = E

        # if output
        if (not filename == None):
            outfile = open(filename,'w')
            outfile.write("# E f\n")
        
        i = 0; f_prev = 0.0
        for E_value in self.E:
            if verbose:
                print "%03i Computing value for E = %.3E"%(i,E_value),
            #self.f[i] = integrate.quad(_integrand, lower_bound, E_value, args=(E_value,))[0] #+ f_prev
            #
            # if E_min / E_value ~ 1, b/c of machine precision it may be
            # just a bit greater than 1 and np.arcisn will return nan...
            # doing error catching to avoid this
            #
            E_ratio = E_min / E_value
            if (E_ratio - 1.0) > 0 and (E_ratio - 1.0) <= 1.0E-10:
                E_ratio = 1.0
            elif (E_ratio > 1.0):
                _my_print("Fatal Error in Energy range in DF calculation.")
                return
            
            lower_bound = np.arcsin(np.sqrt(E_ratio))

            
            self.f[i] = integrate.quad(_recast_integrand, lower_bound, np.pi/2.0, args=(E_value,))[0]# + f_prev

            self.f[i] = self.f[i] + 1.0/np.sqrt(E_value) * (self.dprof.first_derivative(LR_FACTOR*self.dprof.large_r))*(1.0/(self._dPsi_dr(LR_FACTOR*self.dprof.large_r)))

            self.f[i] = self.f[i] * normalization
            
            if verbose:
                print " - f = %0.3E"%(self.f[i])
            elif (i % 500 ==0):
                print "%03i E = %.3E"%(i,E_value),
                print " - f = %0.3E"%(self.f[i])
                

            if (not filename == None):
                outfile.write("%.8E %.8E\n"%(self.E[i],self.f[i]))
                
                
            i = i + 1

        
        if (not filename == None):
            outfile.close()        
        
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
        if self.optimize:
            relative_potential_function = self._interpolate_relative_potential
        else:
            relative_potential_function = self.relative_potential
        
        def _eq_to_solve(x, pval):
            return pval - relative_potential_function(x)


        psi_value = np.asarray(psi_value)
        scalar_input = False
        if psi_value.ndim == 0:
            psi_value = psi_value[None]
            scalar_input = True

        r = np.zeros(np.shape(psi_value))


        i = 0
        for p in psi_value:
            #try:
            # was 0 to large_r
            # AE 11/17 right now this is failing to work for both hernquist and nfw... bounds are the issue here
            # ........
            a = 0.99999
            b = 1.00001
           # print '-----------'
           # print p, self.relative_potential(self.dprof.small_r*a), self.relative_potential(self.dprof.large_r*b), self.relative_potential(self.dprof.large_r)
           # print 'asdfasdfasdfasdf'
           # print p, relative_potential_function(self.dprof.small_r*a), relative_potential_function(self.dprof.large_r*b), relative_potential_function(self.dprof.large_r)
          #  print p - self.relative_potential(self.dprof.small_r*a)
          #  print p - self.relative_potential(self.dprof.large_r*b)
          #  print '-----'
            r[i] = opt.brentq(_eq_to_solve, self.dprof.small_r*a, self.dprof.large_r*b, args=(p,))
            #except:
            #print "Failing in brentq"
            #print psi_value, p, self.relative_potential(1.0E100*self.dprof.large_r), self.relative_potential(0.0)
            #r[i] = 1.0E100*self.dprof.large_r
            i = i + 1

        if scalar_input:
            return np.squeeze(r)
        else:
            return r


    def interpolate_f(self, E, *args, **kwargs):
        """
        Uses cubic spline interpolation in log f and log E to compute f at an arbitrary location.
        """
       
        log_E = np.log10(self.E)
        log_f = np.log10(self.f)
     
        if E < np.min(self.E):
        #try:
            f = 0.0
        else:
            spline = interpolate.interp1d(log_E,log_f, *args, **kwargs)
            f = spline(np.log10(E))
            f = 10.0**(f)

        
        #except:
        #    spline = interpolate.UnivariateSpline(log_E, log_f, *args, **kwargs)
        #    f = spline(np.log10(E))
        #    _my_print("Switching to Univariate Spline in interpolation")

        
    
        return f
        
    def relative_potential(self, r):
        """
        Relative potential (greek letter capital Psi) is defined here as -1.0 times the potential (greek
        eltter lower case phi). Psi written here always refers to relative potential.
        """
        


        return -1.0 * self.dprof.potential(r)

    def _tabulate_relative_potential(self):
        
        rmax  = self.dprof.large_r * 1.01 # some extra splop
        rmin  = self.dprof.small_r * 0.999 
        sub_sample = 0.2
        
        r = np.zeros( self._optimize_npoints + np.ceil(sub_sample*self._optimize_npoints))
        
        r[:self._optimize_npoints] = np.logspace(np.log10(rmin),
                                                np.log10(rmax*0.98), self._optimize_npoints)
        
        r[self._optimize_npoints:] = np.logspace(np.log10(rmax*0.98),
                                                 np.log10(rmax),
                                               np.ceil(sub_sample*self._optimize_npoints))

        # compute relative potential
        relative_potential = self.relative_potential(r)

        # save logged values        
        self._relative_potential_r     = np.log10(r)
        self._relative_potential_psi   = np.log10(relative_potential)
        
        _my_print("Completed potential tabulation")
        return
    
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

def _my_print(string):
    
    print "[DF] : ", string
    return

