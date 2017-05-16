import numpy as np
import cgs 

def stellar_disk_density(r,z,M=2.7e9*cgs.Msun,a=1.7*cgs.kpc,b=.34*cgs.kpc):
        """
                Cylindrically symmetric density profile of stellar disk, a Plum$
                Kuzmin potential as derived in Miyamoto & Nagai 1975.

                Parameters
                ----------
                r : numpy array or ndarray
                        Cylindrical radius from galactic center
                z : numpy array or ndarray
                        Vertical height above disk plane
                a : float, optional
                        radial scale length
                b : float, optional
                        vertical scale height, default .34 kpc
                M : float, mass

        """
        x = (z**2 + b**2)**.5
        rho = a*r**2 + (a+3.0*x)*(a+x)**2
        rho /= ((r**2+(a+x)**2)**(5.0/2.0)*x**3)
        return b**2*M*rho/(4.0*np.pi)


def stellar_surface_density(r,b=.34*cgs.kpc,**kwargs):
        """
                Calculates stellar disk surface density as a function
                of cylindrical radius, r. See Miyamoto & Nagai 1975.

                Parameters
                ----------
                r : numpy array or ndarray
                        Cylindrical radius from galactic center
                a : float, optional
                        radial scale length
                b : float, optional
                        vertical scale height, default .34 kpc
                M : float, mass
        """
        z = np.linspace(-50,50,5000)*b
        R,Z = np.meshgrid(r,z)
        rho = stellar_disk_density(R,Z,b=b,**kwargs)
        return np.sum(rho,axis=0)*(z[1]-z[0])


def gas_disk_density(r,z,a=1.7*cgs.kpc,b=.34*cgs.kpc,M=5.0e8*cgs.Msun,
                 smooth_radius=7.692*cgs.kpc,trunc_radius=10.0*cgs.kpc):
        """
                Computes gas density at given radius and height above plane,
                profile as in Tonnesen & Bryan 2013, defaults for LMC
                disk (see Roediger 2012)

                Parameters
                ----------
                r : numpy array or ndarray
                        cylindrical radius from galactic center (cgs)
                z : numpy array or ndarray
                        height above plane (cgs)
                a : float, optional
                        radial scale height (cgs)
                        default : 1.7 kpc
                b : float, optional
                        vertical scale height (cgs)
                        default : .34 kpc
                M : float, optional
                        total mass (sort of), cgs
                        default : 5e8 Msun
                smooth_radius : float, optional
                        radius where density begins to get cutoff by cosine fcn, cgs
                        default : 7.692 kpc
                trunc_radius : float, optional
                        radius where density is cutoff, cgs
                        default : 10.0 kpc
        """

        # Compute the disk gas density.
        rho = 1.0/(np.cosh(r/a)*np.cosh(z/b))

        rho *= M/(8.0*np.pi*a**2*b)

        # Apply cosine smoothing at disk edge.
        smooth_length = trunc_radius - smooth_radius
        if smooth_length <= 0.0:
                raise ValueError('truncation radius must be larger than smooth radius')
        if type(r) is float or type(r) is np.float64:
                r = np.array([r])
        cut = np.abs(r) > smooth_radius
        rho[cut] *= 0.5*(1.0+np.cos(np.pi*(np.abs(r[cut])-smooth_radius)/(smooth_length)))
        rho[ np.abs(r) > trunc_radius ] = 0.0

        return rho

def gas_surface_density(r,a=1.7*cgs.kpc,b=.34*cgs.kpc,M=5.0e8*cgs.Msun):
        """
                Computes gas surface density at given radius, profile as in
                Tonnesen & Bryan 2013, defaults for LMC disk (see Roediger 2012)

                Parameters
                ----------
                r : numpy array or ndarray
                        cylindrical radius from galactic center (cgs)
                a : float, optional
                        radial scale height (cgs)
                        default : 1.7 kpc
                b : float, optional, unused
                        vertical scale height (cgs)
                        default : .34 kpc
                M : float, optional
                        total mass (sort of), cgs
                        default : 5e8 Msun
        """
        return M/(8.0*a**2*np.cosh(r/a))


