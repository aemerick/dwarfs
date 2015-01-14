"""
galaxy

		description: Galaxy simulation components for plotting / analysis purposes 
		created:		 09.05.2013
		author:			Munier A. Salem, Columbia University
		contact:		 msalem@astro.columbia.edu
"""

import numpy as np
import cgs as cgs

def NFW_mass(r,c=12,M200=1.0e12*cgs.Msun,rho_crit=9.74e-30):
	"""
		Spherical NFW DM profile 
		
		Parameters
		----------
		r : numpy array or ndarray
			spherical radius from galactic center
		c : float, optional
			NFW concentration parameter
			default: 12
		M200 : float, optional
			NFW Halo mass (cgs)
			default : 10^12 Msun
		rho_crit : float, optional
			Critical density (cgs)
			default : 9.74e-30 g/cm^3
	"""
	R200 = ( 3.0*M200/( 4.0*np.pi*200.0*rho_crit ) )**(1.0/3.0)
	fc = np.log(1.0+c) - c/(1.0+c)
	x = r*c/R200
	return M200/fc*( np.log(1.0+x) - x/(1.0+x) )


def burkert_mass(r,r0=3.0*cgs.kpc,d0=3.4e-24):
	"""
		Spherical DM profile of mass enclosed from 
		Mori & Burkert 2000. The density formula 
		corresponding to this M(r) is:

			rho(r) = rho_0 r_0^3  / [ (r+r0)*(r^2+r0^2)]
    
		Parameters
		----------
		r : numpy array or ndarray
			spherical radius from galactic center
		r0 : core radius   
	"""
	x = r/r0
	return np.pi*d0*r0**3* \
		(-2.0*np.arctan(x)+2.0*np.log(1+x)+np.log(1+x**2))


def halo_gas_density(r,n0=(.46+.74),rc=(.35+.29)*cgs.kpc ,beta=(.71-.14)):
	"""
		Beta gas density profile for galactic halo in hydrostatic equilibrium
		See Miller & Bregman 2013
		
		Parameters
		----------
		r : numpy array or ndarray
			spherical radius from galactic center (cgs)
		n0 : float, optional
			core number density (cgs)
			default : .46+.74
		rc : float, optional
			core radius (cgs)
			default : .35 + .29 kpc
		beta : float, optional
			density falloff exponent
			default : .71 - .14

	"""
	return cgs.mp*cgs.mu*n0*(1.0 + (r/rc)**2)**(-3.0*beta/2.0)


def halo_gas_temperature(r,M=None,gamma=1.5,**kwargs):
	"""
		Gas temperature for galactic halo in hydro equ., assumes NFW
		Indpendent of gas density (ignores this)

		Parameters
		----------
		r : numpy array or ndarray
			spherical radius from galactic center, cgs
		M : numpy array or ndarray, optional
			Mass enclosed at a given radius, if none provided, 
			NFW mass enclosed is assumed, with default params
			described below
			default : None
		gamma : float
			Numerical parameter of model (See Miller & Bregman 2013)
			default: 1.5
		c : float, optional
			NFW concentration parameter
			default: 12
		M200 : float, optional
			NFW Halo mass (cgs)
			default : 10^12 Msun
		rho_crit : float, optional
			Critical density (cgs)
			default : 9.74e-30 g/cm^3
	"""
	if M is None:
		M = NFW_mass(r,**kwargs)
	return gamma*cgs.G*cgs.mu*cgs.mp*M/(3.0*r*cgs.kb)


def lmc_halo_gas_density(r,r0=2.0*cgs.kpc,d0=1.8e-27,**kwargs):
	"""
		LMC halo gas density, assuming individual gas particles obey
		KE = 1/2 PE.

		Parameters
		----------
		r : numpy array or ndarray
			spherical radius from galactic center, cgs
		r0 : float, optional 
			core radius, cgs
			default : 2.0 kpc
		d0 : float, optional
			core density, cgs
			default : 1.8e-27 g/cm^3
		kwargs : optional 
			keyword args for temperature profile
			default : NFW w/ MW params (see halo_gas_temperature)
	"""
	T = halo_gas_temperature(r,**kwargs)
	T0 = halo_gas_temperature(r0,**kwargs)
	density = d0*(T0/T)*(r/r0)**(-3)
	density[ density > d0 ] = d0
	return density


def halo_gas_mass(**kwargs):
	"""
		Mass of Beta profile galactic halo
		See Miller & Bregman 2013

		Parameters
		----------
		n0 : float, optional
			core number density (cgs)
			default : .46+.74
		rc : float, optional
			core radius (cgs)
			default : .35 + .29 kpc
		beta : float, optional
			density falloff exponent
			default : .71 - .14

	"""
	r = np.linspace(1,200,1000)*cgs.kpc
	density = halo_gas_density(r,**kwargs)
	
	dr = r[1:]-r[:-1]
	density = density[1:]
	r = r[1:]
	dV = 4.0 * np.pi * r**2 * dr
	
	return np.dot(dV,density)/cgs.Msun


def stellar_disk_density(r,z,M=2.7e9*cgs.Msun,a=1.7*cgs.kpc,b=.34*cgs.kpc):
	"""
		Cylindrically symmetric density profile of stellar disk, a Plummer-
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


def disk_restore_pressure(r,stellar_kwargs={},gas_kwargs={}):
	"""
		Computes the gravitational 'pressure' felt by the disk gas
		due to the disk stars using

			( 2 pi G ) Sigma_gas Sigma_stars

		Typically used to compare to the Ram Pressure experienced 
		by the disk 

		Parameters
		----------
		r : numpy array or ndarray
			cylindrical radius from galactic center (cgs)
		stellar_kwargs : dictionary 
			model parameters for stellar disk (optional)
			default : empty
		gas_kwargs : dictionary
			model parameters for gaseous disk (optional)
			default : empty
	"""
	return 2.0*np.pi*cgs.G*gas_surface_density(r,**gas_kwargs) \
		*stellar_surface_density(r,**stellar_kwargs)
