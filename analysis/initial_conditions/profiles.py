from __future__ import division

import numpy as np
import cgs as cgs
import matplotlib.pyplot as plt
from ic_generator import find_rm_gatto as find_rm
from scipy import optimize as opt
from scipy import integrate

def soliton_density(r, M_h, m_s = 1.5E-22, z = 0):
    """
    Yu et. al. 2014

    r   : radius in cm
    M_h : halo mass in g
    m_s : particle mass in eV
    z   : redshift

    """
    a = 1.0 / (1.0 + z)

    Hubble = 67.0E3  / cgs.Mpc

    A,B = 0.4, 0.1 # constants, at Jerry's suggestion -- very rough
    m22 = m_s / 1.0E-22
    M_s = 0.25 * ( (M_h/cgs.Msun) / (m22*4.7E7))**(1.0/3.0) * m22 * 4.7E7
    R_s = 1.6 / m22 * ( (M_h/cgs.Msun) / 1.0E9)**(-1.0/3.0)
    M_s = M_s * cgs.Msun # now in g
    R_s = R_s * cgs.kpc # now in cm

#
#    old formulae from Jerry
#
#    M_s = A * (cgs.h_bar / (cgs.G * m_s)) * (Hubble * cgs.G / M_h)**(1.0/3.0)
#    R_s = B * (cgs.h_bar**2 / (m_s * cgs.G * M_s))
#    rho = (105/(32.0*np.pi)) * (M_s / R_s**3) * ( 1.0 / (1.0 - (r/R_s)**2))**4
#
    rho = 1.9 * (m22 / 0.1)**(-2) * (R_s/cgs.kpc)**(-4) / (1.0 + 9.1E-2 * (r/R_s)**2)**8
    rho_avg = 4.0 * np.pi * M_s / R_s**3 / 3.0 / cgs.Msun * cgs.kpc**3 
    rho = rho * cgs.Msun / cgs.pc**3

    return rho

def soliton_mass(r, M_h):
    """
    Yu et. al. 2014 - integrate over the DM profile
                      this is unfortunate (analytic seemed difficult)
    """

    integrand =  lambda x : 4.0 * np.pi * x * x * soliton_density(x, M_h)

    result = integrate.quad( integrand, 0.0, r)[0]

    return result

def solve_soliton(r_x, M_x):
    """
    Given a (r,M) pair, solves the soliton mass that matches
    the mass at that radius. 
    """
    
    root_solve = lambda x : soliton_mass(r_x, x) - M_x 
    
    return optimize.brentq(root_solve, 1.0E8 * cgs.Msun, 1.0E11 * cgs.Msun)

def NFW_potential(r, r_s = None, c=12., M200 =1.0E12*cgs.Msun,rho_crit = 9.74E-30):
    """
    Evaluates the NFW potential
    """
 
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    # scale radius is virial divided by c
    if r_s == None:
        r_s = R200 / c
    else:
        c = R200 / r_s

    # the scale density is given as
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    c_nfw = 4.0 * np.pi * cgs.G * rho_s * r_s**2
    r = r / r_s
    if np.size(r) > 1:
        phi = np.zeros(np.size(r))
        phi[r==0] = - c_nfw
        phi[r>0 ] = - c_nfw * np.log(1.0 + r[r>0])/r[r>0]

    else:
        if r > 0:
            phi = np.log(1.0+r)/r
        elif r==0:
            phi = 1.0

        phi = phi * (- c_nfw)

    return phi

def Burkert_potential(r, r_s, M200, rho_crit=9.74E-30):
    """
    Evaluates the burkert potential at a point R
    """
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
    R = R200/r_s

    rho_o = 3.0*M200/(4.0*np.pi*r_s**3) / (1.5*(np.log(1+R) + 0.5*np.log(1+R**2)-\
                               np.arctan(R)))

    #print 'burkert rho_o', rho_o
   
    c_b = 4.0 * np.pi * cgs.G * r_s**2 * rho_o
    r = r / r_s

    r = np.asarray(r)
    scalar_input = False
    if r.ndim == 0:
        r = r[None]
        scalar_input = True
    
    zero = 1.0E-13
    phi = np.zeros(np.size(r))
    phi[r <= zero] = - 0.5 * c_b * (np.pi/2.0)

    R = r[r > zero]

    
    phi[r > zero] = - 0.5 * c_b * ( (1.0/R) * (0.5*np.log(1+R*R) + np.log(1+R) - np.arctan(R)) +\
                                       (np.log(1+R) - 0.5*np.log(1+R*R) - np.arctan(R) + np.pi/2.0) )
    
    #if np.size(r) > 0:
    #     phi = np.zeros(np.size(r))
    #    phi[r==0] = 0.5 * c_b * (-np.pi/4.0 - 1.0)

    #     R = r[r>0] 
    #    phi[r>0] = 0.5 * c_b * (0.5 * np.log(1.0+R*R) + np.arctan(R) -\
    #                            0.5*(np.log(1.0+R*R) + 2.0*np.log(1.0+R))/R -\
    #                            np.log(1.0+R) - np.pi/4.0) 

    if scalar_input:
        return np.squeeze(phi)
    else:
        return phi

def NFW_DM(r, r_s=None, c=12., M200=1.0E12*cgs.Msun, rho_crit = 9.74E-30,
              decay=False, r_decay=None):
    """
        Compute the dark matter density profile for a 
        NFW halo. Taking the scale as R200/c

        Parameters
        -----------
        r : array
            Radius from center of halo in cm
        r_s : float, optional
            The scale paramter. If none is provided, it is computed 
            as R200 / c. If r_s is provided, c is ignored. Default = None
        c : float, optional
            Concentration parameter of the DM halo. This is 
            ignored if r_s is provided. Default c=12 (MW)
        M200 : float, optional
            Virial DM mass of halo. Default 1.0E12 Msun (MW)
        rho_crit : float, optional
            Critical density of the universe. This is a function of 
            redshift. Default 9.74E-30 (z = 0)
        decay : bool, optional
            If True, includes exponential decay at r > r_200 from 
            Springel & White 1999. Default False
        r_decay : float, optional
            Exponential decay scale radius. Default is None and set
            to 0.1 R_vir (only used when decay is True)            
    """

    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    # scale radius is virial divided by c
    if r_s == None:
        r_s = R200 / c
    else:
        c = R200 / r_s

    # the scale density is given as
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    density = rho_s / ( (r/r_s)*(1+r/r_s)**2)

    if decay:
        if r_decay == None:
            density[r > R200] = density[r > R200] * decay_function(r[r>R200],r_decay,R200,r_s,'NFW')

    return density

def burkert_mass(r, r_s, M200, rho_crit=9.74E-30):

    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
    R = R200/r_s

    rho_ds = 3.0*M200 / (4.0 * np.pi * r_s**3) / (1.5*(0.5*np.log(1+R*R) + np.log(1+R) - \
                                                          np.arctan(R)))
    M_ds = 4.0 * np.pi / 3.0 * r_s**3 * rho_ds

    x = r / r_s

    return M_ds * 1.5 * (0.5 * np.log(1.0+x*x) + np.log(1.0+x) - np.arctan(x))

def NFW_mass(r, r_s, M200, rho_crit=9.74E-30):
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
    R = R200/r_s

    rho_ds = 3.0 * M200 / (4.0 * np.pi * r_s**3) / (3.0 * (np.log(1.0+R) - R/(1.0+R)))

    M_ds = 4.0 * np.pi / 3.0 * r_s**3 * rho_ds

    x = r / r_s

    return M_ds * 3.0 * (np.log(1.0 + x) - x/(1.0+x))

def burkert_DM(r, r_s, M200, rho_crit=9.74E-30, decay=False, r_decay=None):
    """
    Given the parameters, calculates the Burkert DM density
    profile from Burkert 1995.
 
    Parameters
    ----------
    r  : array
        Radii to evaluate DM density profile. In cm
    r_s : float
        Burkert profile scaling radius. In cm
    M200 : float
        DM mass at R_200. Or, roughly the virial DM mass. In g
    rho_crit : float, optional
        Critical density of the universe in cgs. Default z = 0
        9.74E-30 g /cm^-3
    decay : bool, optional
        If True, includes exponential decay at r > r_200 from
        Springel & White 1999. Default False
    r_decay : float, optional
        Exponential decay scale radius. Default is None and set
        to 0.1 R_vir (only used when decay is True)

    Returns :
    rho : array
        Array or float of dimensions of r. Dark matter density 
        profile
    """
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)
    R = R200/r_s

 #   rho_o = 3.0*M200/(4.0*np.pi*r_s**3) / (1.5*(np.log(1+R) + 0.5* np.log(1+R**2)-\
 #                              np.arctan(R))
    
    rho_o = 3.0*M200 / (4.0 * np.pi * r_s**3) / (1.5*(0.5*np.log(1+R*R) + np.log(1+R) - \
                                                          np.arctan(R)))

    density = rho_o/((1+r/r_s)*(1+(r/r_s)**2))

    print '[Burkert DM]', rho_o, M200, R200

    
    if decay:
        if r_decay == None:
            density[r > R200] = decay_function(r[r>R200],r_decay,R200,r_s,'Burkert')

            
    return density


def decay_function(r, r_decay, r_vir, r_s, potential_type):
    """
    Exponential decay function from Springel & White 1999 used to 
    truncate the dark matter density function so M(r->infin) -> infin.
    
    Parameters:
    -----------
    r       : float, array
        array of radii
    r_decay : float
        Decay scale radius 
    potential_type : string
        Name of potential type. This is important as scaling depends on
        type of potential. Only NFW and Burkert work at the moment
    """

    if potential_type == 'NFW':
        alpha, beta, gamma = 1.0,3.0,1.0
    elif potential_type == 'Burkert':
        alpha, beta, gamma = 1.0,3.0,1.0

    epsilon = 1.0   

    return (r/r_vir)**epsilon * np.exp(-(r-r_vir)/r_decay)


def Burkert_isothermal_gas(r, r_s, M200, T, n_o, mu=1.31,
                                rho_crit=9.74E-30):
    """
        Returns the gas density profile that is in HSE
        with a Burkert DM potential.

    """
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    # central dark matter denisty
    R = R200/r_s

    rho_DM = 3.0* M200/(4.0*np.pi*r_s**3) /(1.5 * (np.log(1+R) + 0.5* np.log(1+R**2)-\
                                np.arctan(R)))
    rho_o = n_o * cgs.mp * mu
    M200 = 4.0*np.pi/3.0  * (200.0*rho_crit) * R200**3

    print '[Bukert isothermal gas]', rho_DM, rho_o, M200, R200
    # constant f_M(in exponential from gas properties
    C_gas = mu * cgs.mp / (cgs.kb * T)  
    # constant in exp from DM profile
    D_B = 4.0*np.pi*cgs.G*rho_DM*r_s**2

    r = np.asarray(r)
    scalar_input = False
    if r.ndim == 0:
        r = r[None]
        scalar_input = True
    
    
    # R is the unitless radisu
    R = r / r_s

    rho_s = rho_o / (np.exp(-C_gas * Burkert_potential(0.0, r_s, M200)))
    
    zero = 1.0E-13
    rho = np.zeros(np.size(R))
    rho[R <= zero] = rho_o

    R = R[R > zero]

    #rho[R > tolerance] = np.exp(-C_gas * D_B *\
    #      (0.25*np.log(1.0+R**2) + 0.5*np.arctan(R[R>0]) - 0.25*np.log(1.+R[R>0]**2)/R[R>0] -\
    #       0.50*np.log(1.0+R)/R  -0.5*np.log(1.0+R[R>0]) + 0.5))

    rho[ R > zero ] = rho_s * np.exp(-C_gas * Burkert_potential( r[R > zero], r_s, M200))

    
   
    if scalar_input:
        return np.squeeze(rho)
    else:
        return rho



def NFW_isothermal_gas(r, r_s=None, c=None, M200=4.0E7*cgs.Msun,
                          T = 1.0E4, n_o = 0.27,  mu = 1.31, rho_crit=9.74E-30,
                          Pcorona=1.0):
    """
        The isothermal density profile for a gas in HSE with a 
        NFW DM halo. The default parameters are Sextans from
        Gatto et. al. 2013
    """
    # calculate the virial radius
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    if r_s == None and c == None:
        print "must give r_s (scale radius) OR c (concentration parameter)"
        return 0.0
    elif r_s == None:
        r_s = R200 / c
    elif c == None:
        c = R200 / r_s


    # scale density for the DM halo
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    # constant in the exponential
    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu *cgs.mp/(cgs.kb*T)

    # central mass density 
    rho_o = n_o * cgs.mp * mu

    # gas profile
    r = np.asarray(r)
    scalar_input = False
    if r.ndim == 0:
        r = r[None]
        scalar_input = True
        
    zero = 1.0E-13
    rho = np.zeros(np.size(r))
    rho[r <= zero]   = rho_o
    rho[r>zero] = rho_o * np.exp(-C_NFW * (1.0 - np.log(1.0+r[r>zero]/r_s)/(r[r>zero]/r_s)))

    if scalar_input:
        return np.squeeze(rho)
    else:
        return rho


def NFW_isothermal_rmatch(r, r_s=None, c=None, M200=4.0E7*cgs.Msun,
                 T = 1.0E4, mu = 1.31, rho_crit=9.74E-30, Pcorona=1.0,
                 rmatch=300.0*cgs.pc):
    """
        The isothermal density profile for a gas in HSE with a
        NFW DM halo. The default parameters are Sextans from
        Gatto et. al. 2013
    """
    # calculate the virial radius
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    if r_s == None and c == None:
        print "must give r_s (scale radius) OR c (concentration parameter)"
        return 0.0
    elif r_s == None:
        r_s = R200 / c
    elif c == None:
        c = R200 / r_s
        
    r = np.asarray(r)
    scalar_input = False
    if r.ndim == 0:
        r = r[None]
        scalar_input = True
        
        
    # scale density for the DM halo
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    # constant in the exponential
    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu *cgs.mp/(cgs.kb*T)

    # central mass density 
    n_o = (Pcorona/cgs.kb) / T
    n_o = n_o * np.exp(C_NFW * (1.0 - np.log(1.0+rmatch/r_s)/(rmatch/r_s)))
    rho_o = n_o * cgs.mp * mu

    # gas profile
    zero = 1.0E-13
    rho = np.zeros(np.size(r))
    rho[r <= zero]   = rho_o
    rho[r > zero ] = rho_o * np.exp(-C_NFW * (1.0 - np.log(1.0+r[r>zero]/r_s)/(r[r>zero]/r_s)))

    RM = rmatch
    
    if scalar_input:
        return np.squeeze(rho), RM
    else:
        return rho, RM


def solve_burkert(M_DM, r_DM, r_s, M_HI, r_HI, T_dwarf,
                  mu_dwarf=1.31, mu_halo=0.6, T_halo=None, 
                  n_halo = None, rho_crit=9.74E-30):

    """ 
       Given a few properties, solves for the constants
       needed to set up gas in HSE with a burkert potential
    """

    # find central density using notation in Faerman et. al. 2013 
    R = r_DM/r_s
    f_M = lambda x : 1.5* (0.5*np.log(1.0 + x**2) +\
                     np.log(1.0 + x)  -\
                     np.arctan(x) )
    
    

    rho_DM = (M_DM / f_M(R)) * (3.0/(4.0*np.pi)) / (r_s**3)
 
    # Solve the profile for radius with average density equal to
    # 200 * rho_crit (i.e. solve for R200) 
    eq_solve = lambda x : (rho_DM)*f_M(x)/(x**3)-200.0*rho_crit

    R200 = r_s * opt.bisect(eq_solve, .1, 10000.0, xtol=1.0E-12)

    M200 = 4.0*np.pi/3.0  * (200.0*rho_crit) * R200**3
    R = R200/r_s

    print '[solve Burkert]', rho_DM, M200, R200

    
    # we now have M200 and r_s, which defines the DM profile
    # now, find n_o / rho_o to define the gas density profile

    C_g = mu_dwarf * cgs.mp / (cgs.kb * T_dwarf)
    D_B = 4.0 * np.pi * cgs.G * rho_DM * r_s**2

    # integrate density equation to get cumulative mass and solve
    # for the necessary central density.
    def __integrand(r,C_g=C_g,D_B=D_B,r_s=r_s):
        x = r / r_s
 
        if (r>0):
            #val = np.exp(-C_g * D_B *\
              #(0.25*np.log(1.0+x**2) + 0.5*np.arctan(x) - 0.25*np.log(1.+x**2)/x -\
              #0.50*np.log(1.0+x)/x  -0.5*np.log(1.0+x) + 0.5))
            val = np.exp(-C_g * Burkert_potential(r, r_s, M200))

        else:
            val = 1.0

        val = 4.0*np.pi * r * r *val

        return val

    # numerically integrate density profile for mass to get the constant
    # A = 
    A = M_HI / integrate.quad(__integrand, 0.0, r_HI)[0]
    rho_o = A * np.exp( - C_g * Burkert_potential(0.0, r_s, M200))
    
    print ' ===== [solve burkert]', rho_o / np.exp(- C_g * Burkert_potential(0.0, r_s, M200)), rho_o / np.exp(0.25 * np.pi * C_g * D_B)
    
    n_o = rho_o / (mu_dwarf * cgs.mp)

    density = lambda r: Burkert_isothermal_gas(r, r_s, M200, T_dwarf,
                                               n_o, mu=mu_dwarf)

    # now find the pressure equilibrium temperature and density of halo
    n_gas_edge = density(r_HI)/(cgs.mp*mu_dwarf)

    # find the pressure equillibrium
    if (T_halo == None and n_halo == None):
        print 'no halo paramaters provided for pressure eq'
        return c,r_s,M200,n_o

    elif (T_halo == None):
        T_halo = n_gas_edge/n_halo * T_dwarf
    else:
        n_halo     = n_gas_edge * T_dwarf / T_halo

    return r_s, M200, n_o, T_halo, n_halo
 

def solve_NFW_DM(M_DM, r_DM, r_s, rho_crit = 9.74E-30):
    """
        Solves the NFW DM profile for M200, R200, and rho_o given
        the DM mass interior to some radius, and the desired scale
        radius OR the desired concentration parameter
    """

    rho_s = M_DM/(4.0*np.pi*r_s**3) / (np.log(1.0+r_DM/r_s) - r_DM/(r_s+r_DM))

    solve_c = lambda c: (200.0/3.0)*c**3/(np.log(1.0+c)-c/(1.0+c)) *\
                        (rho_crit/rho_s) - 1.0

    c =  opt.bisect(solve_c, 1.0, 40.0)

    R200 = c*r_s

    M200 = (4.0*np.pi/3.0)*200.0*rho_crit * R200**3


    return M200, R200, rho_s, r_s, c


def solve_NFW(M_DM, r_DM, r_s, M_HI, r_HI, T, 
              mu=1.31, mu_halo=0.6, T_halo=None, n_halo=None,
              rho_crit = 9.74E-30, n_o= None, rho_s = None):
    """
        Given some dark matter mass and a scaling parameter
        r_s, solve for the dark matter profile parameters
    """
    if not (rho_s is None):
        func = lambda x : (rho_s * (np.log(1.0 + r_DM/x) - r_DM / (x + r_DM)) * (4.0 * np.pi * x**3) / M_DM) - 1.0
        r_s = opt.bisect(func, 50.0 * cgs.pc, 10.0 * cgs.kpc)
        print 'solved for r_s in NFW', r_s

    M200, R200, rho_s, r_s, c = solve_NFW_DM(M_DM, r_DM, r_s, rho_crit=rho_crit)
    # now solve for the central gas density given M_HI at r_HI

    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu * cgs.mp /(cgs.kb*T)

    if n_o == None: # if no central density is provided, it is solved 
                    # for given r_HI and M_HI. 
        # make the integrand to get the cumulative
        # mass function, or M(r) = integrate( __integrand )
        def __integrand(r,C_NFW=C_NFW,r_s=r_s):
            if (r>0):
                val = np.exp(-C_NFW*(1.0-np.log(1.0+r/r_s)/(r/r_s)))
            else:
                val = 1.0
        
            val = 4.0*np.pi * r * r *val        
  
            return val

        # numerically integrate density profile for mass to get rho_o
        rho_o = M_HI / integrate.quad(__integrand, 0.0, r_HI)[0]
    
        # define the density profile function to get pressure Eq.
        rho = lambda r: rho_o * np.exp(-C_NFW*(1.0-np.log(1.0+r/r_s)/(r/r_s)))

        n_o = rho_o / (mu * cgs.mp)

        # now find the pressure equilibrium temperature and density of halo
        n_gas_edge = rho(r_HI)/(cgs.mp*mu)

        # find the pressure equillibrium
        if (T_halo == None and n_halo == None):
            print 'no halo paramaters provided for pressure eq'
            return c,r_s,M200,n_o

        elif (T_halo == None):
            T_halo = n_gas_edge/n_halo * T
        else:
            n_halo     = n_gas_edge * T / T_halo
    
        rmatch = r_HI
    
    elif not n_o == None and not r_HI == None and r_HI > 0:
        rho_o = n_o * cgs.mp * mu_halo
        
        n_dens = lambda r: n_o * np.exp(-C_NFW*(1.0-np.log(1.0+r/r_s)/(r/r_s)))       

        # now solve for the corona temperature by temperature balance at r_HI
        T_halo = n_dens(r_HI) * T / n_halo
        
        #
        rmatch = r_HI

    else: # if n_o is provided and r_HI and M_HI are not
        rho_o = n_o * cgs.mp * mu_halo

        n_dens = lambda r: n_o * np.exp(-C_NFW*(1.0-np.log(1.0+r/r_s)/(r/r_s)))

        # now, solve for the match radius ... ASSUME HALO PRESSURE IS KNOWN
        Pcorona = n_halo * T_halo * cgs.kb

        Pdwarf = lambda r: n_dens(r) * cgs. kb * T

        eq_solve = lambda r : Pdwarf(r) - Pcorona

        # Now find r such that Pdwarf = Pcorona
        rmatch = opt.bisect(eq_solve, 1.0*cgs.pc, 2000.0*cgs.pc)
    
       
    return c, r_s, M200, n_o, T_halo, n_halo, rmatch

def plot_profile(r, profile, filename=None, persist=False,**kwargs):

    function_dict = {'Isothermal NFW': NFW_isothermal_gas}
    
    rho, RM = NFW_isothermal_gas(r, **kwargs)

    fig, ax1 = plt.subplots()

    ax1.plot(r/cgs.pc, rho, label=profile,lw=1.75,color='black')

    ax1.semilogy()
    ax1.set_xlabel(r'r (pc)')
    ax1.set_ylabel(r'$\rho$ (g cm$^{-3}$)')
    
    if filename == None:
        filename = profile + "_gas_profile.png"

    plt.savefig(filename)
 
    if persist:
        return fig, ax1, rho
    else:
        plt.close()

def cumulative_mass(r, rho):
    """
    computes the cumulative mass profile
    """
    
    mass_enclosed = np.zeros(np.size(r))
    
#    integrand = lambda x : x*x*rho[rx]
    i = 0
    for i in np.arange(np.size(r)):
 #       mass_enclosed[i] = integrate.quad(integrand,0,r[i])[0] * 4.0*np.pi
        #select = r[0:i]
        mass_enclosed[i] = 4.0*np.pi*np.trapz(r[0:i]**2 * rho[0:i], r[0:i])
        i = i + 1
    
    return mass_enclosed

def column_density(r, density_function,R = None, f_H = 0.73, f_ion = 0.0, **kwargs):
    """
       Calculates the column density of a radial number density profile
       assuming spherical symmetry. Assumes f_H = 0.73 and ionization fraction
       of 0.0 unless otherwise supplied (i.e. primordial and neutral). If 
       no radius of object is supplied ('R'), it is assumed to be max(r).
    """
    
    if R == None:
        R = np.max(r)
    
    HI_factor = f_H * (1.0 - f_ion)
    
    integrand = lambda x,b : 2.0*x/np.sqrt(x**2-b**2) * density_function(x,**kwargs)

    NHI = np.zeros(np.size(r)); i = 0
    for bval in r:
        
        igrand = lambda x: integrand(x,bval)
 
        NHI[i] = integrate.quad(igrand, 1.0000000001*bval, R)[1]
        i = i + 1
   
    return NHI * HI_factor

def _run_test():
    Tcorona = 1.8E6
    rho_corona = 1.8E-4 * cgs.mp * cgs.mu# gatto with ionized primordial halo
    Pcorona = Tcorona * rho_corona * cgs.kb / (cgs.mp*cgs.mu)


    r = np.linspace(0.0,10**3.5,10000.0)*cgs.pc
    fig, ax1, rho = plot_profile(r,'Isothermal NFW',persist=True,c=21.5,Pcorona=Pcorona)
    print r[0], r[-1]


    #Tcorona = 1.8E6
    #rho_corona = 1.8E-4 * cgs.mp * cgs.mu# gatto with ionized primordial halo
    #Pcorona = Tcorona * rho_corona * cgs.kb / (cgs.mp*cgs.mu)

    print 'rho ambient', rho_corona

    ax2 = ax1.twinx()
    mu_dwarf = 1.31
    ax2.plot(ax1.get_xlim(),[Pcorona,Pcorona],label='Corona',color='red',lw=1.75,ls='--')
    P_dwarf = cgs.kb * rho/(cgs.mp*mu_dwarf) * 1.0E4

    ax2.plot(r/cgs.pc, P_dwarf,label='Dwarf',color='red',lw=1.75,ls='-')
    ax2.semilogy()
    ax2.set_ylabel('Pressure')
    for t1 in ax2.get_yticklabels():
        t1.set_color('red')
    ax2.legend(loc='best')
    plt.savefig('rho_P.png')
    plt.close()
    


