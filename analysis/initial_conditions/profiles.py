from __future__ import division

import numpy as np
import cgs as cgs
import matplotlib.pyplot as plt
from ic_generator import find_rm_gatto as find_rm
from scipy import optimize as opt
from scipy import integrate

def NFW_DM(r, r_s=None, c=12., M200=1.0E12*cgs.Msun, rho_crit = 9.74E-30):
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
    """

    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    # scale radius is virial divided by c
    if r_s == None:
        r_s = R200 / c
    else:
        c = R200 / r_s

    # the scale density is given as
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    return rho_s / ( (r/r_s)*(1+r/r_s)**2)
    
def burkert_DM(r, r_s, M200, rho_crit=9.74E-30):
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

    Returns :
    rho : array
        Array or float of dimensions of r. Dark matter density 
        profile
    """
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    rho_o = (R200/r_s)**2 * (200.0/3.0) * rho_crit    

    return rho_o/((1+r/r_s)*(1+(r/r_s)**2))


def Burkert_isothermal_gas(r, r_s, M200, T, n_o, mu=1.31,
                                rho_crit=9.74E-30):
    """
        Returns the gas density profile that is in HSE
        with a Burkert DM potential.

    """
    R200 = (3.0*M200/(4.0*np.pi*200.0*rho_crit))**(1.0/3.0)

    # central dark matter denisty
    rho_DM = (R200/r_s)**2 * (200.0/3.0) * rho_crit    

    rho_o = n_o * cgs.mp * mu

    # constant in exponential from gas properties
    C_gas = mu * cgs.mp / (cgs.kb * T)  
    # constant in exp from DM profile
    D_B = 4.0*np.pi*cgs.G*rho_DM*r_s**2

    print "rhoDM   cgas    dB"
    print rho_DM, C_gas, D_B

    # R is the unitless radisu
    R = r / r_s

    rho = np.exp(-C_gas * D_B *\
          (0.25*np.log(1.0+R**2) + 0.5*np.arctan(R) - 0.25*np.log(1.+R**2)/R -\
           0.50*np.log(1.0+R)/R  -0.5*np.log(1.0+R) + 0.5))
    rho = rho_o * rho
   

# OLD PROFILE FROM INCORRECT DM PROFILE
#    rho = rho_o * np.exp(- C_gas * D_B * (0.0 +\
#                              0.5*np.log(1.0 + (r/r_s)**2) + \
#                              np.arctan(r/r_s)/(r/r_s)))

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


    print "c = ", c, "R200 = ", R200, "r_s = ", r_s 
    print "R200 = ", R200/cgs.kpc, " kpc -- r_s = ", r_s/cgs.pc, " pc"
        
    # scale density for the DM halo
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    # constant in the exponential
    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu *cgs.mp/(cgs.kb*T)

    print "pi, G, rho_s, r_s, mu, mp, kb, T"
    print "params:", np.pi, cgs.G, rho_s, r_s, mu, cgs.mp, cgs.kb, T
    print "C_NFW = ", C_NFW
    # central mass density 
    rho_o = n_o * cgs.mp * mu
    print "rho_o = ", rho_o, "M200 = ", M200
    # gas profile

    rho = np.zeros(np.size(r))
    rho[0]   = rho_o
    rho[r>0] = rho_o * np.exp(-C_NFW * (1.0 - np.log(1.0+r[r>0]/r_s)/(r[r>0]/r_s)))

    RM=    find_rm(rho_o, C_NFW, r_s, Pcorona, cgs.kb/(mu*cgs.mp)*T)
    print "RM = ", RM, RM/cgs.pc

    return rho, RM


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

    print "c = ", c, "R200 = ", R200, "r_s = ", r_s
    print "R200 = ", R200/cgs.kpc, " kpc -- r_s = ", r_s/cgs.pc, " pc"

    # scale density for the DM halo
    rho_s = 200.0/3.0 * rho_crit * c**3 / (np.log(1.0+c) - c/(1.0+c))

    # constant in the exponential
    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu *cgs.mp/(cgs.kb*T)

    print "pi, G, rho_s, r_s, mu, mp, kb, T"
    print "params:", np.pi, cgs.G, rho_s, r_s, mu, cgs.mp, cgs.kb, T
    print "C_NFW = ", C_NFW
    # central mass density 
    n_o = (Pcorona/cgs.kb) / T
    n_o = n_o * np.exp(C_NFW * (1.0 - np.log(1.0+rmatch/r_s)/(rmatch/r_s)))
    rho_o = n_o * cgs.mp * mu
    print "n_o   = ", n_o
    print "rho_o = ", rho_o, "M200 = ", M200
    # gas profile

    rho = np.zeros(np.size(r))
    rho[0]   = rho_o
    rho[r>0] = rho_o * np.exp(-C_NFW * (1.0 - np.log(1.0+r[r>0]/r_s)/(r[r>0]/r_s)))

#    RM=    find_rm(rho_o, C_NFW, r_s, Pcorona, cgs.kb/(mu*cgs.mp)*T)
#    print "RM = ", RM, RM/cgs.pc

    RM = rmatch
    return rho, RM


def solve_burkert(M_DM, r_DM, r_s, M_HI, r_HI, T_dwarf,
                  mu_dwarf=1.31, mu_halo=0.6, T_halo=None, 
                  n_halo = None, rho_crit=9.74E-30):

    """ 
       Given a few properties, solves for the constants
       needed to set up gas in HSE with a burkert potential
    """

    # find central density using notation in Faerman et. al. 2013 
    f_M = 1.5 * (0.5*np.log(1+(r_DM/r_s)**2) +\
                     np.log(1.0+(r_DM/r_s))  -\
                     np.arctan(r_DM/r_s))
    rho_DM = (M_DM / f_M) * (3.0/(4.0*np.pi)) / (r_s**3)

    # solve the burkert density profile at (R200,200*rho_crit)
    eq_solve = lambda x : (200.0*rho_crit/rho_DM)*((1.0+x)*(1.0+x*x)) - 1.0

    # find the rot of the Eq. (x = R200 / r_s)
    R200 = r_s * opt.bisect(eq_solve, 0.0, 10.0)
    M200 = 4.0/3.0 * np.pi * R200**3 * rho_crit * 200.0

    # we now have M200 and r_s, which defines the DM profile
    # now, find n_o / rho_o to define the gas density profile

    C_g = mu * cgs.mp / (cgs.kb * T_dwarf)
    D_B = 4.0 * np.pi * cgs.G * rho_DM * r_s**2

    # integrate density equation to get cumulative mass and solve
    # for the necessary central density.
    def __integrand(r,C_g=C_g,D_B=D_B,r_s=r_s):
        x = r / r_s
 
        if (r>0):
            rho = np.exp(-C_g * D_B *\
              (0.25*np.log(1.0+x**2) + 0.5*np.arctan(x) - 0.25*np.log(1.+x**2)/x -\
              0.50*np.log(1.0+x)/x  -0.5*np.log(1.0+x) + 0.5))

 
#            val = np.exp(-C_g * D_B *(1.0 + 0.5*np.log(1.+x**2)-\
#                                      np.arctan(x)/x))
        else:
            val = 1.0

        val = 4.0*np.pi * r * r *val

        return val

    # numerically integrate density profile for mass to get rho_o
    rho_o = M_HI / integrate.quad(__integrand, 0.0, r_HI)[0]
    n_o = rho_o / (mu_dwarf * cgs.mp)

#    density = lambda r: rho_o * np.exp(-C_g * D_B * (1.0+\
#                                0.5*np.log(1.0+(r/r_s)**2)-\
#                                np.arctan(r/r_s)/(r/r_s)))
    density = lambda r: Burkert_isothermal_gas(r, r_s, M200, T_dwarf,
                                               n_o, mu=mu_dwarf)

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

    return r_s, M200, n_o, T_halo, n_halo
 

def solve_NFW(M_DM, r_DM, r_s, M_HI, r_HI, T, 
              mu=1.31, mu_halo=0.6, T_halo=None,n_halo=None,
              rho_crit = 9.74E-30):
    """
        Given some dark matter mass and a scaling parameter
        r_s, solve for the dark matter profile parameters
    """

    rho_s = M_DM/(4.0*np.pi*r_s**3) / (np.log(1.0+r_DM/r_s) - r_DM/(r_s+r_DM))

    solve_c = lambda c: (200.0/3.0)*c**3/(np.log(1.0+c)-c/(1.0+c)) *\
                        (rho_crit/rho_s) - 1.0

    c =  opt.bisect(solve_c, 1.0, 40.0)

    R200 = c*r_s

    M200 = (4.0*np.pi/3.0)*200.0*rho_crit * R200**3

    # now solve for the central gas density given M_HI at r_HI

    C_NFW = 4.0*np.pi*cgs.G*rho_s*r_s**2 * mu * cgs.mp /(cgs.kb*T)


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
    
    return c, r_s, M200, n_o, T_halo, n_halo

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

