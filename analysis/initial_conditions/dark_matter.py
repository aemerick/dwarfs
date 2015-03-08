import profiles as prof
import numpy as np

# this will contain some python test code for setting up IC's of
# live dark matter particles... following work by
# Vijayaraghavan & Ricker 2015 - Kazantzidis et. al. 2004

def particle_positions(particle_mass, total_mass, r_vir, density_function=prof.NFW_DM,
                                      **kwargs):
    """
    Given mass per particle, total halo mass, and density type, calculates initial
    positions for number of partilces.
    """

    num_particles = total_mass / particle_mass
    if (num_particles != floor(num_particles)):
        print "[Dark Matter]: Halo mass not evenly divisible by particle mass, flooring"
        num_particles = floor(num_particles)

    u_array = np.random.random(num_particles)

    integrand = lambda r: density_function(r, **kwargs)*r*r

    rmin = 1.0E-6 * cgs.pc
    rmax = np.inf

    eq_solve  = lambda r: integrate.quad(integrand(r),rmin,r) / integrate.quad(integrand(r),rmin,rmax)
    
    r = np.zeros(num_particles)
    for u in u_array:
        r[i] = opt.bisect(eq_solve,rmin,r_vir)
        i = i + 1

    return r
