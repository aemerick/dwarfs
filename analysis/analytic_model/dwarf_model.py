# this will contain a bunch of stuff to
# solve the anaylytical model we are developing

import numpy as np
import cgs as cgs

# given an orbit of a dwarf galaxy as a function of time
#    - Given: r (position), v_dwarf, n_halo, T_halo 
#    - return r_dwarf and M_gas as a function of that
#             time including: RPS and tidal stripping
#
#
# zeroth order approximation is to calculate the RPS time
# scale at each timestep, and do a linear interpolation of 
# the mass loss i.e.
#  M_i+1 = M_i * (t_i+1 - t_i)/timescale
#
# with a little work, the orbital input can time evolve with
# the RPS and tidal stripping mass loss....
#
#
# 

def evolve_satellite():



    return r, m

def rps_timescale(rho_halo, rho_o, v_dwarf, r_o,
                   n_halo=None, n_dwarf=None,
                   mu_halo=0.61, mu_dwarf=1.30):
    """
        Given the central gas mass density, halo gas mass density,
        dwarf velocity, computes timescale of RPS
    """

    if n_halo is not None:
        rho_halo = n_halo * cgs.mp * mu_halo
    if n_dwarf is not None:
        rho_o    = n_dwarf * cgs.mp * mu_halo

    v_shock = 4.0/3.0 * np.sqrt(rho_halo / rho_o) * v_dwarf

    tau = 2.0 * r_o / v_shock

    return tau

def KH_timescale(M_gas, r_o, rho_halo, v_dwarf, n_halo=None, mu_halo=0.61):
    """
    KH timescale for stripping
    """
    if n_halo is not None:
        rho_halo = n_halo *cgs.mp * mu_halo

    M_rate = np.pi * r_o**2 * rho_halo * v_dwarf

    return M_gas / M_rate
