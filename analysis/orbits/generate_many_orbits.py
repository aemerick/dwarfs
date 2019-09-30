import numpy as np
from dwarfs.analysis import cgs
from astropy import units as u

from dwarfs.analysis.orbits import parabolic

#
# Generate many model orbits for all galaxies, inclinations,
# and pericenters
#
R_vir  = 77.77 * u.kpc # viriral radius in NFW LMC model
M_prim = 5.66931E10 * u.Msun # primary virial mass

r_peri   = np.array([0.5, 0.25, 0.125]) * R_vir  # pericenter passages
a        = np.array([0, 30, 60, 90])             # inclination angles
M_second = np.array([1.0E10, 3.12614E8]) * u.Msun # seconary mass
names    = ['SMC','LeoT']

# aw yess, nested loops in python
for M,name in zip(M_second,names):
    for r in r_peri:
        for angle in a:
            orbit = parabolic.parabolic(M_v = M, M_p = M_prim,
                                        b   = r, t_o = -2000.0 * u.Myr,
                                        alpha = angle,
                                        dt  = 2.0 * u.Myr,
                                        t_final = None)

            outname = name + '_r%.2f_a%i_orbit.dat'%(r.value,angle)
            orbit.save_orbit(outname)
            print outname + " Completed"

            del(orbit)




