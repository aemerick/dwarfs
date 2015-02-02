import cgs as cgs
import orbit
import numpy as np

draco = orbit.orbit_params('Draco')
v     = draco.v_xyz / cgs.km
r     = draco.xyz   / cgs.kpc

print v
print np.dot(v,v)**0.5
print np.dot(v,v)/np.dot(r,r)


draco_orb = orbit.dwarf_orbit('Draco')
draco_orb.set_x0(draco.xyz)
draco_orb.set_v0(draco.v_xyz)
draco_orb.set_acceleration(orbit.spherical_NFW)
draco_orb.integrate_orbit() # use defaults
draco_orb.calculate_density() # MB13_adj
draco_orb.calculate_temperature() # NFW
draco_orb.write_out_wind()
