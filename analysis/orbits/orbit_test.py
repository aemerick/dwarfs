import cgs as cgs
import orbit
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


draco = orbit.orbit_params('Draco')
v     = draco.v_xyz / cgs.km
r     = draco.xyz   / cgs.kpc

print 'velocity', v
print 'mag', np.dot(v,v)**0.5
print 'v_r', np.dot(v,v)/np.dot(r,r)


draco_orb = orbit.dwarf_orbit('Draco')
draco_orb.set_x0(draco.xyz)
draco_orb.set_v0(draco.v_xyz)
draco_orb.set_acceleration(orbit.spherical_NFW)
draco_orb.integrate_orbit() # use defaults
draco_orb.calculate_density() # MB13_adj
draco_orb.calculate_temperature() # NFW
draco_orb.write_out_wind()


# plot some things
x = (draco_orb.x[:,:,0] / cgs.kpc).flatten()
y = (draco_orb.x[:,:,1] / cgs.kpc).flatten()
z = (draco_orb.x[:,:,2] / cgs.kpc).flatten()

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
plt.savefig('3d_orbit.png')
plt.close()

####
t = (draco_orb.t / cgs.Myr).flatten()
r = (draco_orb.r / cgs.kpc).flatten()
plt.plot(t,r)
plt.xlabel('t (Myr)')
plt.ylabel('r (kpc)')
plt.savefig('r_t.png')
plt.close()
