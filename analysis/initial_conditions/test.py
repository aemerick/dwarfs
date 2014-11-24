from astropy import units as u
import numpy as np

kb = 1.380658E-16 * u.erg / u.k
mh = 1.6733E-24 * u.g
G  = 6.67259E-8 * u.cm**3 / u.g / (u.second)**2

R_cl  = (300.0 * u.pc).to(u.cm)
rho_1 = mh * 2.0E-4  / (u.cm**3)
rho_2 = 300.0 * rho_1 

T_1   = 1.8E6 * u.k

P_1   = kb * rho_1 * T_1 / mh
P_2_R = P_1

T_2   = P_2_R * mh / kb / rho_2




x = 4.0*np.pi*G*mh / kb / T_2 / 9.0

rho_core = 1.0 / ((1.0/rho_2) - x * R_cl**2)

r_king = ( 1.0 / (x*rho_core)) ** (0.5)


print "R_cl  = %5.4e"%(R_cl.value)
print "rho_1 = %5.4E"%(rho_1.value)
print "rho_2 = %5.4E"%(rho_2.value)
print "P_1 = %5.4E"%(P_1.value)
print "P_2 = %5.4E"%(P_2_R.value)
print "T_1 = %5.4E"%(T_1.value)
print "T_2 = %5.4E"%(T_2.value)
print "rho_core = %5.4e"%(rho_core.value)
print "R_scale = %5.4e"%(r_king.value)


v_wind = (2.0*np.pi*kb*T_2/(G*mh) * R_cl * ( (rho_2)**2 / (rho_1)**(3.0/2.0)) * G**(3.0/2.0) / (np.pi * (6.0)**0.5) ) ** (1.0/3.0)

print v_wind.to(u.km / u.second)

alpha_1 = 15.0 * kb / (8.0*np.pi*G*mh*mh)
alpha_2 = (5.0/2.0) * (15.0/(8.0*np.pi))**0.5 * (kb/G)**(3.0/2.0) / (mh**2)

M_Jeans = alpha_2 * T_2**(3.0/2.0) * (rho_2/mh)**(-0.5)
R_Jeans = (alpha_1 * T_2 / (rho_2/mh))**0.5

M_Jeans = M_Jeans.to(u.g).value / 1.99E33
R_Jeans = R_Jeans.to(u.pc).value

print "Jeans mass and radius - msun and pc"
print "%5.4e    -     %5.4e"%(M_Jeans, R_Jeans)

M_cloud = 2.0*np.pi * kb * T_2 * R_cl / (G*mh)
M_cloud = M_cloud.to(u.g).value / 1.99E33

print "mass - msun and pc"
print "%5.4e          -        %5.4e"%(M_Jeans, R_cl.to(u.pc).value)
