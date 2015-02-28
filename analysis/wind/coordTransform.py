#	NEW!!! COORDINATE TRANSFORMATION 
#
#   Here we run everything in the LMC coord system, 
#   except rotated 100 CCW about the Y-axis
#
#		Munier Salem
#		May 2013

import numpy as np
from numpy import cos, sin, pi
import matplotlib.pyplot as plt
#import astro.cgs as cgs
import cgs as cgs

# following from VDM&K 2103 II, Table 1, Column 3
r_LMC = 50.1 
ra_LMC = 78.76 * np.pi/180.0 
dec_LMC = -69.19 * np.pi/180.0

def angle(a,b):
	return np.arccos(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))) * 180.0/np.pi

p      = 229.0  # pos angle (of line of nodes wrt x-ax in LOS frame (cntrclckws)
i      = 34.7   # incl angle (of z axis, rotated clockwise about line of nodes)
didt   = 0. #103# rate of change of i, in deg / Gyr

boxSize = 4.18682385e23 # cm
v_wind = np.array([314.84136 , -29.03 , 53.9510]) # in LOS frame
v_unc  = np.array([-453.8,54.4,-262.2]) # in LOS frame
shockSpeed = 56969941. # cm/s

# switch to radins
p *= pi/180.0
i *= pi/180.0
didt *= pi/180.0/(3.15569e16)

# from pstn and inclntn ang's, generate matrices for LMC/LOS frames
LMC_to_LOS = np.array([[ cos(p) ,  -sin(p)*cos(i) ,  -sin(p)*sin(i) ],
                       [ sin(p) ,  cos(p)*cos(i) ,  cos(p)*sin(i) ],
                       [ 0.0    , -sin(i)        ,  cos(i)        ]])
LOS_to_LMC = np.linalg.inv(LMC_to_LOS)

# generate matrices to go between LOS / Galactic coords (from Roeland)
LOS_to_galactic = np.array([[ 0.11638, -0.98270, -0.14410],
                            [ 0.57132, -0.05244,  0.81905],
                            [-0.81243, -0.17765,  0.55533]])
galactic_to_LOS = np.linalg.inv(LOS_to_galactic)

# transform wind into LMC frame
lmc_wind = np.dot( LOS_to_LMC , v_wind )

# make all comps +'ve by rotating 90 CW about z-axis
phi = -100.0*np.pi/180.0
LMC_to_sim = np.array([[ np.cos(phi),-np.sin(phi),0.0  ],
                      [  np.sin(phi), np.cos(phi),0.0  ],
                      [  0.0,         0.0,        1.0 ]])
sim_to_LMC = np.linalg.inv(LMC_to_sim)
sim_wind = np.dot( LMC_to_sim , lmc_wind )
sim_wind *= 1.0e5 # cm/s
wind_speed = np.linalg.norm(sim_wind) 

# Matrices for going directly between SIM and LOS frames
sim_to_LOS = np.dot( LMC_to_LOS , sim_to_LMC )
LOS_to_sim = np.linalg.inv(sim_to_LOS)

# Disk Angular momentum is now simply the -Z-axis
L = np.array([ 0.0, 0.0, -1.0])

# calculate delay differential, compared to old coords
center = np.array([0.5, 0.5, 0.5]) * boxSize
delay = np.dot(sim_wind,center)/np.linalg.norm(sim_wind) - center[0]
delay /= shockSpeed

"""
# print results
print "L          =",L
print "L_galactic =",L_galactic
print "wind       =",sim_wind
print "wind_speed =",wind_speed
print "angle      =",np.arccos(np.dot(L,sim_wind)/np.linalg.norm(sim_wind)) * 180.0/np.pi
print "delay      =",delay
"""

# some vectors for plotting
X_sim = np.array([1.0,0.0,0.0])
Y_sim = np.array([0.0,1.0,0.0])
Z_sim = np.array([0.0,0.0,1.0])

Z_LOS = np.dot( LMC_to_sim , np.dot( LOS_to_LMC , [0.0,0.0,1.0]))
Y_LOS = np.dot( LMC_to_sim , np.dot( LOS_to_LMC , [0.0,1.0,0.0]))

Z_disk = L
Y_disk = np.cross(Z_disk,np.array([1.0,0.0,0.0]))
Y_disk /= np.linalg.norm(Y_disk)
X_disk = np.cross(Y_disk,Z_disk)


#
# ===== LOS VELOCITY TOOLS
#


# compute CM velocity in LOS angular system
v_sys = -v_unc[2]
v_t = np.linalg.norm(v_unc[:2]) 
THETA_t = np.arctan(v_unc[1]/v_unc[0])+np.pi/2

# compute useful trig terms from LMC position 
cd0,sd0 = np.cos(dec_LMC),np.sin(dec_LMC)
ca0,sa0 = np.cos(ra_LMC),np.sin(ra_LMC)

def angular_coords(ra,dec):
	"""
	Given RA and DEC, computes angular coords from
	LMC center (See van der Marel et al 2001 & 2002)
	
	returns:
	--------
		rho: angular distance from LMC CM to given coords

		phi: position angle of this distance, measured
		     CCW from x-hat axis (-RA axis)
	"""

	# Compute useful trig values.
	cd,sd = np.cos(dec),np.sin(dec)
	ca,sa = np.cos(ra),np.sin(ra)

	# Compute angular distance rho
	rho = np.arccos(cd*cd0*np.cos(ra-ra_LMC)+sd*sd0)
	sr = np.sin(rho) 

	# Compute two possible values of phi ... 
	phi = np.arccos(-cd*np.sin(ra-ra_LMC)/sr)
	phi_alt = 2.0*np.pi - phi

	# Use third EQ to break degeneracy (comes down to sign)
	tmp = (sd*cd0-cd*sd0*np.cos(ra-ra_LMC))/sr
	if tmp*np.sin(phi) > 0: return rho,phi
	return rho,phi_alt

def v_los(rho,PHI):
	"""
		Given angular coords in LMC LOS plane,
		computes necessary velocity correction to
		account for CM motion and precession/nutation
		of the LMC disk plane. Does NOT include
		orbital motion of a circularized disk (we
		wish to keep this)

		See Van Der Marel 2002
	"""
	THETA = p - np.pi/2.0
	tmp  = v_sys*np.cos(rho)
	tmp += v_t*np.sin(rho)*np.cos(PHI-THETA_t)
	tmp += r_LMC*didt*np.sin(rho)*np.sin(PHI-THETA)*cgs.kpc/1.0e5
	return tmp
