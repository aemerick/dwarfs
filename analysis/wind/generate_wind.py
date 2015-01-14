"""
	grenerate wind

	From Gurtina's LMC orbit data (Besla et al 2012) generates
	a text file that enzo reads to produce the orbit wind.
"""
import numpy as np
import matplotlib.pyplot as plt
import cgs as cgs
import coordTransform as ct
import galaxy as gal

rampTemp = True
temp_enhance = 4.0  # quadruple it
start_time = 12.8 # gyr

# load and extract data
vecs     = np.loadtxt(open('data/LMC_1_5e12_L18.txt','rb'))
time     = (vecs[::-1,0] - vecs[-1,0])/np.linalg.norm(vecs[-1,0])*13.8
position = np.array( [ vecs[::-1,1], vecs[::-1,2] , vecs[::-1,3] ])
velocity = np.array( [ vecs[::-1,4], vecs[::-1,5] , vecs[::-1,6] ])
distance = (np.sum(position**2,axis=0))**.5

wind_LOS = np.dot( ct.galactic_to_LOS, -velocity )
wind_sim = np.dot( ct.LOS_to_sim , wind_LOS )

time = (time - start_time)
run  = time > 0.0
time = time[run]*cgs.gyr
time[0] = 0.0

# calculate temperature
temperature = gal.halo_gas_temperature(distance*cgs.kpc)

if rampTemp:
# double temperature at start of run, pulling down to real value
# over 200 Myr period
	temperature = temperature[run]
	ramp_time    = .2*cgs.gyr # 200 Myr
	ramp = -(temp_enhance/ramp_time)*time + temp_enhance
	ramp[ramp < 1.0] = 1.0
	temperature *= ramp 

# calculate density
halo_kwargs = { 'n0':.46,'rc':.35*cgs.kpc,'beta':.71 }
density = gal.halo_gas_density(distance[run]*cgs.kpc,**halo_kwargs)

# Dump data to file for enzo to read
data = np.column_stack((time,density,temperature,wind_sim[:,run].transpose()*cgs.km))
header = "NumberOfSteps = " + str(len(time))  \
	+ "\n#time                    Density                  Temp" \
	+ "                     Velocity (x,y,z)"
np.savetxt('ICMinflow_data.in',data,header=header,comments='')
