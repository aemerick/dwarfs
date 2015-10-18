import numpy as np
import cgs as cgs

import df as df
import dm_density_profiles as dm
import particle_distribution as pd

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time


NFW = dm.general_dm_profile('LT_NFW')
NFW.set_params(profile_shape_params=[1.0,3.0,1.0])
NFW.set_params(M_vir = 3.1E8 * cgs.Msun, r_vir = 1370.0 * cgs.pc)
NFW.set_params(r_decay = 0.1*NFW.r_vir, r_s = 795.0 * cgs.pc )

df_filename = "./LT_NFW_500.dat"

NFW_DF = df.DF(NFW)
f = NFW_DF.load_df(df_filename)



Ncore = 4
Npart = 1.0E4

start = time.time()
NFW_PD = pd.particle_distribution(NFW_DF,Npart, optimize=True)
pos, vel = NFW_PD.parallel_generate_particle_distribution(Npart*100.0,Ncore=Ncore,outfile='lt_100k_optimize_parallel.dat')

#pos, vel = NFW_PD.generate_particle_distribution(Npart*100.0,outfile='lt_100_parallel.dat')

end = time.time()

print end - start


# plt.plot
#plt.plot(NFW_PD._cumulative_mass_r - np.log10(cgs.pc), NFW_PD._cumulative_mass_m)
#plt.plot(NFW_PD._cumulative_mass_r - np.log10(cgs.pc), np.log10(NFW_PD._interpolate_cumulative_mass(10.0**(NFW_PD._cumulative_mass_r))),color='green')
#plt.savefig('cum_mass.png')
#plt.close()
#plt.scatter(NFW_PD._relative_potential_r - np.log10(cgs.pc), NFW_PD._relative_potential_psi)
#plt.plot(NFW_PD._relative_potential_r - np.log10(cgs.pc), np.log10(NFW_PD._interpolate_relative_potential(10.0**(NFW_PD._relative_potential_r))),color='green')



#plt.savefig('relative_potential.png')
