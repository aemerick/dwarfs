import numpy as np
import cgs as cgs

from dark_matter_ic import df as df
from dark_matter_ic import dm_density_profiles as dmd
from dark_matter_ic import particle_distribution as pd

# settings for generating particles in parallel
df_filename = 'KMM_model_f-newbounds.dat'
outfile     = 'KMM-F-11_9-parallel.dat'
Npart       = 1.0E3
Ncore       = 4
# --------------------------------------------

# Setting up the DM profile (nicer if this can load from a file)
NFW = dmd.general_dm_profile("Model_F")
NFW.set_params(profile_shape_params=[1.0,3.0,1.0])

# assuming h = 0.7
h = 0.7
c = 15.0

NFW.set_params(M_vir = 1.0E10 /h * cgs.Msun)
NFW.set_params(r_s   = 2.94 / h * cgs.kpc)
NFW.set_params(r_vir = NFW.r_s * c)
NFW.set_params(r_decay = 0.1 * NFW.r_vir)
# -------------------------------------------------

# Code running below:

# Load DF
NFW_DF = df.DF(NFW)
f      = NFW_DF.load_df(df_filename)

# Initialize PDF distribution
NFW_PD = pd.particle_distribution(NFW_DF, Npart, optimize=True)
max_loop = np.inf

# generate particles 
pos, vel = NFW_PD.parallel_generate_particle_distribution(max_loop, Ncore=Ncore, outfile=outfile)
#pos, vel = NFW_PD.generate_particle_distribution(max_loop, outfile = outfile)
