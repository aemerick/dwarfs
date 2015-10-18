import numpy as np
import cgs   as cgs
import dm_density_profiles as dmd
import df as df
import time

# output filename and number of points to generate
npoints     = 750
df_filename = "LT_NFW_%0004i.dat"%(npoints)

# set up the halo
# using the dm  profile
NFW = dmd.general_dm_profile('LT_NFW')             # give it a name
NFW.set_params(profile_shape_params=[1.0,3.0,1.0]) # NFW profile
NFW.set_params(M_vir = 3.0934E8 * cgs.Msun, r_vir = 13.70* cgs.kpc) # params
NFW.set_params(r_decay = 0.1*NFW.r_vir, r_s = 795.0 * cgs.pc ) # more params

# make DF object
NFW_DF = df.DF(NFW)

# generate DF
start = time.time()

f = NFW_DF.compute(npoints, filename = df_filename)

end   = time.time()

print "Generated DF in ", (end - start)/60.0, " minutes."
