from initial_conditions import ic_list as icl
import numpy as np
import cgs as cgs

LT_names = []
LT_list  = []
for ic in icl.ic_object_dict.keys():
  if ic.startswith('LT_n'):
      LT_names.append(ic)
      LT_list.append(icl.ic_object_dict[ic])

print "n_o n_halo v_halo T_halo M_HI M_200 r_HI R200"

for LT_icl in LT_list:

    ic = LT_icl.ic

    print "%.3f %.2E %.3f %.2E %5.5E %5.5E %5.5E %5.5E"%(ic['n_o'], ic['n_halo'], ic['v_halo']/cgs.km, ic['T_halo'], ic['M_HI']/cgs.Msun, ic['M200']/cgs.Msun, ic['r_HI']/cgs.pc, ic['R200']/cgs.pc)
      


