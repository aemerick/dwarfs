import profiles as prof
import cgs as cgs
import numpy as np

class dwarf_ic:

    def __init__(self, name):
        self.name = name

    def set_ic(self, ic_dict):
        self.ic = ic_dict

        if not 'rho_crit' in self.ic.keys():
            self.ic['rho_crit'] = 9.74E-30 # z=0 crit. dens. of U

        self.Pcorona = self.ic['n_halo']*cgs.kb*self.ic['T_halo']
                                        
    def find_density_profile(self, r, type='NFW_isothermal'):
               
        if type == 'NFW_isothermal':
            function = prof.NFW_isothermal_gas
        

            rho, RM = function(r, r_s=self.ic['b'],
                           M200=self.ic['M200'], T=self.ic['T_dwarf'], 
                           n_o = self.ic['n_o'], mu = self.ic['mu_dwarf'], 
                           rho_crit = self.ic['rho_crit'],
                           Pcorona = self.Pcorona)

        elif type == 'NFW_isothermal_rmatch':
            rho, RM = prof.NFW_isothermal_rmatch(r, r_s=self.ic['b'],
                           M200=self.ic['M200'], T=self.ic['T_dwarf'],
                           rmatch = self.ic['rmatch'], mu = self.ic['mu_dwarf'],
                           rho_crit = self.ic['rho_crit'],
                           Pcorona = self.Pcorona)

        self.rvals  = r
        self.rho    = rho
        self.radius = RM

    



known_initial_conditions = {'Sextans_test':
                               {'T_dwarf' : 1.0E4, 'T_halo': 1.8E6,
                                'mu_halo' : 0.6  , 'mu_dwarf': 1.31,
                                'M200'    : 2.0E7*cgs.Msun,
                                'b'       : 170.0*cgs.pc,
                                'n_o'     : 0.27, 'n_halo': 1.8E-4},

                            'Leo_T' :
                               {'T_dwarf' : 1.0E4, 'T_halo': 8.45955E5,
                                'mu_halo' : 0.6  , 'mu_dwarf':1.31,
                                'M200': 4.5166E8*cgs.Msun,
                                'b'   : 800.0*cgs.pc,
                                'n_o' : 0.084692897, 'n_halo': 6.0E-5},
                            # from jana's ic.F90
                            'Leo_test'    :
                               {'T_dwarf' : 9.182E3, 'T_halo': 7.40914E5,
                                'mu_halo' : 0.6  , 'mu_dwarf': 1.31,
                                'M200'    : 1.0E7*cgs.Msun,
                                'b'       : 500.0*cgs.pc,
                                'rmatch'  : 300.0*cgs.pc,
#                                'n_o'     : 1.0,
                                'n_halo'  : 4.6E-5} 
                           }
                                        
num_dwarfs = len(known_initial_conditions.keys())
dwarf_dict = {}

for known_dwarf in known_initial_conditions:
    dwarf_dict[known_dwarf] = dwarf_ic(known_dwarf)
    dwarf_dict[known_dwarf].set_ic(known_initial_conditions[known_dwarf])
    
 
print "Loaded IC's for ", num_dwarfs, " dwarf galaxies"
