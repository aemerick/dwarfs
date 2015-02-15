import profiles as prof
import cgs as cgs
import numpy as np

class dwarf_ic:

    def __init__(self, name):
        self.name = name

    def set_ic(self, ic_dict):
        """
            Provided a dict containing the initial conditions, 
            makes sure all the IC's are there and valid to generate
            a legitimate initial condition. If not, attempts to solve
            for them assuming the correct additional IC's are provided.

            More info to come.....
        """

        self.ic = ic_dict

        if not 'rho_crit' in self.ic.keys():
            self.ic['rho_crit'] = 9.74E-30 # z=0 crit. dens. of U


        if not 'potential_type' in self.ic.keys():
            print "MUST SET A 'potential_type'"
            return

        # check to make sure correct parameters are in dict
        # require: 1) r_s 2) T 
        if not 'b' in self.ic.keys() and not 'r_s' in self.ic.keys():
            print "You must set the scale radius as 'b' or 'r_s'"
        elif not 'b' in self.ic.keys():
            self.ic['b'] = self.ic['r_s']
        else:
            self.ic['r_s'] = self.ic['b']

        # if no mu's are set, assume them
        if not 'mu_dwarf' in self.ic.keys():
            print "Assuming primordial neutral for dwarf mu = 1.31"
            self.ic['mu_dwarf'] = 1.31
        if not 'mu_halo' in self.ic.keys():
            print "Assuming primordial ionized for halo mu = 0.6"
            self.ic['mu_halo'] = 0.6

        
        # M200 and n_o MUST BE SET!!!!! If not, try and solve for them
        if not 'M200' in self.ic.keys() or not 'n_o' in self.ic.keys():
            print "If M200 or n_o are not set, they will be solved for"
            print "using the profile choice"

            if not 'M_DM' in self.ic.keys() or\  
               not 'M_HI' in self.ic.keys() or\
               not 'r_DM' in self.ic.keys() or\
               not 'r_HI' in self.ic.keys():
                   print "MUST SET 'M_DM' 'M_HI' 'r_DM' 'r_HI' to solve"
                   return

            if 'T_halo' in self.ic.keys() and 'n_halo' in self.ic.keys():
                print "ERROR: Both n_halo and T_halo canno be specified"
                return
            elif 'T_halo' in self.ic.keys():
                T_halo = self.ic['T_halo']
                n_halo = None
            elif 'n_halo' in self.ic.keys():
                n_halo = self.ic['n_halo']
                T_halo = None
            else:
                print "ERROR: Either n_halo and T_halo must be specified"
                return

 
            if self.ic['potential_type'] == 'NFW':
    return r_s, M200, n_o, T_halo, n_halo
                

                c, r_s, M200, n_o, T_halo, n_halo =\
                               prof.solve_NFW(self.ic['M_DM'], self.ic['r_DM'],
                               self.ic['r_s'] , self.ic['M_HI'],
                               self.ic['r_HI'], self.ic['T'],
                               mu=self.ic['mu_dwarf'],
                               mu_halo=self.ic['mu_halo'],
                               T_halo = T_halo, n_halo= n_halo,
                               rho_crit = self.ic['rho_crit'])
                self.ic['M200'] = M200
                self.ic['n_o' ] = n_o
                self.ic['T_halo'] = T_halo
                self.ic['n_halo'] = n_halo
                self.ic['c']      = c

            elif self.ic['potential_type'] == 'Burkert':

                r_s, M200, n_o, T_halo, n_halo=\
                      prof.solve_Burkert(self.ic['M_DM'], self.ic['r_DM'],
                               self.ic['r_s'] , self.ic['M_HI'],
                               self.ic['r_HI'], self.ic['T'],
                               mu_dwarf=self.ic['mu_dwarf'],
                               mu_halo=self.ic['mu_halo']
                               T_halo = T_halo, n_halo= n_halo,
                               rho_crit = self.ic['rho_crit'])

                self.ic['M200'] = M200
                self.ic['n_o' ] = n_o
                self.ic['T_halo'] = T_halo
                self.ic['n_halo'] = n_halo



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


        elif type == 'Burkert_isothermal':
            rho, RM = prof.Burkert_isothermal_gas(r, r_s=self.ic['b'],
                                 M200=self.ic['M200'], T=self.ic['T_dwarf'],
                               n_o = self.ic['n_o'], mu = self.ic['mu_dwarf'],
                               rho_crit=self.ic['rho_crit'])  

        self.rvals  = r
        self.rho    = rho
        self.radius = RM


    def FLASH_readable_ic(filename=None):
        """
           Saves to file or prints to screen the IC's
           in a fashion that can be copy pasted into
           the flash.par file.
        """
        
        
        print "currently does nothing"
        return


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
