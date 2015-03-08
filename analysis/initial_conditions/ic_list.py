import profiles as prof
import cgs as cgs
import numpy as np
import copy

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
               not 'r_DM' in self.ic.keys():
                   print "MUST SET 'M_DM', 'r_DM' to solve for DM halo prof"
                   return 
            elif (not 'M_HI' in self.ic.keys() or\
                  not 'r_HI' in self.ic.keys()) and\
             (not ('n_halo' in self.ic.keys() and 'T_halo' in self.ic.keys())):
                  print "MUST SET EITHER 'M_HI' and 'r_HI' OR 'n_halo' and 'T_halo'"
                  return


            
            if ('T_halo' in self.ic.keys() and 'n_halo' in self.ic.keys()) and\
               ('M_HI' in self.ic.keys() or 'r_HI' in self.ic.keys()):
                print "ERROR: Both n_halo and T_halo cannot be specified"
                print "       while M_HI and r_HI are also specified"
                return

            elif 'M_HI' in self.ic.keys() and 'r_HI' in self.ic.keys():
                M_HI = self.ic['M_HI'] ; r_HI = self.ic['r_HI']
                if 'T_halo' in self.ic.keys():
                    T_halo = self.ic['T_halo']
                    n_halo = None
                elif 'n_halo' in self.ic.keys():
                    n_halo = self.ic['n_halo']
                    T_halo = None

                else:
                    print "ERROR: Either n_halo and T_halo must be specified"
                    return

            elif 'T_halo' in self.ic.keys() and 'n_halo' in self.ic.keys():
                T_halo = self.ic['T_halo']
                n_halo = self.ic['n_halo']
                M_HI   = -1.0
                r_HI   = -1.0

          
            else:
                print "T_halo and n_halo must be set"
                return

 
            if self.ic['potential_type'] == 'NFW':


                if 'M_HI' in self.ic.keys() and 'r_HI' in self.ic.keys():                             
                    n_o = None
                elif 'n_o' in self.ic.keys():
                    n_o = self.ic['n_o']
                else:
                    print 'n_o must be set as an initial condition'
                    return
          
 
                c, r_s, M200, n_o, T_halo, n_halo, rmatch=\
                                   prof.solve_NFW(self.ic['M_DM'], self.ic['r_DM'],
                                   self.ic['r_s'] , M_HI, r_HI,
                                   self.ic['T_dwarf'],
                                   mu=self.ic['mu_dwarf'],
                                   mu_halo=self.ic['mu_halo'],
                                   T_halo = T_halo, n_halo= n_halo,
                                   rho_crit = self.ic['rho_crit'],
                                   n_o = n_o)
             

                self.ic['M200'] = M200
                self.ic['n_o' ] = n_o
                self.ic['T_halo'] = T_halo
                self.ic['n_halo'] = n_halo
                self.ic['c']      = c
                self.ic['RM']     = rmatch

            elif self.ic['potential_type'] == 'Burkert':

                r_s, M200, n_o, T_halo, n_halo=\
                      prof.solve_burkert(self.ic['M_DM'], self.ic['r_DM'],
                               self.ic['r_s'] , self.ic['M_HI'],
                               self.ic['r_HI'], self.ic['T_dwarf'],
                               mu_dwarf=self.ic['mu_dwarf'],
                               mu_halo=self.ic['mu_halo'],
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
            rho  = prof.Burkert_isothermal_gas(r, r_s=self.ic['b'],
                                 M200=self.ic['M200'], T=self.ic['T_dwarf'],
                               n_o = self.ic['n_o'], mu = self.ic['mu_dwarf'],
                               rho_crit=self.ic['rho_crit'])  

        self.rvals  = r
        self.rho    = rho
       # self.radius = RM


    def FLASH_readable_ic(self,filename=None):
        """
           Saves to file or prints to screen the IC's
           in a fashion that can be copy pasted into
           the flash.par file.
        """
        
        params = {'T_dwarf' : ['sim_TCloud',self.ic['T_dwarf']],
                  'T_halo'  : ['sim_TAmbient',self.ic['T_halo']],
                  'M200'    : ['sim_M200',self.ic['M200']],
                  'n_halo'  : ['sim_rhoAmbient',self.ic['mu_halo']*cgs.mp*self.ic['n_halo']],
                  'n_o'     : ['sim_rhoCenter',self.ic['mu_dwarf']*cgs.mp*self.ic['n_o']],
                  'mu_dwarf': ['sim_mu_dwarf',self.ic['mu_dwarf']],
                  'mu_halo' : ['sim_mu_halo',self.ic['mu_halo']],
                  'b'       : ['sim_bparam',self.ic['b']],
                  'rho_crit': ['sim_rho_crit',self.ic['rho_crit']]}
        
 
        output_form = "{:<18} = {:8.8E}"
        if (filename == None):        
            for p in params:
                print output_form.format(params[p][0],params[p][1])
        else:
           format = format + "\n"
           f = open(filename, 'w')
           for p in params:
               f.write(output_form.format(params[p][0],params[p][1]))
 
           f.close()
      
      
   # def SF_ic(sfr=1.0):
    
    
    
      


known_initial_conditions = {'CarinaMidMed': # see Table 4 in Gatto et. al.
                            {'n_o'  : 0.4, 'mu_dwarf' : 1.31,
                             'T_halo': 1.8E6, 'n_halo' : 1.7E-4,
                             'M_DM'  : 3.7E7 * cgs.Msun, 
                             'r_DM'  : 0.87  * cgs.kpc,
                             'T_dwarf': 1.0E4,
                             'b'     :  795.0*cgs.pc,### Walker et. al. 2009 +erratum ###,
                             'potential_type' : 'NFW'},

                            'Leo_T_obs':
                            {'T_dwarf' : 6000.0, 'M_DM' :7.3E6*cgs.Msun,
                             'r_DM': 300.*cgs.pc, 'r_HI':300.0*cgs.pc,
                             'M_HI' : 2.8E5*cgs.Msun,
                             'b'    : 795.0*cgs.pc,
                             'n_halo' : 4.5E-4, 
                             #'T_halo' : 7.5E5,
                             'potential_type':'NFW'},
                      #      'Leo_T_obs_burkert':
                      #      {'T_dwarf' : 1.0E4, 'M_DM' : 1.0E7*cgs.Msun,
                      #       'r_DM': 300.*cgs.pc, 'r_HI':300.0*cgs.pc,
                      #       'M_HI' : 2.8E5*cgs.Msun,
                      #       'b'    : 500.0*cgs.pc,
                      #       'n_halo' : 6.E-5,
                      #       'potential_type':'Burkert'},
                            # from Faerman et. al. 2013
                            'Leo_T_burkert':
                            {'T_dwarf' : 6.0E3,
                             'M200'   : 0.979796E7 * cgs.Msun,
                             'b'    : 472.4 * cgs.pc, 
                             'T_halo' : 8.0E5,
                             'n_halo' : 4.6E-5,
                             'n_o'    : 1.6954E-3,
                             #'n_o'    : 2.816185897E-3,
                             'mu_dwarf' : 1.31,
                             'mu_halo'  : 0.6,
                             'potential_type':'Burkert'},
                           'Leo_T_solve_burkert':
                            {'T_dwarf' : 6000.0, 'M_DM' : 7.3E6*cgs.Msun,
                             'r_DM': 300.*cgs.pc, 'r_HI':300.0*cgs.pc,
                             'M_HI' : 2.8E5*cgs.Msun,
                       #      'b'    : 708.627233*cgs.pc,
                             'r_s' : 472.4*cgs.pc,
                             'n_halo': 4.5E-4,
                             #'T_halo' : 7.5E5, 
                             'potential_type':'Burkert'},  
                             

                             'Sextans_test':
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
ic_object_dict = {}

for known_dwarf in known_initial_conditions:
    print "Loading IC for ", known_dwarf
    ic_object_dict[known_dwarf] = dwarf_ic(known_dwarf)
    ic_object_dict[known_dwarf].set_ic(copy.deepcopy(known_initial_conditions[known_dwarf]))
    
 
print "Loaded IC's for ", num_dwarfs, " dwarf galaxies"
