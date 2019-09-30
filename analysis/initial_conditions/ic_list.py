import profiles as prof
import cgs as cgs
import numpy as np
import copy

from scipy.optimize import minimize

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
            print self.name,  "MUST SET A 'potential_type'"
            return

        # check to make sure correct parameters are in dict
        # require: 1) r_s 2) T 
        if not 'b' in self.ic.keys() and not 'r_s' in self.ic.keys():
            print self.name,  "You must set the scale radius as 'b' or 'r_s'"
        elif not 'b' in self.ic.keys():
            self.ic['b'] = self.ic['r_s']
        else:
            self.ic['r_s'] = self.ic['b']

        # if no mu's are set, assume them
        if not 'mu_dwarf' in self.ic.keys():
            print self.name,  "Assuming primordial neutral for dwarf mu = 1.31"
            self.ic['mu_dwarf'] = 1.31
        if not 'mu_halo' in self.ic.keys():
            print self.name,  "Assuming primordial ionized for halo mu = 0.6"
            self.ic['mu_halo'] = 0.6

        
        # M200 and n_o MUST BE SET!!!!! If not, try and solve for them
        if not 'M200' in self.ic.keys() or not 'n_o' in self.ic.keys():
#            print self.name,  "If M200 or n_o are not set, they will be solved for"
#            print self.name,  "using the profile choice"

            if not 'M_DM' in self.ic.keys() or\
               not 'r_DM' in self.ic.keys():
                   print self.name,  "MUST SET 'M_DM', 'r_DM' to solve for DM halo prof"
                   return 
            elif (not 'M_HI' in self.ic.keys() or\
                  not 'r_HI' in self.ic.keys()) and\
                 (not ('n_halo' in self.ic.keys() and 'T_halo' in self.ic.keys())) and\
                 (not ('n_o' in self.ic.keys() and 'n_halo' in self.ic.keys() and 'r_HI' in self.ic.keys())):
                print self.name,  "MUST SET EITHER 'M_HI' and 'r_HI' OR 'n_halo' and 'T_halo' OR 'n_o' and 'n_halo' and 'r_HI'"
                return


            
            if ('T_halo' in self.ic.keys() and 'n_halo' in self.ic.keys()) and\
               ('M_HI' in self.ic.keys() or 'r_HI' in self.ic.keys()):
                print self.name,  "ERROR: Both n_halo and T_halo cannot be specified"
                print self.name,  "       while M_HI and r_HI are also specified"
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
                    print self.name,  "ERROR: Either n_halo or T_halo must be specified"
                    return

            elif 'T_halo' in self.ic.keys() and 'n_halo' in self.ic.keys():
                T_halo = self.ic['T_halo']
                n_halo = self.ic['n_halo']
                M_HI   = -1.0
                r_HI   = -1.0
                print self.name,  't and n are both set for the halo', T_halo, n_halo

            elif 'n_o' in self.ic.keys() and 'n_halo' in self.ic.keys() and\
                 'r_HI' in self.ic.keys():
                T_halo = None
                n_halo = self.ic['n_halo']
                M_HI   = -1.0
                r_HI   = self.ic['r_HI']
                n_o    = self.ic['n_o']
          
            else:
                print self.name,  "T_halo and n_halo must be set"
                return

 
            if self.ic['potential_type'] == 'NFW':

                if not ('rho_s' in self.ic.keys()):
                    self.ic['rho_s'] = None
                

                if 'M_HI' in self.ic.keys() and 'r_HI' in self.ic.keys():                             
                    n_o = None
                elif 'n_o' in self.ic.keys():
                    n_o = self.ic['n_o']
                else:
                    print self.name,  'n_o must be set as an initial condition'
                    return

                if 'T_halo' in self.ic.keys():
                    T_halo = self.ic['T_halo']
                else:
                    T_halo = None

                if 'n_halo' in self.ic.keys():
                    n_halo = self.ic['n_halo']
                else:
                    n_halo = None

                c, r_s, M200, n_o, T_halo, n_halo, rmatch=\
                                   prof.solve_NFW(self.ic['M_DM'], self.ic['r_DM'],
                                   self.ic['r_s'] , M_HI, r_HI,
                                   self.ic['T_dwarf'],
                                   mu=self.ic['mu_dwarf'],
                                   mu_halo=self.ic['mu_halo'],
                                   T_halo = T_halo, n_halo= n_halo,
                                   rho_crit = self.ic['rho_crit'],
                                   n_o = n_o, rho_s = self.ic['rho_s'])

                if (self.ic['b'] is None) or (self.ic['r_s'] is None):
                    self.ic['b'] = r_s
                    self.ic['r_s'] = r_s

                self.ic['M200'] = M200
                self.ic['n_o' ] = n_o
                self.ic['T_halo'] = T_halo
                self.ic['n_halo'] = n_halo
                self.ic['c']      = c
                self.ic['RM']     = rmatch

                self.Pcorona = self.ic['n_halo']*cgs.kb*self.ic['T_halo']

                if not 'M_HI' in self.ic.keys() and 'r_HI' in self.ic.keys():
                    r = np.linspace(0.0, self.ic['r_HI'],1.0E3)
                    rho, RM = prof.NFW_isothermal_rmatch(r, r_s=self.ic['b'],
                           M200=self.ic['M200'], T=self.ic['T_dwarf'],
                           rmatch = self.ic['r_HI'], mu = self.ic['mu_dwarf'],
                           rho_crit = self.ic['rho_crit'],
                           Pcorona = self.Pcorona)# solve the profile for total mass

                    M_HI = prof.cumulative_mass(r,rho)
                    self.ic['M_HI'] = M_HI[-1]

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

        if not 'R200' in self.ic.keys():
            R200 = (self.ic['M200'] * (3.0 / (4.0*np.pi))/(200.0 * self.ic['rho_crit']))**(1.0/3.0)
            self.ic['R200'] = R200

        if not 'V_max' in self.ic.keys():
#            res = minimize( lambda x : -self.circular_velocity(x), self.ic['b'])
#            self.ic['V_max' ] = -res.fun[0] # result in cm/s
#            self.ic['R_vmax'] = res.x[0]
            r = np.logspace(np.log10(0.01 * self.ic['b']), np.log10(10*self.ic['b']), 100000)
            vc = self.circular_velocity(r)
            self.ic['V_max'] = np.max(vc)
            self.ic['R_vmax'] = r[np.argmax(vc)]

        if not 'V_vir' in self.ic.keys():
            self.ic['V_vir'] = self.circular_velocity(self.ic['R200'])
        if not 'V_s' in self.ic.keys():
            self.ic['V_s'] = self.circular_velocity(self.ic['b'])

        if not 'T_vir' in self.ic.keys():
            self.ic['T_vir'] = (1.0 / 3.0) * (3.0/2.0) * cgs.mu_ionized * cgs.mp / cgs.kb * self.ic['V_vir']**2



        if self.ic['potential_type'] == 'Burkert':
            if (not 'T_halo' in self.ic.keys()):
                rho = self.find_density_profile(self.ic['r_HI'], type='Burkert_isothermal')

                self.ic['T_halo'] = rho / (self.ic['mu_dwarf'] * cgs.mp) * self.ic['T_dwarf'] / self.ic['n_halo']

            if (not 'RM' in self.ic.keys()):
                self.ic['RM'] = self.ic['r_HI']

        self.Pcorona = self.ic['n_halo']*cgs.kb*self.ic['T_halo']

        self.SF_ic()

    def DM_density(self, r, type='NFW'):
        if 'potential_type' in self.ic.keys():
            type = self.ic['potential_type']

        if type == 'NFW':
            rho_DM = prof.NFW_DM(r, r_s=self.ic['b'],
                                M200=self.ic['M200'],
                                rho_crit = self.ic['rho_crit'])

        elif type == 'Burkert' or  type =='Burkert_isothermal':
            rho_DM = prof.burkert_DM(r, r_s=self.ic['b'], M200 = self.ic['M200'], rho_crit=self.ic['rho_crit'])
        self.rho_DM = rho_DM
        return rho_DM

    def M_r(self, r, type='NFW'):

        if 'potential_type' in self.ic.keys():
            type = self.ic['potential_type']

        if type == 'Burkert' or type == 'Burkert_isothermal':
            M = prof.burkert_mass(r, self.ic['b'], self.ic['M200'], self.ic['rho_crit'])

        elif type == 'NFW':
            M = prof.NFW_mass(r, self.ic['b'], self.ic['M200'], self.ic['rho_crit'])

        return M

    def circular_velocity(self, r):
        return np.sqrt(cgs.G * self.M_r(r) / r)

    def find_density_profile(self, r, type='NFW_isothermal'):

        if 'potential_type' in self.ic.keys():
            type = self.ic['potential_type']

        if type == 'NFW_isothermal' or type == 'NFW':
            function = prof.NFW_isothermal_gas

            rho = function(r, r_s=self.ic['b'],
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


        elif type == 'Burkert_isothermal' or type == 'Burkert':
            rho  = prof.Burkert_isothermal_gas(r, r_s=self.ic['b'],
                                 M200=self.ic['M200'], T=self.ic['T_dwarf'],
                               n_o = self.ic['n_o'], mu = self.ic['mu_dwarf'],
                               rho_crit=self.ic['rho_crit'])  

        self.rvals  = r
        self.rho    = rho
       # self.radius = RM
        return rho

    def column_density(self, r, type = 'NFW_isothermal', **kwargs):

        if 'potential' in self.ic.keys():
            type = self.ic['potential']

        if type == 'Burkert_isothermal' or type == 'Burkert':
            density_function = lambda x : prof.Burkert_isothermal_gas(x, r_s=self.ic['b'],
                           M200=self.ic['M200'], T=self.ic['T_dwarf'],
                           n_o = self.ic['n_o'], mu = self.ic['mu_dwarf'],
                           rho_crit = self.ic['rho_crit'])/ (self.ic['mu_dwarf']*cgs.mp)

            NHI = prof.column_density(r, density_function=density_function, **kwargs)


        elif type == 'NFW_isothermal' or type == 'NFW':
            density_function = lambda x : prof.NFW_isothermal_gas(x, r_s=self.ic['b'],
                           M200=self.ic['M200'], T=self.ic['T_dwarf'],
                           n_o = self.ic['n_o'], mu = self.ic['mu_dwarf'],
                           rho_crit = self.ic['rho_crit'],
                           Pcorona = self.Pcorona) / (self.ic['mu_dwarf']*cgs.mp)

            NHI = prof.column_density(r, density_function=density_function, **kwargs)



        self.NHI = NHI


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
    def SF_ic(self):

        if not 'r_HI' in self.ic.keys():
            r_HI = self.ic['RM']/cgs.pc
        else: 
            r_HI = self.ic['r_HI'] / cgs.pc

        if not 'M_HI' in self.ic.keys():
            return 
        else:
            M_HI = self.ic['M_HI'] / cgs.Msun

        SFR = ((M_HI/(np.pi*r_HI*r_HI))** (2.47) * 2.13E-5)

        SFR = (SFR * (np.pi*r_HI*r_HI/1000.0/1000.0))
        SNR = 6.0E-3 * SFR

        self.ic['SFR'] = SFR # solar masses per year
        self.ic['SNR'] = SNR # supernova per year

        # SNR in Myr per kpc^2        
        SNR_surface = SNR * 1.0E6 / (np.pi * (r_HI/1000.0) * (r_HI/1000.0)) 

        self.ic['Gamma_pe'] = 8.5E-26 * (SNR_surface/11.0)**(2.0/7.0)

   
 
    
    
      


known_initial_conditions = {'CarinaMidMed': # see Table 4 in Gatto et. al.
                            {'n_o'  : 0.4*1.36, 'mu_dwarf' : 1.297,
                             'mu_halo' : 0.62,
                             'T_halo': 1.8E6, 'n_halo' : 1.7E-4,
                             'M_DM'  : 3.7E7 * cgs.Msun, 
                             'r_DM'  : 0.87  * cgs.kpc,
                             'T_dwarf': 1.0E4,
                             'b'     :  795.0*cgs.pc,### Walker et. al. 2009 +erratum ###,
                             'potential_type' : 'NFW'},
                            'SextansNFW' : {
                               'n_o' : 1.0, 'n_halo' :1.0E-4, 'T_halo' : 2.0E6,
                               'T_dwarf' : 1.0E4, 'M_DM' : 2.5E7 * cgs.Msun, 'r_DM' : 682.0*cgs.pc,
                               #'rho_s' : 1.9*cgs.Msun/cgs.pc**3, 
                               'potential_type' : "NFW", 'r_s' : 800*cgs.pc},
                            'DracoNFW' : {
                               'n_o' : 1.0, 'n_halo' :1.0E-4, 'T_halo' : 2.0E6,
                               'T_dwarf' : 1.0E4, 'M_DM' : 9.4E6 * cgs.Msun, 'r_DM' : 196.0*cgs.pc,
                               #'rho_s' : 3.0*cgs.Msun/cgs.pc**3, 
                               'potential_type' : "NFW", 'r_s' : 800*cgs.pc},
                            'FornaxNFW' : {
                               'n_o' : 1.0, 'n_halo' :1.0E-4, 'T_halo' : 2.0E6,
                               'T_dwarf' : 1.0E4, 'M_DM' : 5.3E7 * cgs.Msun, 'r_DM' : 668.0*cgs.pc,
                               #'rho_s' : 4.2*cgs.Msun/cgs.pc**3, 
                               'potential_type' : "NFW", 'r_s' : 800.0*cgs.pc},
                            'CarinaNFW' : {
                               'n_o' : 1.0, 'n_halo' :1.0E-4, 'T_halo' : 2.0E6,
                               'T_dwarf' : 1.0E4, 'M_DM' : 6.1E6 * cgs.Msun, 'r_DM' : 241.0*cgs.pc,
                               #'rho_s' : 1.0*cgs.Msun/cgs.pc**3, 
                                'potential_type' : "NFW", 'r_s' : 800*cgs.pc},
                            'LeoTNFW' : {
                               'n_o' : 1.0, 'n_halo' :1.0E-4, 'T_halo' : 2.0E6,
                               'T_dwarf' : 1.0E4, 'M_DM' : 5.8E6 * cgs.Msun, 'r_DM' : 178.0*cgs.pc,
                               #'rho_s' : 2.5 *cgs.Msun/cgs.pc**3, 
                               'potential_type' : "NFW", 'r_s' : 800*cgs.pc},
                            'LMCNFW' : {
                               'n_o' : 1.0, 'n_halo' :1.0E-4, 'T_halo' : 2.0E6,
                               'T_dwarf' : 1.0E4, 'M_DM' : 1.4E10 * cgs.Msun, 'r_DM' : 8.7*cgs.kpc,
                               #'rho_s' : 2.5 *cgs.Msun/cgs.pc**3, 
                               'potential_type' : "NFW", 'r_s' : None, 'rho_s' : 3.4E-24},
                            'SextansMidMed' :  {
                             'n_o' : 1.36 * 0.27, 'mu_dwarf' : 1.297,
                             'mu_halo' : 0.62, 
                             'n_halo' : 1.8E-4, 'T_halo' : 1.8E6,
                             'T_dwarf' : 1.0E4,
                             'b'  : 795.0 * cgs.pc,
                             'M_DM' : 2.0E7 * cgs.Msun, 'r_DM' : 1.0*cgs.kpc, 
                             'potential_type' : 'NFW'},
                            'SMC_bubble': {
                             'T_dwarf' : 10000.0, 'M_DM' : 1.95E10 * cgs.Msun,
                             'r_DM'    : 57.1 * cgs.kpc,
                             'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                             'b' : 2.5 * cgs.kpc, 'potential_type' : 'Burkert',
                             'M_HI' : 7.9E8*cgs.Msun, 'r_HI' : 2000.0*cgs.pc,
                             'n_halo' : 1.0E-7 },
                            'CarinaMidHigh': 
                             {'n_o' : 0.4*1.36, 'mu_dwarf' : 1.297,
                             'mu_halo' : 0.62,
                             'n_halo' : 6.8E-4, 'T_halo' : 0.45E6,
                             'M_DM' : 3.7E7 * cgs.Msun,
                             'r_DM' : 0.87 * cgs.kpc,
                             'T_dwarf' : 1.0E4,
                             'b' : 795.0*cgs.pc,
                             'potential_type' : 'NFW'},        
                            'Leo_T_obs':
                            {'T_dwarf' : 6000.0, 'M_DM' :7.3E6*cgs.Msun,
                             'r_DM': 300.*cgs.pc, 'r_HI':300.0*cgs.pc,
                             'M_HI' : 2.8E5*cgs.Msun,
                             'b'    : 795.0*cgs.pc,
                             'n_halo' : 4.5E-4, 
                             #'T_halo' : 7.5E5,
                             'potential_type':'NFW'},

                             'LT_n020_v2_nh5':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.20, 'n_halo' : 1.0E-5, 'v_halo' : 200.0E5},
                            'LT_n075_v2_nh5':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 1.0E-5, 'v_halo' : 200.0E5},
                             'LT_n150_v2_nh5':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-5, 'v_halo' : 200.0E5},
                             'LT_n020_v2_nh3':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.20, 'n_halo' : 1.0E-3, 'v_halo' : 200.0E5},
                            'LT_n075_v2_nh3':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 1.0E-3, 'v_halo' : 200.0E5},
                             'LT_n150_v2_nh35':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 10.0**(-3.5), 'v_halo' : 200.0E5},
                             'LT_n150_v2_nh4_exp':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 400.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-4, 'v_halo' : 200.0E5},
                             'LT_n150_v2_nh4_exp_2':
                            {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 400.0 * cgs.pc,
                              'M_HI'    : 4.74E5 * 1.99E33, 'n_halo' : 1.0E-4, 'v_halo' : 200.0E5},
                             'LT_n150_v2_nh3':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-3, 'v_halo' : 200.0E5},
                             'LT_n020_v2_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.20, 'n_halo' : 1.0E-4, 'v_halo' : 200.0E5},
                            'LT_n075_v2_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 1.0E-4, 'v_halo' : 200.0E5},
                            'LT_n075_v2_nh4_4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 4.0E-4, 'v_halo' : 200.0E5},
                            'LT_n150_v2_nh4_4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 4.0E-4, 'v_halo' : 200.0E5},
                             'LT_n150_v2_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-4, 'v_halo' : 200.0E5},
                            'LT_n192_v2_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.92, 'n_halo' : 1.0E-4, 'v_halo' : 200.0E5},
                            'LT_n192_v4_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.92, 'n_halo' : 1.0E-4, 'v_halo' : 400.0E5},
                            
                            ##### same as the above but with v = 400.0E5
                             'LT_n020_v4_nh5':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.20, 'n_halo' : 1.0E-5, 'v_halo' : 400.0E5},
                            'LT_n075_v4_nh5':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 1.0E-5, 'v_halo' : 400.0E5},
                             'LT_n150_v4_nh5':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-5, 'v_halo' : 400.0E5},
                             'LT_n020_v4_nh3':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.20, 'n_halo' : 1.0E-3, 'v_halo' : 400.0E5},
                            'Hu':
                             {'T_dwarf' : 6000.0, "M_DM" : 1.0E10 * cgs.Msun,
                              'r_DM' : 44.0 * cgs.kpc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 4400.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 1000.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'Hu-NFW':
                             {'T_dwarf' : 6000.0, "M_DM" : 1.0E10 * cgs.Msun,
                              'r_DM' : 44.0 * cgs.kpc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 4400.0 * cgs.pc, 'potential_type' : "NFW",
                              'r_HI' : 2000.0 * cgs.pc, 'M_HI' : 1.4E8, 'n_halo':1.0E-6},
                            'Leo_P':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 369.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'WLM':
                             {'T_dwarf' : 1.0E4, "M_DM" : 1.0E10 * cgs.Msun,
                               'r_DM' : 45 * cgs.kpc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                               'b'    : 3 * cgs.kpc, 'potential_type' : 'NFW', 'r_HI' : 6.0*cgs.kpc,
                               'M_HI' : 7.0E7 * cgs.Msun, 'n_halo' : 1.0E-5},
                            'Leo_P_2rs':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 738.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'Leo_P_30km':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 985.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'Leo_P_40km':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 1415.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'Leo_P_50km':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, "mu_halo" : 0.6, "mu_dwarf":1.31,
                              'b' : 1850.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo' : 1.0E-4},
                            'Leo_P_4':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 900.0 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'Leo_P_3':
                             {'T_dwarf' : 6000.0, "M_DM" : 2.67E7 * cgs.Msun,
                              'r_DM' : 500.0 * cgs.pc, 'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 645.75 * cgs.pc, 'potential_type' : "Burkert",
                              'r_HI' : 500.0 * cgs.pc, 'M_HI' : 1.4E6, 'n_halo':1.0E-4},
                            'Forbes':
                             {'T_dwarf' : 6000.0, "M200" : 1.0E10*cgs.Msun, "n_o" : 100.0,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b' : 2000.0 * cgs.pc, 'potential_type' : "Burkert",
                              'M_HI' : 1.0E7*cgs.Msun, 'r_HI' : 8.0*cgs.kpc, 'n_halo':1.0E-4},
                            'LT_n075_v4_nh3':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 1.0E-3, 'v_halo' : 400.0E5},
                             'LT_n150_v4_nh3':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-3, 'v_halo' : 400.0E5},
                             'LT_n020_v4_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.20, 'n_halo' : 1.0E-4, 'v_halo' : 400.0E5},
                            'LT_n075_v4_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 0.75, 'n_halo' : 1.0E-4, 'v_halo' : 400.0E5},
                             'LT_n150_v4_nh4':
                             {'T_dwarf' : 6000.0, "M_DM" : 7.3E6 * cgs.Msun,
                              'r_DM'    : 300.0 *cgs.pc,
                              'mu_halo' : 0.6, 'mu_dwarf' : 1.31,
                              'b'       : 795.0 * cgs.pc,
                              'potential_type' : 'NFW',
                              'r_HI'    : 300.0 * cgs.pc,
                              'n_o' : 1.50, 'n_halo' : 1.0E-4, 'v_halo' : 400.0E5},                            
                                                            
                              ##############################
                             
                      #      'Leo_T_obs_burkert':
                      #      {'T_dwarf' : 1.0E4, 'M_DM' : 1.0E7*cgs.Msun,
                      #       'r_DM': 300.*cgs.pc, 'r_HI':300.0*cgs.pc,
                      #       'M_HI' : 2.8E5*cgs.Msun,
                      #       'b'    : 500.0*cgs.pc,
                      #       'n_halo' : 6.E-5,
                      #       'potential_type':'Burkert'},
                            # from Faerman et. al. 2013
                            'Leo_T_Burkert':
                            {'T_dwarf' : 6.0E3,
                             'M_DM'   : 7.3E6 * cgs.Msun, 'r_DM' : 300 * cgs.pc,
                             'b'      : 105.4229* cgs.pc, 
                             #'T_halo' : 8.0E5,
                             'n_halo' : 1.0E-4,
                             'M_HI' : 4.74E5 * cgs.Msun,
                             'r_HI'   : 300.0*cgs.pc,
                             #'n_o'    : 2.816185897E-3,
                             'mu_dwarf' : 1.31,
                             'mu_halo'  : 0.6,
                             'potential_type':'Burkert'},
#
#                           'Leo_T_solve_burkert':
#                            {'T_dwarf' : 6000.0, 'M_DM' : 7.3E6*cgs.Msun,
#                             'r_DM': 300.*cgs.pc, 'r_HI':300.0*cgs.pc,
#                             'M_HI' : 2.8E5*cgs.Msun,
#                       #      'b'    : 708.627233*cgs.pc,
#                             'r_s' : 472.4*cgs.pc,
#                             'n_halo': 4.5E-4,
#                             #'T_halo' : 7.5E5, 
#                             'potential_type':'Burkert'},  
                             

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
    ic_object_dict[known_dwarf] = dwarf_ic(known_dwarf)
    ic_object_dict[known_dwarf].set_ic(copy.deepcopy(known_initial_conditions[known_dwarf]))
    
print "Loaded IC's for ", num_dwarfs, " dwarf galaxies"
