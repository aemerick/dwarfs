import os
import yaml
import numpy as np
import yt
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


class dwarf:

    def __init__(self, ds, param_file, raw = True):
        """
        If raw = true, do some grooming first to import the file
        """
    
    
        self.ds = ds
        self.param_file = param_file
        
        if raw:
            newfile = param_file + ".mod"
           # os.system("cp " + param_file + " " + newfile)
           
            # convert = to : and save to new file
            bash_command = "sed 's/=/:/g' " + param_file + " > " + newfile
            os.system(bash_command)
            
            # remove all comments
            bash_command = "sed -i 's:#.*$::g' " + newfile
            os.system(bash_command)
            
            # remove all blank lines
            bash_command = "sed -i '/^$/d' " + newfile
            os.system(bash_command)
            
        else:
            newfile = param_file           
            
            
        stream = open(newfile, 'r')
        self.params     = yaml.load(stream, Loader=Loader)
        
        
        self.center = np.array( [ np.float(self.params['sim_xctr']),
                                  np.float(self.params['sim_yctr']),
                                  np.float(self.params['sim_zctr']) ])

        self.center = self.center * yt.units.cm   
                                       
        self.radius = np.float( self.params['sim_RL'] ) * yt.units.cm
        
        
       # self._define_default_params()
        
    def _define_default_params(self):
        """
        make sure to define any default parameters that may be needed later
        but aren't necessarily specified in the flash.par file
        """
        
        defaults = {'gamma': 1.4}
        
        for d in defaults:
            if not d in self.params.keys():
                self.params[d] = defaults[d]
        
        
        
    def param_contained(self, field, weight = None, calc_type = 'average',
                              r = None):
        """
        computes the average or total quantity contained within some radius.
        The default is to compute an average or total quantity, whichever makes,
        sense, and to use the dwarf radius. Default is unweighted.
        """
        
        if r == None:
            r = self.radius
        
        sp = self.ds.sphere( self.center, r)
        
        if field == 'total_mass':
            return sp.quantities.total_mass()
        
        
        if calc_type == 'total':
            result = sp.quantities.total_quantity(field)
        elif calc_type == 'average':
        
            result = sp.quantities.weighted_average_quantity(field, weight)
            
            
        return result


    def rmatch(self, nProfile, eps = 1.0E-7, nmax = 2000):
        """
        Calculates the match radius of the dwarf galaxy 
        This is just a copy paste of the Fortran code (almost) in FLASH sim
        """
        G = 6.67259E-8 ! grav constant in cgs
        k = 1.380658E-16
        mh = 1.6733E-24
        
        gamma = self.params['gamma']
        
        cs1 = (gamma * k * self.params['sim_TCloud'] / mh)
        cs2 = (gamma * k * self.params['sim_TAmbient' / mh)
        
        cPhi = 4.0 * np.pi * G * self.params['sim_rhoCenter'] *\
                    self.params['sim_bParam'] ** 3
                    
        cRho2 = self.params['sim_rhoRL'] * np.exp((-cPhi/\
                     (cs2*cs2*self.params['sim_RL']))*\
                     np.log(1.0 + self.params['sim_RL'] /\
                     self.params['sim_bParam']))
                     
        ###
        rhi = 2.0 * self.params['sim_yCenter']
        rlo = 0.0
        ic  = 0
        
        rhomid = 1.0E-32
        
        
        while (((np.abs(rhomid - self.params['sim_rho2rm'])\
                /self.params['sim_rho2rm'] > eps) and (ic <= nmax))):
                
            rmid = 0.5 * (rhi + rlo)
            rhomid = cRho2 * np.exp( (cPhi/(cs2*cs2*rmid))*\
                                    np.log(1.0+rmid/self.params['sim_bParam']))
                                    
            if (rhomid .gt. self.params['sim_rho2rm']):
                rlo = rmid
            else:
                rhi = rmid
            
        
        self.RM = rmid
        
        return self.RM
                



def rps_param():
    """ 
    Calculate the ram pressure parameter
    """
    rps = 20
    
    return rps
    
        
        
    
        



        
