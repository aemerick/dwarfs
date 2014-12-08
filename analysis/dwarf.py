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
        self._define_param_units()
        
        
        self.center = np.array( [ np.float(self.params['sim_xctr']),
                                  np.float(self.params['sim_yctr']),
                                  np.float(self.params['sim_zctr']) ])

        self.center = self.center    
                                       
        self.radius = np.float( self.params['sim_RL'] )
        

        
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
    
    
    def _define_param_units(self):
        """
        Damnit yt
        """
        
        length_unit = yt.units.cm
        temp_unit   = yt.units.Kelvin
        time_unit   = yt.units.s
        mass_unit   = yt.units.g
        speed_unit  = yt.units.km / time_unit
        pressure_unit = mass_unit / length_unit / time_unit**2
        density_unit  = mass_unit / length_unit**3
        
        length_params = ['sim_xctr','sim_yctr', 'sim_zctr', 'sim_RL',
                         'sim_bParam', 'sim_wTaper', 'sim_rScale',
                         'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
   
   
        density_params = ['sim_smallRho', 'sim_rhoAmbient', 'sim_rhoCloud',
                          'sim_rhoCenter', 'sim_rho1rm', 'sim_rho2rm',
                          'sim_rhoRL']
   
   
        temperature_params = ['sim_TAmbient', 'sim_TCloud']

        pressure_params    = ['sim_smallP', 'sim_pAmbient']
         
        speed_params       = ['sim_windVel', 'sim_flowSpeed']
   
        param_dict = {'length': length_params, 'density': density_params,
                      'pressure': pressure_params,
                      'temperature': temperature_params,
                      'speed'      : speed_params}
                      
        units_dict = {'length': length_unit, 
                      'density': density_unit,
                      'temperature': temp_unit,
                      'speed': speed_unit,
                      'pressure': pressure_unit} 
                         
        for ptype in param_dict:
        
            # would like to do without this loop, but is really the only way
            # to be save with the exceptions... 
            for pname in param_dict[ptype]:
            
                try:
                    self.params[pname] = np.float(self.params[pname])\
                                                             * units_dict[ptype]
                except KeyError:
                    pass


                
        
        
        
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
        G = 6.67259E-8    * yt.units.cm**3 / yt.units.g / yt.units.s**2
        k = 1.380658E-16  * yt.units.erg / yt.units.Kelvin
        mh = 1.6733E-24   * yt.units.g
        
        gamma = self.params['gamma']

        cs1 = (gamma * k * self.params['sim_TCloud'] / mh) ** 0.5
        cs2 = (gamma * k * self.params['sim_TAmbient'] / mh) ** 0.5
        
        cPhi = 4.0 * np.pi * G * self.params['sim_rhoCenter'] *\
                    self.params['sim_bParam'] ** 3
                    
        cRho2 = self.params['sim_rhoRL'] * np.exp((-cPhi/\
                     (cs2*cs2*self.params['sim_RL']))*\
                     np.log(1.0 + self.params['sim_RL'] /\
                     self.params['sim_bParam']))
                     
        ###
        rhi = 2.0 * self.params['sim_yctr']
        rlo = 0.0 * yt.units.cm
        ic  = 0
        
        rhomid = 1.0E-32 * yt.units.g / yt.units.cm**3
        
        
        while (((np.abs(rhomid - self.params['sim_rho2rm'])\
                /self.params['sim_rho2rm'] > eps) and (ic <= nmax))):
                
            rmid = 0.5 * (rhi + rlo)
            rhomid = cRho2 * np.exp( (cPhi/(cs2*cs2*rmid))*\
                                    np.log(1.0+rmid/self.params['sim_bParam']))
                                    
            if (rhomid > self.params['sim_rho2rm']):
                rlo = rmid
            else:
                rhi = rmid
            
            ic = ic + 1
            
        self.RM = rmid 
        
        return self.RM
                
    def profile(self, field, nbin = 10, weight = None, data_source = None,
                      xfield='radius', xmin=None, xmax=None):
        """
        If data source is none, computes profile from rmatch
        """

        if data_source == None:
            if hasattr(self, 'RM'):
                r = self.RM
            else:
                self.rmatch(2000)
                r = self.RM
                
                
            data_source = self.ds.sphere(self.center, (2.0*r,'cm'))
       
        if xmin == None or xmax == None: 
            if xfield == 'radius':
                 xmin = 0.0 * yt.units.cm
                 xmax = 2.0 * r
            
            else:
                xmin = np.min(sphere[xfield])
                xmax = np.max(sphere[xfield])
                
            
        
        prof = yt.Profile1D(data_source, xfield, nbin, xmin, xmax, weight)
            
        prof.add_fields(field)
        
        return prof.x, prof.field_data[field]



def rps_param():
    """ 
    Calculate the ram pressure parameter
    """
    rps = 20
    
    return rps
    
        
        
    
        



        
