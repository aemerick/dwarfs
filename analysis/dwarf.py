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
        
