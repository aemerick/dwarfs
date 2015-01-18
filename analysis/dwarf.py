import os
import yaml
import numpy as np
import yt
import glob

from yt import units as u

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper



class simulation: # need a better name

    def __init__(self, ds_prefix, param_file = "flash.par", ds_dir="./"):
        """
        Initiate simulation class

        ds_prefix : string
            Prefix name to checkpoint and particle files

        param_file : string
            Name of parameter file (assumed "flash.par")

        ds_dir : string
            Directory where all the data is located (assumed "./")
        """

        self.param_file = param_file
        self.ds_dir     = ds_dir

        # make the ds_prefix prettier and just the meat        
        if "chk" in ds_prefix:
            ds_prefix = ds_prefix.replace("chk","")
            if "__" in ds_prefix:
                ds_prefix = ds_prefix.replace("__","")

        if ds_prefix[-1] == "_":
            ds_prefix = ds_prefix[:-1]
 
        self.ds_prefix = ds_prefix

        # get the checkpoint file names into a list
        ds_list = glob.glob(ds_dir + ds_prefix + "*chk*")
        ds_list.sort()
        self.ds_list = ds_list
    
        # check if there are particle files. Do the same
        # if there are         
        part_list = glob.glob(ds_dir + ds_prefix + "*part*")
        part_list.sort()
        self.part_list = part_list 

        # load SN and SB files
        self._load_SN_data()
        self._load_SB_data()
        
        # load the flash.par parameter file
        self.params = _load_param_file(self)

        # compute (roughly) t at each checkpoint file
        self.times = _get_ds_index()\
                     * self.params['checkpointFileIntervalTime']
  
    def _get_ds_index(self):
        """
        Returns the integer number associated will all 
        known checkpoint files as an integer np array
        """
    
        values = np.zeros(np.size(self.ds_list))
        i = 0
        for dsname in ds_list:
            values[i] = int(dsname[-4:])
            i = i + 1

        return values
          
        
    def _load_param_file(self):
    
        newfile = self.ds_dir + self.param_file + ".mod"
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
    
        stream = open(newfile, 'r')
        self.params     = yaml.load(stream, Loader=Loader)
        self._define_param_units()

        self.center = np.array( [ np.float(self.params['sim_xctr']),
                                  np.float(self.params['sim_yctr']),
                                  np.float(self.params['sim_zctr']) ])


    def _load_SN_data(self):
  
        print "this function will load supernova info if it exists"

        sn_path = self.ds_dir + "SNfeedback.dat"

        if not os.path.isfile(sn_path):
            self.SNdata = None
            print "No supernova feedback file found at " + sn_path
            return
 
        # load the file here
        #SN_data = np.genfromtxt(sn_path, names=True)
        self.SN = SN(sn_path)
        


    def _load_SB_data(self):

        sb_path = self.ds_dir + "SBfeedback.dat"

        if not os.path.isfile(sb_path):
            self.SBdata = None
            print "No SB feedback file found at " + sb_path
            return

        # load superbubble data
        self.SB = SB(sb_path)     

        
        
        

    def _define_param_units(self):

        print "currently does nothing"
    
    
class SNSB:

    def __init__(self, file_path):
        self.fname = file_path

        data = np.genfromtxt(file_path, names=True)

        # do some bookkeeping on header names to make them
        # nicer. Currently, headers are 00xxx, 01yyy, 02zzz, etc.
        # want to remove the "00", "01", "02" ...
        i = 0
        names_list = list(data.dtype.names)
        for name in names_list:
            num_str = "%02i"%(i)           
            names_list[i] = name.replace(num_str,'')
            i = i + 1
        
        # now, reassign
        data.dtype.names = names_list

        # save
        self.data = data
        
    def plot_positions(self):
        print "does nothing"

    def _set_units(self):
        self.data['time']   *=  u.s
        self.data['posx']   *= u.cm
        self.data['posy']   *= u.cm
        self.data['posz']   *= u.cm
#        self.data['radius'] *= u.cm
#        self.data['mass']   *=  u.g

    



class SN(SNSB):
 
    def __init__(self, file_path):

        SNSB.__init__(self, file_path)

        self._set_units()

    def _set_units(self):
        """
        sets yt units for the data values that have units. 
        calls parent function for shared values 
        """
 #       self.data['time']   *=  u.s
 #       self.data['posx']   *= u.cm
 #       self.data['posy']   *= u.cm
 #       self.data['posz']   *= u.cm
        SNSB._set_units(self)
        self.data['radius'] *= u.cm
        self.data['mass']   *=  u.g
    


class SB(SNSB):

    def __init__(self,file_path):

        SNSB.__init__(self, file_path)

        self._set_units()

        self._load_creation_data()

    def _set_units(self):
        """
        """
        SNSB._set_units()
        
        self.data['velx'] *= u.cm / u.s
        self.data['vely'] *= u.cm / u.s
        self.data['velz'] *= u.cm / u.s

    def _load_creation_data(self):
        """
        """
      
        sb_create_path = self.fname.replace('feedback.dat','create.dat')
    
        if (not os.path.isfile(sb_create_path)):
            self.creation_data = None
            print "No SBcreate.dat found at " + sb_create_path
            return

        data = np.genfromtxt(sb_create_path, names=True)

        # adjust header names
        i = 0
        names_list = list(data.dtype.names)
        for name in names_list:
            num_str = "%02i"%(i)
            names_list[i] = name.replace(num_str,'')
            i = i + 1

        # now, reassign
        data.dtype.names = names_list

        data['time'] *= u.s
        data['posx'] *= u.cm; data['poxy'] *= u.cm; data['posz'] *= u.cm
        data['velx'] *= u.cm/u.s ; data['vely'] *= u.cm / u.s; 
        data['velz'] *= u.cm/u.s
        data['SNinvT'] *= u.s

        self.creation_data = data
   
           


class dwarf:

    def __init__(self, ds, param_file, raw = True, rm=None):
        """
        If raw = true, do some grooming first to import the parameter
        file
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
        self.time = self.ds.current_time.value * yt.units.s
        
        self.center = np.array( [ np.float(self.params['sim_xctr']),
                                  np.float(self.params['sim_yctr']),
                                  np.float(self.params['sim_zctr']) ])

        self.center = self.center    
                                       
        self.radius = np.float( self.params['sim_RL'] )

        if not rm == None:
            self.RM = rm

        
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
                if not self.RM == None:
                    r = self.RM
                else:
                    self.rmatch(2000)
                    r = self.RM
            else:
                self.rmatch(2000)
                r = self.RM
                
                
            data_source = self.ds.sphere(self.center, (2.0*r,'cm'))
        else:
            r = data_source.radius

       
        if xmin == None or xmax == None: 
            if xfield == 'radius':
                 xmin = 0.0 * yt.units.cm
                 xmax = 2.0 * r
            
            else:
                xmin = np.min(sphere[xfield])
                xmax = np.max(sphere[xfield])
                
        print xmin, xmax, 'asdfasdfa'    

#yt.create_profile(sphere, 'radius', ['Hot_Gas_Mass'],
 #                               n_bins=nbins, extrema = {'radius': (0.0, R)},un$
  #                              weight_field = None, logs={'radius':False})

        prof = yt.create_profile(data_source,xfield, [field], n_bins=nbin,
                          extrema={xfield:(xmin,xmax)},units={'radius':'cm'},
                          weight_field = weight, logs={'radius':False})

        
#        prof = yt.Profile1D(data_source, xfield, nbin, xmin, xmax, weight)
            
 #       prof.add_fields(field)
  
       
        return prof.x_bins, prof[field]



    


def rps_param():
    """ 
    Calculate the ram pressure parameter
    """
    rps = 20
    
    return rps
    
        
        
    
        



        
