import os
import yaml
import numpy as np
import yt
import glob

from yt import units as u
from plotting import plotTools as myplot # some convenience plotting tools

from initial_conditions import ic_generator as ic

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


# some functions for analysis
from initial_conditions import profiles as prof


class simulation: # need a better name

    def __init__(self, ds_prefix, param_file = "flash.par", ds_dir="./",
                       exact_times=True):
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
        # try to guess if user wants plot files or checkpoint files
        if "chk" in ds_prefix:
            ds_prefix = ds_prefix.replace("chk","")
            if "__" in ds_prefix:
                ds_prefix = ds_prefix.replace("__","")

        if "plt" in ds_prefix:
            ds_prefix = ds_prefix.replace("plt","")
            if "__" in ds_prefix:
                ds_prefix = ds_prefix.replace("__","")
        
        if "cnt" in ds_prefix:
            ds_prefix = ds_prefix.replace("cnt","")
            if "__" in ds_prefix:
                ds_prefix = ds_prefix.replace("__","")


        # remove any trailing underscores
        if ds_prefix[-1] == "_":
            ds_prefix = ds_prefix[:-1]
 
        self.ds_prefix = ds_prefix

        # get the checkpoint file names into a list
        chk_list = glob.glob(ds_dir + ds_prefix + "*chk*")
        chk_list.sort()
        self.chk_list = chk_list
        
        plt_list = glob.glob(ds_dir + ds_prefix + "*plt*")
        plt_list.sort()
        self.plt_list = plt_list
    
        # check if there are particle files. Do the same
        # if there are         
        part_list = glob.glob(ds_dir + ds_prefix + "*part*")
        part_list.sort()
        self.part_list = part_list 


        self._load_param_file()
        self.center = np.array( [ np.float(self.params['sim_xctr']),
                                  np.float(self.params['sim_yctr']),
                                  np.float(self.params['sim_zctr']) ])
        self.center = self.center * yt.units.cm


        self._load_SN_data()
        self._load_SB_data()

        # define the radius using the cloud temperature
        self._find_radius()

        

        # compute (roughly) t at each checkpoint file

#        print self._get_ds_index()
#        print self.params['checkpointFileIntervalTime']

        self._load_times(exact_times)

    def _find_radius(self):
        """
        Finds an estimate for the initial radius of the dwarf galaxy
        as the farthest point away from the center of the box that
        has the dwarf initial temperature
        """
        ds = yt.load(self.plt_list[0])
        data = ds.all_data()        

        x = data['x'].convert_to_units('cm')
        y = data['y'].convert_to_units('cm')
        z = data['z'].convert_to_units('cm')
        r = self.dist_from_center(x,y,z)

        T = data['temp'].convert_to_units('K')
        tcloud = self.params['sim_TCloud']

        rdwarf = r[np.abs(T-tcloud)/tcloud < 0.1]
        rdwarf = rdwarf.convert_to_units('pc')

        self.radius = np.max(rdwarf)

        return

    def evaluate_potential(self, r):
        """
        Uses potential given in profiles to calculate the static 
        analytic potential value at some point r
        """

        functions_lookup = {'3': prof.NFW_potential,
                            '4': prof.NFW_potential,
                            '5': prof.Burkert_potential}

        function = functions_lookup[str(self.params['density_profile'])]

        # load params and sanitize units
        r_s  = self.params['sim_bparam'].value
        M200 = self.params['sim_M200'].value
        rho_crit = self.params['sim_rho_crit'].value
        r = r.convert_to_units('cm').value

        return function(r, r_s=r_s, M200=M200, rho_crit=rho_crit) * yt.units.cm**2 / yt.units.s**2


    def dist_from_center(self, x, y, z):
        """ 
        Returns radial distance of cartesian coordinate values
        from simulation center
        """
        
        return ((x - self.center[0])**2 + (y - self.center[1])**2 +\
               (z - self.center[2])**2 )**0.5
     
        

    def _load_times(self, exact):
        """
        Roughly calculates the times of every loaded file based
        upon the parameter file interval time. THIS WILL NOT 
        BE COMPLETELY ACCURATE IF dt DURING SIMULATION IS 
        LARGE.
        """
#        print "approximating time stamps from .par file"

        self.times = {}

        if exact == False:
            print "approximating time stamps from .par file"
            self.times['plt'] = self._get_ds_index(ds_type='plt')*\
                            self.params['plotfileIntervalTime']
            self.times['chk'] = self._get_ds_index(ds_type ='chk')*\
                            self.params['checkpointFileIntervalTime']
        
            
        else:
            # do this exactly by loading every file and reading timestamp
            self.times['plt'] = np.zeros(np.size(self.plt_list))
            self.times['chk'] = np.zeros(np.size(self.chk_list))
 
            i = 0
            for plt in self.plt_list:
                ds = yt.load(plt)
                self.times['plt'][i] = ds.current_time.value
                i = i + 1
            i = 0
            for chk in self.chk_list:
                ds  = yt.load(chk)
                self.times['chk'][i] = ds.current_time.value
                i = i + 1
      
        # convert to Myr
        self.times['plt'] = (self.times['plt'] * yt.units.s).convert_to_units('Myr')
        self.times['chk'] = (self.times['chk'] * yt.units.s).convert_to_units('Myr')
        return

    def _get_ds_index(self, ds_type='plt'):
        """
        Returns the integer number associated will all 
        known plt / checkpoint files as an integer np array
        """
        if ds_type == 'plt':
            ds_list = self.plt_list
        else:
            ds_list = self.chk_list    


        if np.size(self.ds_list) == 0:
            values = np.array([0.0])
        else:

            values = np.zeros(np.size(self.ds_list))
            i = 0
            for dsname in self.ds_list:
                values[i] = int(dsname[-4:])
                i = i + 1

        return np.array(values)
          
        
    def _load_param_file(self):
    
        newfile = self.ds_dir + self.param_file + ".mod"
        # os.system("cp " + param_file + " " + newfile)

        # convert = to : and save to new file
        bash_command = "sed 's/=/:/g' " + self.param_file + " > " + newfile
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

    def _define_param_units(self):
        """
        Damnit yt..... loads in all probably relevant parameter file
        names and assigns units to make working with yt things a little
        bit easier...
        """

       # print "in define param units function"

        length_unit = yt.units.cm
        temp_unit   = yt.units.Kelvin
        time_unit   = yt.units.s
        mass_unit   = yt.units.g
        speed_unit  = yt.units.km / time_unit
        pressure_unit = mass_unit / length_unit / time_unit**2
        density_unit  = mass_unit / length_unit**3

        length_params = ['sim_xctr','sim_yctr', 'sim_zctr', 'sim_RL',
                         'sim_bparam', 'sim_wTaper', 'sim_rScale',
                         'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']


        density_params = ['sim_smallRho', 'sim_rhoAmbient', 'sim_rhoCloud',
                          'sim_rhoCenter', 'sim_rho1rm', 'sim_rho2rm',
                          'sim_rhoRL','sim_rho_crit']

        time_params   = ['checkpointFileIntervalTime',
                         'particleFileIntervalTime',
                         'plotfileIntervalTime',
                         'tmax', 'dtmin', 'dtmax',
                         'sblife','tsn1','tsn2','tsb']

        temperature_params = ['sim_TAmbient', 'sim_TCloud']
        pressure_params    = ['sim_smallP', 'sim_pAmbient']

        speed_params       = ['sim_windVel', 'sim_flowSpeed']

        mass_params        = ['sim_M200']

        param_dict = {'length': length_params, 'density': density_params,
                      'pressure': pressure_params,
                      'temperature': temperature_params,
                      'speed'      : speed_params, 'time' : time_params,
                      'mass': mass_params}

        units_dict = {'length': length_unit,
                      'density': density_unit,
                      'temperature': temp_unit,
                      'speed': speed_unit,
                      'pressure': pressure_unit, 'time' : time_unit,
                      'mass'    : mass_unit}

               
        # now that everything is defined, load and assign all units

        for ptype in param_dict:

            # would like to do without this loop, but is really the only way
            # to be safe with catching exception if param is not in flash.par
            for pname in param_dict[ptype]:
                try:
                   self.params[pname] = np.float(self.params[pname])\
                                                             * units_dict[ptype]
                except KeyError:
                    print "Did not find parameter: " + pname    
                    pass
               


    def _load_SN_data(self):
  
        print "this function will load supernova info if it exists"

        sn_path = self.ds_dir + "SNfeedback.dat"

        if not os.path.isfile(sn_path):
            self.SNdata = None
            print "No supernova feedback file found at " + sn_path
            return
 
        # load the file here
        #SN_data = np.genfromtxt(sn_path, names=True)
        self.SN = SN(sn_path, center = self.center)
        


    def _load_SB_data(self):

        sb_path = self.ds_dir + "SBfeedback.dat"

        if not os.path.isfile(sb_path):
            self.SBdata = None
            print "No SB feedback file found at " + sb_path
            return

        # load superbubble data
        self.SB = SB(sb_path, center = self.center)     

        
    def get_initial_conditions(self, r=[]):
        """
        Wrapper to ic generator to get initial density, pressure,
        and temperature profiles
        """
   
        if np.size(r) == 0:
            rmin = 0.0 
            rmax = (2000.0 * u.kpc).convert_to_units('cm').value
            r = np.linspace(rmin, rmax, 1.0E4)
            
       
        
        RM, rho, pressure, T =  ic.spherical_NFW(r, self.params['sim_TCloud'].value,
                               self.params['sim_TAmbient'].value,
                               self.params['sim_bParam'].value,
                               #nself.params['mu'],
                               self.params['sim_RL'].value,
                               self.params['sim_rhoRL'].value,
                               self.params['sim_rhoCenter'].value,
                               self.params['sim_rho1rm'].value,
                               self.params['sim_rho2rm'].value)
                           
        return RM*u.cm, r * u.cm, (rho * u.g / (u.cm**3)),\
               pressure * (u.g/(u.cm*u.s*u.s)), T * u.Kelvin

    
    
    
class SNSB:
    """
    General class to read in SN and SB feedback data and do some 
    analysis with them. This is intended to be a parent class
    to the SN and SB classes, that will do things that are general / 
    applicable to both (since they are similar). Many of these functions
    will have overwritten version in the SN and SB classes that will
    be called from there. There should really be no point in 
    initializing just a SNSB object. 
    """


    def __init__(self, file_path, center = np.zeros(3)):
        """
        Initializes SNSB feedback data files. Does general things
        for both, but leaves details to separate SN and SB functions
        
        """
        self.fname = file_path

        self.center = center
 
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
        
    def plot_positions(self, **kwargs):
        """
        returns a 3d matplotlib plot of the 3d positions
        of the 
        """
        plt.close() # clear  

        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d')

        posx, posy, posz = self._recenter()
        ax.scatter(posx.value, posy.value, posz.value, kwargs)

        return fig, ax

    def _recenter_pos(self):
        """
        returns coordinates in the frame of sim center (the dwarf)
        """
  
        return self.data['posx'] - self.center[0],\
               self.data['posy'] - self.center[1],\
               self.data['posz'] - self.center[2]

    def _set_units(self):
        self.data['time']   *=  u.s
        self.data['posx']   *= u.cm
        self.data['posy']   *= u.cm
        self.data['posz']   *= u.cm
    



class SN(SNSB):
    """
        Derived class of the supernova and superbubble feedback
        class (SNSB). Contains functions specific to the supernova
        feedback data.
    """

 
    def __init__(self, file_path, center = np.zeros(3)):

        SNSB.__init__(self, file_path, center)

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
    

    def plot_positions(self, draw_radius = False, **kwargs):

        fig, ax = SNSB.plot_positions(self,kwargs)
 
        if draw_radius:

            sn_spheres = []

            posx, posy, posz = self._recenter()
            i = 0
            for r in self.data['radius']:
                x,y,z = posx[i], posy[i], posz[i]
                myplot.draw_sphere(ax, r.value,
                        center = np.array([x.value,y.value,z.value]))
                i = i + 1
        
            return fig, ax, sn_spheres

        else:

            return fig, ax
            

class SB(SNSB):
    """
       Derived class of supernova and superbubble feedback class.
       This contains functions and things specifically for the 
       superbubble data
    """


    def __init__(self,file_path, center = np.zeros(3)):
        
 
        SNSB.__init__(self, file_path, center)

        self._set_units()

        self._load_creation_data()

    def _set_units(self):
        """
        """
        SNSB._set_units(self)
        
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
    """
        This is an old class that is here for backwards compatability 
        sakes until I update everything. The replacement will be the 
        simulation class, which is not tied to a single output file,
        but rather the simulation in the abstract.
    """

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
                                       
        #self.radius = np.float( self.params['sim_RL'] )

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



def dwarf_mass(sim, out_file, tmin=None, tmax=None, mode='grav', T_range=[]):
    """
       Calculate the mass profile of the dwarf as a function of 
       time. Can do this by either mass contained in the radius or
       computing all gravitationally bound gas in the box
    """
    
    ds_list = sim.plt_list
    
    if tmin == None:
        tmin = 0.0 * yt.units.Myr
    if tmax == None:
        tmax = 1.0E4 * yt.units.Myr

    if len(T_range) == 2:
        T_range = np.array(T_range)*yt.units.Kelvin

    
    sim.dwarf_mass = {}
    sim.dwarf_mass[mode] = np.ones(np.size(sim.times['plt']))*-1

    
    # do only over select time range
    ds_min  = np.argmin(np.abs((sim.times['plt'] - tmin).value))
    ds_max  = np.argmin(np.abs((sim.times['plt'] - tmax).value))


    file = open(out_file, 'w')
    file.write("# t m\n")
    format = "%8.8e %8.8e\n"
    print mode
    if mode == 'grav':
        print 'calculating'
        # calculate dwarf gas mass as total bound gas mass of dwarf
        # 1) --- do over whole box! -- 
        # 2) get thermal energy and kinetic energy in each cell.. sum..
        # 3) get position in each cell (r) and compute the value 
        #    of the potential
        # 4) get mass in each cell
        # 5) sum cell mass where   E_th + E_kin - m*Phi(r) < 0
 
        print ds_min, ds_max
        i = 0
        for dsname in ds_list[ds_min:ds_max]:
            ds = yt.load(dsname); data = ds.all_data()
            x = data['x'].convert_to_units('cm')
            y = data['y'].convert_to_units('cm')
            z = data['z'].convert_to_units('cm')  
 
            mass      = data['dens'] * data['dx'] * data['dy'] * data['dz']        
            mass = mass.convert_to_units('g')

            E_kin = 0.5 * mass * (data['velx']**2 + data['vely']**2 + data['velz']**2)
            E_kin = E_kin.convert_to_units('erg')

            r = sim.dist_from_center(x,y,z)

            phi       = sim.evaluate_potential(r)
            U_grav    = mass * phi

            if len(T_range) == 2:
                T = data['temp'].convert_to_units('K')
                total_mass = np.sum( mass[(E_kin<U_grav)*(T>T_range[0])*(T<T_range[1])] )
           
            else:
                total_mass = np.sum(mass[(E_kin<U_grav)])

            total_mass = total_mass.convert_to_units('Msun')
           

            sim.dwarf_mass[mode][i + ds_min] = total_mass
            file.write(format%(ds.current_time.convert_to_units('Myr'),total_mass))
            i = i + 1
    
    elif mode == 'contained':
        # calculate the total dwarf mass as just the total mass contained
        # within the initial dwarf radius (with optional temperature cuts)
        i = 0
        for dsname in ds_list[ds_min:ds_max]:
            ds = yt.load(dsname); data = ds.all_data()

            sp = ds.sphere(sim.center,sim.radius)
            mass = sp['dens'] * sp['dx'] * sp['dy'] * sp['dz']
            
            if len(T_range) == 2:
                T = sp['temp'].convert_to_units('K')
                total_mass = np.sum(mass[(T >= T_range[0])*(T<=T_range[1])])


            else:
                total_mass = np.sum(mass)

            sim.dwarf_mass[mode][i + ds_min] = total_mass
            file.write(format%(ds.current_time.convert_to_units('Myr').value,total_mass.value))
            i = i + 1

    file.close()
    return 
        



        
