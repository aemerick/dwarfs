#"""
#dwarf
#
#    Description: Dwarf Galaxy Simulation Analysis
#               - Purpose of code is to develop object oriented
#                 method to quickly and easily analyze large sets
#                 of data from FLASH simulations of dwarf galaxy
#                 ram pressure stripping. Relies substantially on
#                 yt 3.x
#    Created : Nov. 2014
#    Author  : Andrew Emerick 
#    Affil   : Columbia University - American Museum of Natural History
#    Contact : emerick@astro.columbia.edu
#

import os
import yaml
import numpy as np
import yt
import glob

import matplotlib.pyplot as plt

from yt import units as u
from plotting import plotTools as myplot # some convenience plotting tools

from initial_conditions import ic_generator as ic

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


# now import some of my own functions 
from initial_conditions import profiles as prof
import cgs as cgs
#from dwarf_yt_functions import derived_fields




class simulation: # need a better name
    """
    The simulation class is meant to be an easy way to abstract some analysis
    and interface with the entire dwarf galaxy simulation as a single object.
    Rather than loading a slew of individual plot / checkpoint files into a 
    list, for example, loading a simulation class will automatically look
    for all data outputs, the corresponding flash.par file, and any
    ancillary outputs like SNfeedback.dat, SBfeedback.dat. For the data files,
    it stores the filepaths to all to make easy looping over the data sets.
    It also checks through all outputs and records the time of each output.
   
    This is convenient and abtracts some analysis, so one could just load
    the simulation object and pass it to a function, and do the same operation
    over all contained data sets fairly easily. For example, calculating the
    bound mass of the dwarf galaxy over time can be done in two lines:

    sim = dwarf.simulation('data_prefix_')
    dwarf.dwarf_mass(sim, 'outputfile.dat')

    A.E. 6/1/15:
    A lot of code is written with no current obvious usefullness (such as
    loading in some of the supernova data, for example). This will be 
    changed ultimately as work progresses from FLASH coding / simulation to
    analysis.
    """


    def __init__(self, ds_prefix, param_file = "flash.par", ds_dir="./",
                       exact_times=True, reload_times=False):
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

        # load the simulation times 
        filename = ds_dir + ds_prefix + "times.dat"
        self._load_times(exact_times, filename, reload_times)

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
     
        

    def _load_times(self, exact, filename, reload_times):
        """
        Roughly calculates the times of every loaded file based
        upon the parameter file interval time. THIS WILL NOT 
        BE COMPLETELY ACCURATE IF dt DURING SIMULATION IS 
        LARGE.
        """

        self.times = {}

        if exact == False:
            _myprint("WARNING - approximating time stamps from .par file")
            self.times['plt'] = self._get_ds_index(ds_type='plt')*\
                            self.params['plotfileIntervalTime']
            self.times['chk'] = self._get_ds_index(ds_type ='chk')*\
                            self.params['checkpointFileIntervalTime']
        
            
        else:
            
            if reload_times or not os.path.isfile(filename):
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
                self.times['plt'] = (self.times['plt'] *\
                                     yt.units.s).convert_to_units('Myr')
                self.times['chk'] = (self.times['chk'] *\
                                     yt.units.s).convert_to_units('Myr')
                
                f = open(filename,'w')
                f.write("#time\n")
                for time in self.times['plt']:
                    f.write("%8.8e\n"%(time.value))
                f.close()
            else:
                data = np.genfromtxt(filename,names=True)
                self.times['plt'] = data['time']*yt.units.Myr
      
        # convert to Myr
       
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
        oldfile = self.ds_dir + self.param_file
        
        # convert = to : and save to new file
        bash_command = "sed 's/=/:/g' " + oldfile+ " > " + newfile
        os.system(bash_command)

        # remove all comments
        bash_command = "sed -i 's:#.*$::g' " + newfile
        os.system(bash_command)

        # remove all blank lines
        bash_command = "sed -i '/^$/d' " + newfile
        os.system(bash_command)
    
        # load the params from the new file and define units
        stream = open(newfile, 'r')
        self.params     = yaml.load(stream, Loader=Loader)
        self._define_param_units()

    def _define_param_units(self):
        """
        Damnit yt..... loads in all probably relevant parameter file
        names and assigns units to make working with yt things a little
        bit easier...
        """

        # for cleanliness, compile unit types first
        length_unit = yt.units.cm
        temp_unit   = yt.units.Kelvin
        time_unit   = yt.units.s
        mass_unit   = yt.units.g
        speed_unit  = yt.units.km / time_unit
        pressure_unit = mass_unit / length_unit / time_unit**2
        density_unit  = mass_unit / length_unit**3

        # make lists of parameters to assign units to, grouped by unit type
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

        # dicts to associate units with above lists of variables
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
                    _myprint("[INIT] Did not find parameter: " + pname)
                    pass
               


    def _load_SN_data(self):
  
        _myprint("[INIT] Looking for supernova files")

        sn_path = self.ds_dir + "SNfeedback.dat"

        if not os.path.isfile(sn_path):
            self.SNdata = None
            _myprint("[INIT] No supernova feedback file found at " + sn_path)
            return
 
        # load the file here.. use the supernova class
        self.SN = SN(sn_path, center = self.center)
        

    def _load_SB_data(self):

        sb_path = self.ds_dir + "SBfeedback.dat"

        if not os.path.isfile(sb_path):
            self.SBdata = None
            _myprint("[INIT] No SB feedback file found at " + sb_path)
            return

        # load superbubble data using superbubble class
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

#
# -----    end simulation class
    
    
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

        posx, posy, posz = self._recenter_pos()
        ax.scatter(posx, posy, posz, **kwargs)

        return fig, ax

    def _recenter_pos(self):
        """
        returns coordinates in the frame of sim center (the dwarf)
        """
  
        return self.data['posx'] - self.center[0].value,\
               self.data['posy'] - self.center[1].value,\
               self.data['posz'] - self.center[2].value

    def _set_units(self):
        self.data['time']   = self.data['time']*  u.s
        self.data['posx']   = self.data['posx']*  u.cm
        self.data['posy']   = self.data['posy']*  u.cm
        self.data['posz']   = self.data['posz']*  u.cm
    

# end SNSB general class

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

        SNSB._set_units(self)
        _myprint( "trying to set units in supernova")
        self.data['radius'] = self.data['radius']* yt.units.cm
        self.data['mass']   = self.data['mass']  * yt.units.g
    

    def plot_positions(self, draw_radius = False, **kwargs):

        fig, ax = SNSB.plot_positions(self, **kwargs)
 
        if draw_radius:

            sn_spheres = []

            posx, posy, posz = self._recenter_pos()
            i = 0
            for r in self.data['radius']:
                x,y,z = posx[i], posy[i], posz[i]
                myplot.draw_sphere(ax, r.value,
                        center = np.array([x.value,y.value,z.value]))
                i = i + 1
        
            return fig, ax, sn_spheres

        else:

            return fig, ax

    
        #
    # end plot positions
    def rate(self, sntype = 2, units="Myr"):
        """
            Calculates and returns the supernova rate.

            Assumes rate for type II supernova. 

            sntype : int
               Type of supernova rate to calculate. Default is 2
               0 = Total Rate (all SN)
               1 = Type 1a
               2 = Type II
               3 = Type II from SB
        """
        # do a numerical derivative
        total_SN = np.size(self.data['time'])        

        if sntype == 0:
            type_select = np.arange(0, total_SN)
        else:
            type_select = self.data['type'] == sntype

        

        t = self.data['time'][type_select]
        N = np.arange(1, total_SN + 1)[type_select]

        dNdt = np.zeros(t.shape, np.float)

        dNdt[1:-1] = (N[2:] - N[0:-2])/(t[2:] - t[0:-2])
        dNdt[0]    = (N[1]  - N[0]   )/(t[1]  - t[0]   )
        dNdt[1]    = (N[-1] - N[-2]  )/(t[-1] - t[-2]  )

        

        return t, dNdt
  
    def rate_from_gas(units="1/Myr"):
        """
        Calculates the total supernova rate from
        """

        return 'does nothing right now'



# - end SN class

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
        
        self.data['velx'] = self.data['velx'] * u.cm / u.s
        self.data['vely'] = self.data['vely'] * u.cm / u.s
        self.data['velz'] = self.data['velz'] * u.cm / u.s

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

        data['time'] = data['time'] * u.s
        data['posx'] = data['posx'] * u.cm
        data['posy'] = data['posy'] * u.cm
        data['posz'] = data['posz'] * u.cm
        data['velx'] = data['velx'] * u.cm/u.s
        data['vely'] = data['vely'] * u.cm / u.s; 
        data['velz'] = data['velz'] * u.cm/u.s
        data['SNinvT'] = data['SNinvT'] * u.s

        self.creation_data = data
        #
    #- end load creation data
# - end SB  class
   
## ------------------------------------------------------------           

def _select_ds_list(sim, tmin, tmax, ftype='plt', dt=None):
    """
    Given time range, get list of plt files to operate over. This 
    needs to be improved upon by quite a bit.
    """

    # set tmin/tmax if not already
    # if tmin/tmax set by user without units, assume Myr
    if tmin == None:
        tmin = 0.0 * yt.units.Myr
    elif not hasattr(tmin,'value'):
        tmin = tmin * yt.units.Myr

    if tmax == None:
        tmax = 1.0E4 * yt.units.Myr
    elif not hasattr(tmax,'value'):
        tmax = tmax * yt.units.Myr
        
    ds_min  = np.argmin(np.abs((sim.times[ftype] - tmin).value))
    ds_max  = np.argmin(np.abs((sim.times[ftype] - tmax).value)) + 1

    dn = 1
    if not dt == None:
        if not hasattr(dt, 'value'): # assume Myr if not provided
            dt = dt * yt.units.Myr
        print 'if not dt None'
        dt_sim = sim.times['plt'][1:] - sim.times['plt'][0:-1]
        dt_avg = np.median(dt_sim)

        dn = dt / dt_avg
        if hasattr(dn, 'value'):
            dn = dn.value
        print dt_sim, dt_avg, dt, dn
    
    dn = int(dn)
    ds_selection = np.arange(ds_min, ds_max + dn, dn)
    print dn, ds_selection
    ds_selection = ds_selection[ds_selection <= ds_max]

    if np.size(ds_selection) > np.size(sim.times['plt']):
        ds_selection = ds_selection[0:np.size(sim.times['plt'])]
    return ds_selection

def dwarf_mass(sim, out_file, tmin=None, tmax=None, dt=None, mode='grav', T_range=[],
               neutral_temp = 2.0E4 * yt.units.Kelvin):
    """
       Calculate the gas mass of the dwarf as a function of 
       time. Method set by mode to 'grav' or 'contained' for gravitationally
       bound gas vs. total mass within dwarf radius. 'grav' assumes
       neutral, primordial gas below neutral_temp and ionized, primordial
       above. Assumes static, analytic potential for the DM halo with NO
       self gravity of gas.

    Parameters
    ----------
    sim : simulation class object
        The simulation data set(s) to operate over
    out_file : string
        Filename/path to write results to. Printed out in two colums
        time [Myr], mass [Msun]
    tmin : float, optional
        Beginning of time range to calculate over. If none, set to 0.0 Myr
    tmax : float, optional
        End of time range to calculate over. If none, set to 14 Gyr
    mode : string, optional
        Mass calulation mode to operate with. Options are 'grav' and
        'contained'.
        'grav': total mass of graviationally bound gas in the entire box
        (defined as E_thermal + E_kinetic < U_grav). Computes gravitational
        energy from analytic potentials implemented in "profiles.py". No
        support for self gravitating / time variable potentials.
        'contained': sums mass contained within initial dwarf radius
        'contained_evolve': sums mass contained within evolving radius
                     computed with dwarf.dwarf_radius function
    T_range : 2 element list or array, optional
        If supplied, performs a temperature cut, ignoring gas outside
        the supplied range. Default, no cuts

    neutral_temp : float with yt units kelvin, optional
        Only used in 'grav'. Thermal energy depends on number density,
        requiring value for mean molecular weight. Best case assumtion
        is gas below neutral_temp is neutral (mu = 1.31), above ionized
        mu = 0.6. Default 2.0E4 Kelvin.
    """

    ds_list = sim.plt_list
    if len(T_range) == 2:
        T_range = np.array(T_range)*yt.units.Kelvin

    
    sim.dwarf_mass = {}
    sim.dwarf_mass[mode] = np.ones(np.size(sim.times['plt']))*-1

    
    # do only over select time range
    ds_selection = _select_ds_list(sim, tmin, tmax, dt=dt)


    file = open(out_file, 'w')
    file.write("# t m\n")
    format = "%8.8e %8.8e\n"
    if mode == 'grav':

        # calculate dwarf gas mass as total bound gas mass of dwarf
        # 1) --- do over whole box! -- 
        # 2) get thermal energy and kinetic energy in each cell.. sum..
        # 3) get position in each cell (r) and compute the value 
        #    of the potential
        # 4) get mass in each cell
        # 5) sum cell mass where   E_th + E_kin - m*Phi(r) < 0
 

        i = 0
        for k in ds_selection:

            dsname = ds_list[k]


            _myprint("Calculating mass for file %4i of %4i"%(k,np.size(ds_selection)))
            ds = yt.load(dsname); data = ds.all_data()
            x = data['x'].convert_to_units('cm')
            y = data['y'].convert_to_units('cm')
            z = data['z'].convert_to_units('cm')  
 
            rho = data['dens']
            mass      = rho * data['dx'] * data['dy'] * data['dz']        
            mass = mass.convert_to_units('g')

            T = data['temp'].convert_to_units('K')

            # kinetic energy
            E_kin = 0.5 * mass * (data['velx']**2 + data['vely']**2 + data['velz']**2)
            E_kin = E_kin.convert_to_units('erg')

            # thermal energy is calculated in an approximate matter since
            # mu is NOT tracked / evolved in the simulation.... 
            # if T < neutral temp, assume neutral and mu = 1.31
            # if T > neutral temp, assume ionized and mu = 0.60
#            E_therm = 1.5 * mass * T 
#            E_therm = E_therm / yt.physical_constants.mass_hydrogen_cgs * yt.physical_constants.kboltz
#            E_therm[T <= neutral_temp] *= 1.0/cgs.mu_neutral
#            E_therm[T >  neutral_temp] *= 1.0/cgs.mu_ionized
#
#            E_therm = E_therm.convert_to_units('erg')
            E_therm = mass * data['eint'].convert_to_units('cm**2/s**2')

            # total energy
            E_tot = E_therm + E_kin

            # calcualte gravitational potential energy
            r = sim.dist_from_center(x,y,z)
            phi       = sim.evaluate_potential(r)
            U_grav    = -1.0 * mass * phi

            if len(T_range) == 2:
                total_mass = np.sum( mass[(E_tot<U_grav)*(T>T_range[0])*(T<T_range[1])] )
           
            else:
                total_mass = np.sum(mass[(E_tot<U_grav)])

            total_mass = total_mass.convert_to_units('Msun')
           

            sim.dwarf_mass[mode][k] = total_mass
            file.write(format%(ds.current_time.convert_to_units('Myr'),total_mass))
            i = i + 1
    
    elif mode == 'contained':
        # calculate the total dwarf mass as just the total mass contained
        # within the initial dwarf radius (with optional temperature cuts)
        i = 0
        for dsname in ds_list[ds_min:ds_max]:
            _myprint("Calculating mass for file %4i of %4i"%(i+ds_min+1,ds_max-ds_min+1))
            
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
        
def dwarf_radius(sim, outfile, tmin=None, tmax=None, 
                 mode = 'slice', density_limit=1.0E-26, nsample=12.0,
                 dr_cell = 16.0, f_bound = 0.95, T_bound_limit = 2.0E4,
                 ds_selection=None):
    """
    
        Compute the estimated dwarf galaxy radius at the supplied time
        steps. Two modes are available 'slice' and 'grav'. 

        'slice' takes 2D slices through dwarf center in each axis
        ('x','y','z') and draws nsample equal length rays in each,
        evenly spaced in angle. Radius is determined in each ray as the
        farthest cell with density >= density_limit. Final radius is 
        average over every ray in every axis (i.e. 3 * nsample rays).

        'grav' computes the gas mass in a series ofspherical annuli
        bound to the dwarf as a fraction of total gas mass in each annuli.
        Radial seperation of annuli is controlled by dr_cell (dr_cell * median
        grid cell size) to scale with resolution. Annuli are examined from 
        closest to dwarf center outward, and radius is determined by first 
        annulus with bound_mass / total_gas_mass < f_bound. The radius 
        scales with chosen f_bound as r = (f_bound)*(r_upper - r_lower) + r_lower
        where r_upper is outer edge of annulus, r_lower is lower edge. This way,
        f_bound = 1 gives r_upper as radius, f_bound = 0 gives r_lower.
        

    Parameters
    ----------
    
    """

    if not hasattr(density_limit,'convert_to_units'):
        # assume units are supplied in cgs if not provided
        density_limit = density_limit* yt.units.g / yt.units.cm**3

    if not hasattr(T_bound_limit,'convert_to_units'):
        T_bound_limit = T_bound_limit * yt.units.Kelvin

    ds_list = sim.plt_list

    if tmin == None:
        tmin = 0.0 * yt.units.Myr
    if tmax == None:
        tmax = 1.0E4 * yt.units.Myr


    sim.dwarf_radius = {}
    sim.dwarf_radius[mode] = np.ones(np.size(sim.times['plt']))*-1 * yt.units.pc


    # do only over select time range
    if ds_selection == None:
        ds_selection = _select_ds_list(sim, tmin, tmax, dt=dt)
 
    
    file = open(outfile, 'w')
    file.write("# t r\n")
    format = "%8.8e %8.8e\n"

    if mode == 'slice': # do x y and z slices of dwarf
        rmax = 2.0 * sim.radius.convert_to_units('cm') # maximum ray distance

        theta = np.arange(0.0, 2.0*np.pi, 2.0*np.pi / (1.0*nsample))

        i = 0
        for dsname in ds_selection:
            _myprint("Calculating radius for file %4i of %4i"%(i+ds_min+1,ds_max-ds_min+1))
            ds = yt.load(dsname)

            sim.dwarf_radius[mode][i + ds_min] = 0.0 * yt.units.pc # initialize
            for slice_axis in ['x','y','z']:
                i_vals = np.ones(nsample); j_vals = np.ones(nsample); k_vals=np.ones(nsample)
            
                if slice_axis == 'x':
                    j_vals = np.cos(theta)
                    k_vals = np.sin(theta)
                elif slice_axis == 'y':
                    i_vals = np.cos(theta)
                    k_vals = np.sin(theta)
                elif slice_axis == 'z':
                    i_vals = np.cos(theta)
                    j_vals = np.sin(theta)

                # make array of final positions
                final_pos = rmax * np.array([i_vals, j_vals, k_vals])
                final_pos = final_pos.transpose()
  
                # now loop over ray drawing
                for fp in final_pos:
                    ray = ds.ray(sim.center, fp)
                    ray_sort = np.argsort(ray['t'])

                    r_ray = ray['t'][ray_sort] * np.sqrt(np.sum(ray.vec**2))
                    r_ray = r_ray.convert_to_units('cm')
                    dens  = ray['dens'][ray_sort].convert_to_units('g/cm**3')
    
                    rvals = r_ray[dens >= density_limit]
                    if np.size(rvals) == 0:
                        radius = 0.0 * yt.units.pc
                    else:
                        radius = np.max(rvals)

                    radius = radius.convert_to_units('pc')
                    sim.dwarf_radius[mode][i+ds_min] += radius
                # end ray loop            
            
                
            # end slice loop
            # average over number of slices
            sim.dwarf_radius[mode][i + ds_min] *= 1.0 / (nsample*3.0)
            file.write(format%(ds.current_time.convert_to_units('Myr').value,sim.dwarf_radius[mode][i+ds_min].value))
            i = i + 1

    elif mode == 'grav': # do x y and z slices of dwarf
        # compute shell thickness 
        rmax = 2.0 * sim.radius.convert_to_units('cm') # maximum ray distance

        r_prev = sim.radius.convert_to_units('cm')

        if f_bound == 1.0:
            _myprint("f_bound CANNOT BE 1.... changing to 0.95")
            f_bound = 0.95

        i = 0
        for dsname in ds_selection:
           

            _myprint("Calculating radius for file %4i of %4i"%(i+ds_min+1,ds_max-ds_min+1))
            ds = yt.load(dsname)

            # compute shell thickness as a function of resolution
            # if resolution is too low, fix dr
            sp = ds.sphere(sim.center,sim.radius)
            dr = np.min([sim.radius.convert_to_units('cm').value * 0.01, dr_cell*np.median(sp['dx'].convert_to_units('cm').value)])
            if not hasattr(dr, 'convert_to_units'):
                dr = dr * yt.units.cm

            # calculate radii of spheres
            rmin = np.max([0.5*r_prev.convert_to_units('cm').value,np.min(sp['dx'].convert_to_units('cm').value)])
            
            sphere_radii  = np.arange(rmin*yt.units.cm, 2.0*r_prev.convert_to_units('cm'), dr)
            if not hasattr(sphere_radii, 'convert_to_units'):
                sphere_radii  = sphere_radii * yt.units.cm
            
            print dr.convert_to_units('pc')
            print sphere_radii.convert_to_units('pc')
            # initialize values of inner sphere....
            # this method calculate bound mass in spherical shells...
            #  outer sphere -inner sphere
            # precalculating and saving inner sphere prevents
            # calculating each sphere twice in the loop
            sp = ds.sphere(sim.center, sphere_radii[0])
            x = sp['x'].convert_to_units('cm')
            y = sp['y'].convert_to_units('cm')
            z = sp['z'].convert_to_units('cm')

            rho = sp['dens']
            mass      = rho * sp['dx'] * sp['dy'] * sp['dz']
            mass = mass.convert_to_units('g')

            T = sp['temp'].convert_to_units('K')

            # kinetic energy
            E_kin = 0.5 * mass * (sp['velx']**2 + sp['vely']**2 + sp['velz']**2)
            E_kin = E_kin.convert_to_units('erg')
            E_therm = mass * sp['eint'].convert_to_units('cm**2/s**2')

            # total energy
            E_tot = E_therm + E_kin

            # calcualte gravitational potential energy
            r = sim.dist_from_center(x,y,z)
            phi       = sim.evaluate_potential(r)
            U_grav    = -1.0 * mass * phi
            T_cold = (T<=T_bound_limit)
            mass = mass[T_cold]
            E_tot = E_tot[T_cold]
            U_grav = U_grav[T_cold]
            # total mass and bound mass
            prev_total_mass       = np.sum(mass)
            prev_total_bound_mass = np.sum(mass[E_tot<U_grav]) 

            j = 1
            grav_bound = True
            # for max interations or until some fraction of total mass
            # is not bound
            while (j < np.size(sphere_radii)) and (grav_bound):

                sp = ds.sphere(sim.center, sphere_radii[j])
                x = sp['x'].convert_to_units('cm')
                y = sp['y'].convert_to_units('cm')
                z = sp['z'].convert_to_units('cm')
  
                rho = sp['dens']
                mass      = rho * sp['dx'] * sp['dy'] * sp['dz']
                
                mass = mass.convert_to_units('g')

                T = sp['temp'].convert_to_units('K')
   
                # kinetic energy
                E_kin = 0.5 * mass * (sp['velx']**2 + sp['vely']**2 + sp['velz']**2)
                E_kin = E_kin.convert_to_units('erg')
                E_therm = mass * sp['eint'].convert_to_units('cm**2/s**2')

                # total energy
                E_tot = E_therm + E_kin
    
                # calcualte gravitational potential energy
                r = sim.dist_from_center(x,y,z)
                phi       = sim.evaluate_potential(r)
                U_grav    = -1.0 * mass * phi

                T_cold = (T<=T_bound_limit)
                mass = mass[T_cold]
                E_tot = E_tot[T_cold]
                U_grav = U_grav[T_cold]

                total_mass       = np.sum(mass)
                total_bound_mass = np.sum(mass[E_tot<U_grav])
 
                fractional = np.abs(total_bound_mass - prev_total_bound_mass) / prev_total_bound_mass
                #print total_bound_mass , total_mass
                # check to see what fraction is bound
                if ((total_bound_mass - prev_total_bound_mass) <\
                       f_bound * (total_mass - prev_total_mass)) or\
                       (fractional < 1.0E-6):
                    grav_bound = False
#                    print '-----',radius, j
                    # instead of taking midpoint of unbound shell, choose 
                    # dwarf radius depending on how strict f_bound is (i.e. if 
                    # f_bound = 1.0, dwarf radius is radius of inner sphere, if 
                    # f_bound = 0.5, radius is halfway between inner and outer
                    if fractional < 1.0E-6:
                        radius = sphere_radii[j-1]
                    else:
                        radius     = (f_bound)*(sphere_radii[j] - sphere_radii[j-1]) + sphere_radii[j-1]

#                    print '-----',radius,j
                prev_total_bound_mass = total_bound_mass
                prev_total_mass       = total_mass
                j = j + 1
            # end shell loop

            # save to array and dict
            sim.dwarf_radius[mode][i + ds_min] = radius.convert_to_units('pc')
            file.write(format%(ds.current_time.convert_to_units('Myr').value,radius.value))

            r_prev = radius
            i = i + 1
        # end loop over ds files
    #------------------------------------------------------------------------------------------    
    elif mode == 'grav2': 
        # compute shell thickness 
        r_prev = sim.radius.convert_to_units('cm')

        i = 0
        ds_list = sim.plt_list
        for k in ds_selection:

        # load file and create sphere on simulation center
            ds_name = ds_list[k]

            _myprint("Calculating radius for file %4i of %4i"%(i+1,len(ds_selection)))
            ds = yt.load(ds_name)
            data = ds.all_data()
            
            # get positions of all points
            x = data['x'].convert_to_units('cm')
            y = data['y'].convert_to_units('cm')
            z = data['z'].convert_to_units('cm')

            rho = data['dens']
            mass      = rho * data['dx'] * data['dy'] * data['dz']
            mass = mass.convert_to_units('g')

            T = data['temp'].convert_to_units('K')

            # kinetic energy and thermal energy
            E_kin = 0.5 * mass * (data['velx']**2 + data['vely']**2 + data['velz']**2)
            E_kin = E_kin.convert_to_units('erg')
            E_therm = mass * data['eint'].convert_to_units('cm**2/s**2')

            # total energy
            E_tot = E_therm + E_kin

            # calcualte gravitational potential energy
            r = sim.dist_from_center(x,y,z)
            r = r.convert_to_units('cm')
            phi       = sim.evaluate_potential(r)
            U_grav    = -1.0 * mass * phi
            
            # select out gas using tempeature cutoff
            T_cold = (T<=T_bound_limit)
            r       = r[T_cold]  
            mass    = mass[T_cold]
            E_tot   = E_tot[T_cold]
            U_grav  = U_grav[T_cold]
            
            # find bound cells:
            bound_cells = E_tot < U_grav
            
            # radii and masses of bound cells only
            
            R_bound = r[bound_cells]
            M_bound = mass[bound_cells]
            total_bound_mass = np.sum(M_bound)
            print np.median(R_bound.convert_to_units('pc')), np.average(R_bound.convert_to_units('pc'))
            print np.max(R_bound.convert_to_units('pc')), np.min(R_bound.convert_to_units('pc'))
            # now sort R and M by R:
            R_argsort = np.argsort(R_bound)
            R_bound = R_bound[R_argsort]
            M_bound = M_bound[R_argsort]
            
            # get array of cumulative bound mass fraction
            cum_bound_fraction = np.cumsum(M_bound) / total_bound_mass

            # now get the index where bound fraction is greater than or equal to f_bound
            index = ((cum_bound_fraction >= f_bound).tolist()).index(True)
            
            # now retrieve the dwarf radius
            dwarf_radius = R_bound[index]
            
            # save to array and dict
            sim.dwarf_radius[mode][k] = dwarf_radius.convert_to_units('pc')
            file.write(format%(ds.current_time.convert_to_units('Myr').value,dwarf_radius.value))

            i = i + 1
        # end loop over ds files

    file.close()
    return 

def profile_1D(sim, x_field, y_field, nbins=10, weight_field='cell_mass', tmin = None, tmax = None,
               dt = None, accumulation=False, ds_selection=None):
    """
    
    """
    
    ds_list = sim.plt_list
    
    # set tmin/tmax if not already
    # if tmin/tmax set by user without units, assume Myr
    if tmin == None:
        tmin = 0.0 * yt.units.Myr
    elif not hasattr(tmin,'value'):
        tmin = tmin * yt.units.Myr

    if tmax == None:
        tmax = 1.0E4 * yt.units.Myr
    elif not hasattr(tmax,'value'):
        tmax = tmax * yt.units.Myr
    
    if not hasattr(x_field, '__iter__'):
        x_field = [x_field]
    if not hasattr(y_field, '__iter__'):
        y_field = [y_field]
    
    if ds_selection == None:
        ds_selection = _select_ds_list(sim, tmin, tmax, dt=dt)

    np.ones(np.size(sim.times['plt']))*-1
    
    # need to add profiles to dict
    #  ... if sim doesn't have profiles attribute, add it
    if not hasattr(sim, 'profiles'):
        sim.profiles = {}
    if not hasattr(sim, 'profile_bins'):
        sim.profile_bins = {}
    if not hasattr(sim, 'profile_times'):
        sim.profile_times = [None]*np.size(sim.times['plt'])
    
    for x in x_field:
        if not x_field in sim.profiles.keys():
            sim.profiles[x] = {}
        if not x_field in sim.profile_bins.keys():
            sim.profile_bins[x] = [None]
        
        for y in y_field: # initialize lists to contain data for all times
            sim.profiles[x][y] = [None]*np.size(sim.times['plt'])

    # loop over the selected data files
    for k in ds_selection:

        # load file and create sphere on simulation center
        ds_name = ds_list[k]
        ds = yt.load(ds_name)
        sp = ds.sphere(sim.center,sim.radius)
        
        # calculate all of the desired profiles and save them
        for x in x_field:
            profile = yt.Profile1D(sp, x, nbins, 0.0*yt.units.pc, 
                                   sp.radius,x_log=False, weight_field=weight_field)
    
            for y in y_field:
                profile.add_fields(y)
                
                if profile.field_data.has_key(y):
                    store_profile = profile.field_data[y]
                else:
                    store_profile = profile.field_data[('gas',y)]

                if accumulation:
                    store_profile = np.cumsum(store_profile)
                
                sim.profiles[x][y][k] = store_profile
            sim.profile_times[k] = ds.current_time.convert_to_units('Myr')
        sim.profile_bins[x] = profile.x
                
                

        
def _myprint(string):
    """
        special print function just to prepend some text before every print
    """
    print "DWARF ANALYSIS : " + string
    return


#def analytic_RPS():
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

        
            
  
       
        return prof.x_bins, prof[field]
