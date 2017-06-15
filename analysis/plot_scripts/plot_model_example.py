import matplotlib.pyplot as plt
import numpy as np
import cgs as cgs
from analytic_model import dwarf_model as dw_model
from initial_conditions import ic_list as icl ;
from matplotlib import rc

line_width = 2.5

#
simulation_dir = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'
model_name = "LT_n075_v2_nh3"
sim_mass   = model_name + '_cont_mass_1-0.dat'

data = np.genfromtxt(simulation_dir + model_name + '/' + sim_mass,names=True)
t, m = data['t'], data['m']
plt.plot(t,m/m[0])


sim_mass   = 'grav_bound_mass.dat'


#
initial_conditions = icl.ic_object_dict[model_name]
anal_model = dw_model.analytical_dwarf(model_name, initial_conditions.ic)
anal_model.load_simulation_mass(simulation_dir +model_name +'/' + sim_mass)
anal_model.setup_orbit(0.0, initial_conditions.ic['n_halo'], initial_conditions.ic['v_halo'])


tmax = anal_model.simulation_data['time'][anal_model.simulation_data['mass'] == 0.0][5]

# desired timestep size and start / end times... in Myr
dt = 0.5
tmin = 0.0
npoints = np.ceil((tmax - tmin)/dt)

t = np.linspace(tmin, tmax, npoints)*cgs.Myr

shock_kwargs = {'RPS':{'alpha': 0.25,'beta': 0.5152,'method':'shock'},'KH':{'beta': 1.0}}

#shock_kwargs = {'RPS':{'alpha': 1.475,'beta': 0.800,'method':'shock'},'KH':{'beta': 6.0}}
#shock_kwargs = {'RPS':{'alpha': 3.6591,'beta': 1.0909,'method':'shock'},'KH':{'beta': 3.8788}}

# do the full plot
M , R = anal_model.evolve(t, ['RPS','KH'], RPS_KH_exclusive = True, physics_kwargs=shock_kwargs)
plt.plot(anal_model.simulation_data['time'], anal_model.simulation_data['mass']/anal_model.simulation_data['mass'][0], label='Simulation', color='black',ls='-',lw=line_width)
plt.plot(t/cgs.Myr, M/M[0], label = 'Analytic - RPS + KH', color='red',ls='-',lw=line_width)

# just RPS
M , R = anal_model.evolve(t, ['RPS'], RPS_KH_exclusive = True, physics_kwargs=shock_kwargs)
plt.plot(t/cgs.Myr, M/M[0], label = 'Analytic - RPS', color='red',ls=':',lw=line_width)

# just KH
M , R = anal_model.evolve(t, ['KH'], physics_kwargs=shock_kwargs)
plt.plot(t/cgs.Myr, M/M[0], label = 'Analytic - KH', color='red',ls='--',lw=line_width)

plt.xlim(0,tmax)
plt.ylim(0,1)

plt.legend(loc='best',fancybox=True)

plt.savefig('LT_example_adiabatic_fit.png')

plt.close()
