import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt

k = 1.380658E-16 
mh = 1.6733E-24
G  = 6.67259E-8


rho_center = 8.3665E-24 # 5E2 * 1E-2 / mh
T1         = 9.34E3
T2         = 7.4E5      # 
bparam     = 2.46854206E21   # 8.993E-2 * 8.896E3 pc -> cm  = 800 pc
RL         = 9.27267343E20   # 3.378E-2 * 8.896E3 pc -> cm  = 300 pc
rhoRL      = 7.72027154E-29  # 4.6138E-3 * 1.0E-2 * mh
rho2rm     = 7.72027154E-29   # 4.6138E-3 * 1E-2 * mh
rho1rm     = 6.22132940E-27  # 3.718E-1 * 1.0E-2 * mh
wTaper     = (5.0*u.pc).to(u.cm).value

cs1 = (1.4 * k * T1 / mh)**0.5
cs2 = (1.4 * k * T2 / mh)**0.5
cPhi = 4.0 * np.pi * G * rho_center * bparam**3

cp2 = rhoRL / np.exp( (cPhi/(cs2*cs2*RL)) * np.log(1.0 + RL/bparam))

print 'cs2 - cp2 - cPhi'
print cs2, cp2, cPhi

#r = np.logspace(-10, 22, 1.0E6)

#rhomid = cp2 * np.exp( (cPhi/(cs2*cs2*r)) * np.log(1.0 + r/bparam))

rhi = 6.0E21
rlo = 0.0
ic  = 0
eps = 1.0E-7
nmax = 2000
rhomid = 0.0
while (abs(rhomid - rho2rm)/rho2rm > eps) and (ic <= nmax):

    rmid = 0.5 * (rhi + rlo)
    rhomid = cp2 * np.exp((cPhi/(cs2*cs2*rmid))*np.log(1.0 + rmid/bparam))

    if (rhomid > rho2rm):
        rlo = rmid
    else:
        rhi = rmid
 
    ic = ic + 1


RM = rmid

cp1 = rho1rm * np.exp((-cPhi/(cs1*cs1*RM))* np.log(1.0 + (RM/bparam)))

rvals = np.linspace(0.0, (1.5*u.kpc).to(u.cm).value, 1.0E5)


density = np.zeros(np.size(rvals))
pressure = np.zeros(np.size(rvals))
rho1   = np.zeros(np.size(rvals))
rho2   = np.zeros(np.size(rvals))

rho1[0] = cp1*np.exp((cPhi/(cs1*cs1)) / bparam)
rho2[0] = cp2*np.exp((cPhi/(cs2*cs2)) / bparam)

rho1[1:] = cp1 * np.exp( cPhi / (cs1*cs1*rvals[1:]) *np.log(1+rvals[1:]/bparam)) 
rho2[1:] = cp2 * np.exp( cPhi / (cs2*cs2*rvals[1:]) *np.log(1+rvals[1:]/bparam))

tanh2 = 5.0 - 1.0*( np.tanh((rvals - RM)/wTaper) + 1.0)
tanh1 = 1.0 - tanh2

density = (tanh1*rho1) + (tanh2*rho2)

pressure[rvals < RM] = rho1[rvals<RM] * k * T1 / mh
pressure[rvals > RM] = rho2[rvals>RM] * k * T2 / mh

rvals = (rvals*u.cm).to(u.pc).value

plt.plot(rvals, rho1, label='rho1')
plt.plot(rvals, rho2, label='rho2')
#plt.plot(rvals, density, label='smoothed density')

#plt.plot([ (RM*u.pc).to(u.cm).value, (RM*u.pc).to(u.cm).value],plt.ylim(),color='purple')
#plt.plot([300.0,300.0],plt.ylim(),color='black')
#plt.plot(plt.xlim(),[rho2rm,rho2rm],color='black',ls='-.')
plt.legend(loc='best')
#plt.plot(plt.xlim(),[cp2,cp2], color='red')
plt.savefig('rho.png')
plt.semilogy()
plt.savefig('rho_log.png')
plt.close()

plt.plot(rvals,tanh1*rho1,label='tan rho1')
plt.plot(rvals,tanh2*rho2,label='tan rho2')
plt.plot(rvals,density,label='smoothed')
plt.legend(loc='best')
plt.savefig('rho_smooth.png')

#plt.savefig('rho.png')


print 'rmatch = %5.4e'%((RM*u.cm).to(u.pc).value)

