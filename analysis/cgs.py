# some useful CGS units we always need ...

print "using cgs py"

mu_neutral = 1.31
mu_ionized = 0.60

mu   = 0.6
mp   = 1.67262E-24
#mp   = 1.67262178e-24 # mass of proton, grams
kpc  = 3.08567758e21  # kiloparsec, cm
pc   = kpc*1e-3
Mpc  = kpc*1e3        # megaparsec, cm
Msun = 1.9891e33      # mass of sun, grams
G    = 6.67428e-8       # Gravity, cgs
kb   = 1.38065E-16
#kb   = 1.380658e-16		# Boltzmann, cgs
yr   = 3.15569e7      # year in seconds
Myr  = yr*1.0e6       # Megayear ""
gyr  = yr*1.0e9       # gigayear ""
km   = 1.0e5					# kilometer in cm
h100 = 3.24077929e-18 # hertz
c    = 29979245800.0  # speed of light [ cm/s ]
eV   = 1.60217657e-12 # 1 eV in erg

arcsec = 0.000277777777778 # 1 arcsec to degrees
mas    = arcsec / 1000.0

def P(rho,T):
	return rho*kb*T/(mu*mp)
