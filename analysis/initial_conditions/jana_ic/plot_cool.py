import numpy as np
import matplotlib.pyplot as plt

lw=1.75

data = np.genfromtxt('plotcool.dat',names=True,delimiter=",")

fig,ax1 = plt.subplots()

ax1.plot(data['n'],data['T'],label="Temperature",lw=lw,color='black')
ax1.loglog()
ax1.set_xlabel(r'n (cm$^{-3}$)'); ax1.set_ylabel(r'T (K)')

ax2 = ax1.twinx()
ax2.plot(data['n'], data['nT'], label='Pressure',lw=lw,color='red')
ax2.set_ylabel(r'Pressure (cm$^{-3}$ K)')
ax2.semilogy()

plt.legend(loc='best',fancybox=True)
plt.savefig('plotcool.png')
plt.close()
