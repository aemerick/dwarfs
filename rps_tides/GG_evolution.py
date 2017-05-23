import numpy as np
import deepdish as dd
import yt
import cgs
import glob

import matplotlib.pyplot as plt

# internal
import profiles
import GG

def GG_evolution(dslist = [], data = None, outname = 'GG_rstrip.dat'):


    ds0 = yt.load(dslist[0])
    M = ds0.parameters['DiskGravityStellarDiskMass'] * cgs.Msun
    a = ds0.parameters['DiskGravityStellarDiskScaleHeightR'] * cgs.Mpc
    b = ds0.parameters['DiskGravityStellarDiskScaleHeightz'] * cgs.Mpc

    star_profile = lambda x : profiles.stellar_surface_density(x, b = b, M =M, a =a)

    Pram = np.logspace(-15,-10, 100) # ram pressure in cgs

    Rstrip = np.zeros( (np.size(dslist), np.size(Pram)) )
    t      = np.zeros(np.size(dslist))

    f = open(outname, 'w')
    # header
    f.write("# t")
    for P in Pram:
        f.write(" %3.3f"%(np.log10(P)))
    f.write("\n")

    i = 0
    for d in dslist:
        ds = yt.load(d)

        tnow  = ds.current_time.convert_to_units('Myr').value
        t[i]  = round(tnow)
        index = np.argmin( np.abs(data['t'] - tnow))

        x     = (data['r'][1:] + data['r'][:-1])*0.5
        gas   = data['surface_density'][index]

        gas_profile = profiles.gas_profile_function(x * cgs.kpc, gas, unit = 'cgs')

        j = 0
        f.write("%i"%(round(tnow)))
        for P in Pram:
            Rstrip[i][j] = GG.stripping_radius(star_profile, gas_profile, Pram = P,
                                      lu = cgs.kpc)

            f.write(" %3.3E"%(Rstrip[i][j]))
            j = j + 1

        f.write("\n");

        i = i + 1

    f.close()

    data = {}

    data['t']    = t
    data['Pram'] = Pram

    data['Rstrip'] = Rstrip
    data['stats']  = {}

    data['stats']['avg'] = np.average(Rstrip, 0)
    data['stats']['min'] = np.min(Rstrip, 0)
    data['stats']['max'] = np.max(Rstrip, 0)
    data['stats']['std'] = np.std(Rstrip, 0)

    dd.io.save('GG_evolution.h5', data)

    return


def plot_evolution(filename):

    data = dd.io.load(filename)

    fig, ax = plt.subplots(1)

    i = 0
    pram = np.log10(data['Pram'])

    ax.plot(pram, data['stats']['avg'], label = "Average", lw = 3, color = 'black', ls = '--')

    ax.fill_between(pram, data['stats']['min'], data['stats']['max'], facecolor = 'black', interpolate=True, alpha = 0.25)
    ax.fill_between(pram, data['stats']['avg'] - data['stats']['std'],
                          data['stats']['avg'] + data['stats']['std'], facecolor = 'blue', interpolate=True, alpha = 0.75)

    i = i + 1

    ax.set_xlabel(r'log [P$_{\rm ram}$ (cgs)]')
    ax.set_ylabel(r'R$_{\rm strip}$ (kpc)')
    fig.set_size_inches(8,8)
    plt.tight_layout()
    fig.savefig('GG_evolution.png')
    plt.close()


    return

if __name__ == "__main__":

    dslist = np.sort(glob.glob('DD????/DD????'))
    data   = dd.io.load('gas_profile_evolution.h5')

    GG_evolution(dslist, data)

    plot_evolution('GG_evolution.h5')
