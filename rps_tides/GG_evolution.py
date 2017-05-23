import numpy as np
import deepdish as dd
import yt
import cgs

def GG_evolution(dslist = [], data = None, outname = 'GG_rstrip.dat'):


    ds0 = yt.load(dslist[0])
    M = ds0.parameters['DiskGravityStellarDiskMass'] * cgs.Msun
    a = ds0.parameters['DiskGravityStellarDiskScaleHeightR'] * cgs.Mpc
    b = ds0.parameters['DiskGravityStellarDiskScaleHeightz'] * cgs.Mpc
    

    star_profile = lambda x : prof.stellar_surface_density(x, b = b, M =M, a =a)

    Pram = np.logspace(-15,-10, 100) # ram pressure in cgs

    Rstrip = np.zeros( (np.size(ds_list), np.size(Pram)) )

    f = open(outname, 'w')
    # header
    f.write("# t")
    for P in Pram:
        f.write(" %3.3f"%(np.log10(Pram)))
    f.write("\n")



    i = 0
    for d in dslist:
        ds = yt.load(d)

        tnow  = ds.current_time.convert_to_units('Myr').value
        index = np.argmin( np.abs(data['t'] - tnow))

        x     = (data['r'][1:] + data['r'][:-1])*0.5
        gas   = data['surface_density'][index]

        gas_profile = profile.gas_profile_function(r, sigma, unit = 'cgs')

        j = 0
        f.write("%5.5f"%(tnow))
        for P in Pram:
            Rstrip[i][j] = stripping_radius(star_profile, gas_profile, Pram = Pram,
                                      lu = cgs.kpc)

            f.write(" %3.3E"%(Rstrip[i][j]))
            j = j + 1

        f.write("\n");

        i = i + 1

    f.close()

    return



if __name__ == "__main__":

    dslist = np.sort(glob.glob('DD????/DD????'))
    data   = dd.io.load('gas_profile_evolution.h5')

    



