# Note, uses yt to read Halogen's tipsy dataset output

import yt
import numpy as np
import cgs   as cgs

# mutliplicative factors to convert from Halogen's
# code units to cgs
unit_conversion = {'MU': 2.222962E5 * cgs.Msun,
                   'VU': cgs.kpc / (1.0E9 * cgs.yr),
                   'TU': 1.0E9 * cgs.yr,
                   'LU': cgs.kpc}
                  

def convert_to_ascii(datafile, units = None, fmt=None, outfile = None):

    ds   = yt.load(datafile)
    data = ds.all_data()

    m    = data[('DarkMatter', 'Mass')].value
    pos  = data[('DarkMatter', 'Coordinates')].value
    vel  = data[('DarkMatter', 'Velocities')].value

    if units == 'cgs': # convert from code units to cgs
        pos = pos * unit_conversion['LU']
        vel = vel * unit_conversion['VU']
        m   = m   * unit_conversion['MU']


    N = np.size(m)

    if fmt == None:
        fmt = "%8.8E "

    if outfile == None:
        outfile = datafile.replace('.tipsy.std','.ascii.dat')


    fmt_string = fmt + fmt + fmt + fmt + fmt + fmt + fmt + "\n"  

    f = open(outfile, 'w')   
    for i in np.arange(N):
        f.write(fmt_string%(m[i], pos[i,0], pos[i,1], pos[i,2],
                                 vel[i,0], vel[i,1], vel[i,2]))

    f.close()
    return

            

