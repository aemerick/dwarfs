# Note, uses yt to read Halogen's tipsy dataset output

import yt
import numpy as np
import cgs   as cgs
import sys

# mutliplicative factors to convert from Halogen's
# code units to cgs
unit_conversion = {'MU': 2.222962E5 * cgs.Msun,
                   'VU': cgs.kpc / (1.0E9 * cgs.yr),
                   'TU': 1.0E9 * cgs.yr,
                   'LU': cgs.kpc}
                  

def convert_to_ascii(datafile, units = 'sim', fmt=None, outfile = None):

    ds   = yt.load(datafile)
    data = ds.all_data()

    m    = data[('DarkMatter', 'Mass')].value
    pos  = data[('DarkMatter', 'Coordinates')].value
    vel  = data[('DarkMatter', 'Velocities')].value

    if units == 'cgs': # convert from code units to cgs
        pos = pos * unit_conversion['LU']
        vel = vel * unit_conversion['VU']
        m   = m   * unit_conversion['MU']

    if units == 'sim':
        pos = pos * unit_conversion['LU'] / cgs.pc
        vel = vel * unit_conversion['VU'] / 1.0E5
        m   = m * unit_conversion['MU'] / cgs.Msun


    N = np.size(m)

    if fmt is None:
        fmt = "%8.8E "


    if outfile is None:
        outfile = datafile.replace('.tipsy.std','.ascii.dat')
        outfile = outfile.replace('.tipsy.dpp.std','.ascii.dat')

    fmt_string = fmt + fmt + fmt + fmt + fmt + fmt + fmt + "\n"  

    f = open(outfile, 'w')   
    for i in np.arange(N):
        f.write(fmt_string%(m[i], pos[i,0], pos[i,1], pos[i,2],
                                 vel[i,0], vel[i,1], vel[i,2]))

    f.close()
    return

            

if __name__=='__main__':
    convert_to_ascii(sys.argv[1],outfile = sys.argv[2])    
