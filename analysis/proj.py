import yt
import numpy as np
import os
from astropy import units as u

def _v_r(field, data):
    

    center = data['blob_center']#.ds.domain_center.value * yt.units.cm
    
    vx = data['velx'].value * yt.units.cm / yt.units.s
    vy = data['vely'].value * yt.units.cm / yt.units.s
    vz = data['velz'].value * yt.units.cm / yt.units.s
    
    x = data['x'] - center[0]
    y = data['y'] - center[1]
    z = data['z'] - center[2]
    
    top = vx * x + vy * y + vz * z
    bottom = (x**2 + y**2 + z**2)**0.5
    #
    return -1.0 * (top/bottom)
#    return np.abs(top/bottom)


def _g_r(field, data):
    

    center = data['blob_center']#.ds.domain_center.value * yt.units.cm
    
    accx = data['accx'].value * yt.units.cm / yt.units.s**2
    accy = data['accy'].value * yt.units.cm / yt.units.s**2
    accz = data['accz'].value * yt.units.cm / yt.units.s**2
    
    x = data['x'] - center[0]
    y = data['y'] - center[1]
    z = data['z'] - center[2]
    
    top = accx * x + accy * y + accz * z
    bottom = (x**2 + y**2 + z**2)**0.5
    #
    return top/bottom
#    return np.abs(top/bottom)
    
def _blob_center(field, data):
    return   np.array([8.00E21,3.0857E21,3.0857E21]) * yt.units.cm


yt.add_field('v_r', function=_v_r, units="cm/s", display_name = r'v$_{\rm{r}}$', take_log=False)
yt.add_field('g_r', function=_g_r, units="cm/s**2", display_name = r'g$_{\rm{r}}$')
yt.add_field('blob_center', function=_blob_center, units="cm")

name = 'blobtest_hdf5_chk_'
fields_list = ['v_r']

fields_list = ['density','temperature','pressure','density']
blob_center = np.array([8.00E21,3.0857E21,3.0857E21])

ilist = range(0,400,1)
#ilist = range(0,1001,100)


axis = 'z'
outfolder = 'slicez'
project = False
velocity = False
scale_field = False# if true, let yt decide how to scale color bar
#i = imin
#for field in fields_list:
extraStr = ''
#while i <= imax:
for i in ilist:
    extraStr = ''
    j = 0
    for field in fields_list:
        extraStr = ''
        ds = yt.load(name + "%0004i"%(i))
#	data = ds.all_data()

        time = (ds.current_time.value * u.s).to(u.Myr).value    
        if project:
            prsl = yt.ProjectionPlot(ds, 2, field, weight_field='density')
#            prsl.annotate_velocity()          
            
#            prsl.save('projections/')
        else:

            if (axis == 'x'):
                prsl = yt.SlicePlot(ds, axis, field, center = blob_center)
            else:
                prsl = yt.SlicePlot(ds, axis, field)
 #           slc.annotate_velocity()
 #           slc.annotate_velocity()
#        slc.save('slice/' + field)
       
        if velocity or j == 3:
            prsl.annotate_velocity()
            extraStr = '_vel'

        if not scale_field:
            if field == 'temperature':
                prsl.set_zlim(field,3E4, 1E6)
            elif field == 'density':
                prsl.set_zlim(field,1.0E-28,1.0E-25)
            elif field == 'pressure':
                prsl.set_zlim(field,5E-15,5E-12)
 
        if not os.path.isdir(outfolder + "/" ):
            os.makedirs(outfolder + "/")
        if not os.path.isdir(outfolder + "/" + field + extraStr):
            os.makedirs(outfolder + "/" + field + extraStr)
 

        if not field == 'g_r':
            prsl.annotate_title('t = %0.0f Myr'%(time))
         
        if field == 'velx' or field =='vely' or field == 'velz' or\
           field == 'v_r':
            data = prsl.data_source
            vmin,vmax = np.min(data[field]), np.max(data[field])

            vlimit = np.max([np.abs(vmin),np.abs(vmax)])
            prsl.set_log(field, False)
            prsl.set_zlim(field, -vlimit,vlimit)
            prsl.set_cmap(field, 'seismic')
 
        j = j + 1
        prsl.save(outfolder + "/" + field + extraStr)
 #           slc.save('slice/velocity/') 
#    i = i + 1

