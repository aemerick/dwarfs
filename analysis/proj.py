import yt
import numpy as np
import os
from astropy import units as u
import dwarf as danal
import derived_fields as mydf
import cgs as cgs

square_box = False
imin = 84
imax = 1000
di   = 1
name = 'dwarf_fullp_hdf5_plt_cnt_'
ilist = range(imin,imax,di)

square_width = ((15.,'kpc'),(15.,'kpc'))

cmap = 'RdYlBu_r' #'viridis'
temperature_cmap = 'RdYlBu_r' #'algae'

axis = 'z'

full_box = True
#ilist = [34,41,53]
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


#yt.add_field('v_r', function=_v_r, units="cm/s", display_name = r'v$_{\rm{r}}$', take_log=False)
#yt.add_field('g_r', function=_g_r, units="cm/s**2", display_name = r'g$_{\rm{r}}$')
yt.add_field('blob_center', function=_blob_center, units="cm")

mydf.add_field('gas_column_density')
mydf.add_field('gas_number_density')
#name = 'dwarf_fullp_hdf5_chk_'
#name = 'blobtest_hdf5_chk_'
#fields_list = ['velx','g_r']
#fields_list = ['mach_number','sound_speed']
#fields_list = ['velx']
#fields_list = ['density','temperature','pressure','density','velx','vely','velx']#,'v_r','velx']
#fields_list = ['density','temperature','pressure','density']
#fields_list = ['gas_number_density']
fields_list = ['density','temperature']#,'pressure','velx','vely','velz']
weight_field = 'density'
#blob_center = np.array([8.00E21,3.0857E21,3.0857E21])
#blob_center  = np.array([3.0857E21,3.0857E21,3.0857E21])
#fields_list = ['density']


ds = yt.load(name + "%0004i"%(0))


dw = danal.dwarf(ds, 'flash.par', rm=None)
blob_center = dw.center

#left_edge  = blob_center - np.array([2.0,2.0,2.0])*cgs.kpc
left_edge = blob_center - np.array([2.0,2.0,2.0])*cgs.kpc
left_edge = (left_edge * yt.units.cm)
right_edge = blob_center + np.array([2.0,2.0,2.0])*cgs.kpc
#right_edge = blob_center + np.array([6.0,2.0,2.0])*cgs.kpc
right_edge = (right_edge * yt.units.cm)

print blob_center, left_edge, right_edge

RM = (300.0 * yt.units.pc)
add_radius = True
heating = False
if dw.params.has_key('useHeat'):
    if dw.params['useHeat'] == '.true.':
        heating = True

#ilist = range(imin,imax,di)
#ilist = range(0,1001,100)

weight_field='dens'
#axis = 'z'
project = False
velocity = False
scale_field = False # if true, let yt decide how to scale color bar


if project:
    outfolder = 'proj'
else:
    outfolder = 'slice'
outfolder = outfolder + axis

#i = imin
#for field in fields_list:
extraStr = ''
#while i <= imax:
heating = True
for i in ilist:

    extraStr = ''
    j = 0
    for field in fields_list:
        extraStr = ''
        ds = yt.load(name + "%0004i"%(i))
        print ds.field_list
        if full_box:
            left_edge = ds.domain_left_edge.convert_to_units('cm')
            right_edge = ds.domain_right_edge.convert_to_units('cm')

        reg = ds.region(blob_center*yt.units.cm,left_edge, right_edge)

        time = (ds.current_time.value * u.s).to(u.Myr).value    
        if project:
#            prsl = yt.ProjectionPlot(ds, 2, field, weight_field=weight_field,data_source=reg)
#            center = 0.5*(ds.domain_right_edge - ds.domain_left_edge).convert_to_units('cm')
            
#            width = ((2.,'kpc'),(2.,'kpc'))
#            width = ((2.0,'kpc'),(2.0,'kpc'))
#            width  = ((4.5,'kpc'),(2,'kpc'))

#            center = np.array([3,5,5])*yt.units.kpc
#            width = ((1.5,'kpc'),(1.5,'kpc'))
                   # center = np.array([4.0,4.0,4.0])*yt.units.kpc
#            center = np.array([1.234271032E22,1.234271032E22,1.1234271032E22])*yt.units.cm
           # center = np.array([5,5,5])*yt.units.kpc
           # width  = ((4.5,'kpc'),(4,'kpc'))
            if square_box:
                width = square_width
                 # center = np.array([4.0,4.0,4.0])*yt.units.kpc
                 #  center = np.array([4.0,4.0,4.0])*yt.units.kpc
                center = np.array([ ds.parameters['sim_xctr'],
                                      ds.parameters['sim_yctr'],
                                      ds.parameters['sim_zctr']])*yt.units.cm
            else:
                center = np.array([ ds.parameters['sim_xctr'],
                                        ds.parameters['sim_yctr'],
                                        ds.parameters['sim_zctr']])*yt.units.cm
                center[0] = center[0] + (1.5 * yt.units.kpc).convert_to_units('cm')


                width  = ((4.5,'kpc'),(2,'kpc'))
   

            prsl = yt.ProjectionPlot(ds, axis, field, weight_field=weight_field,data_source=reg,center=center,width=width)


#            prsl.annotate_velocity()                      
#            prsl.save('projections/')
        else:

            if (axis == 'x'):
                prsl = yt.SlicePlot(ds, axis, field, center = blob_center)
            else:
                #center = blob_center
#                width  = 
                center = 0.5*(ds.domain_right_edge - ds.domain_left_edge).convert_to_units('cm')

           #    center = np.array([1.222518940E+22, 1.265911281E+22, 1.222518940E+22])*yt.units.cm
                if square_box:
                    width = square_width
                   # center = np.array([4.0,4.0,4.0])*yt.units.kpc
                  #  center = np.array([4.0,4.0,4.0])*yt.units.kpc
                    center = np.array([ ds.parameters['sim_xctr'],
                                        ds.parameters['sim_yctr'],
                                        ds.parameters['sim_zctr']])*yt.units.cm
                else:
                    center = np.array([ ds.parameters['sim_xctr'],
                                        ds.parameters['sim_yctr'],
                                        ds.parameters['sim_zctr']])*yt.units.cm
                    center[0] = center[0] + (1.5 * yt.units.kpc).convert_to_units('cm')
                    width  = ((4.5,'kpc'),(2,'kpc'))
                 #   center = np.array([2.25,5.3,5])*yt.units.kpc
                #    width  = ((1.2,'kpc'),(0.6,'kpc'))
               

                prsl = yt.SlicePlot(ds, axis, field, center = center,width=width)
            
            if axis == 'x':
                ii,jj = 2,3
            elif axis == 'y':
                ii,jj = 1,3
            else:
                ii,jj = 1,2

            if not full_box:
                width = ((right_edge[jj] - left_edge[jj])**2 + (right_edge[ii]-left_edge[ii])**2)**0.5
                width.convert_to_units('kpc')
                prsl.set_width(width)
 #           slc.annotate_velocity()
 #           slc.annotate_velocity()
#        slc.save('slice/' + field)
       
        if velocity:
            prsl.annotate_velocity()
            extraStr = '_vel'
  
        #if add_radius:
        #    prsl.annotate_line( RM.value*np.cos(np.linspace(0.0,2.0*np.pi,200.0)) ,
        #                        RM.value*np.sin(np.linspace(0.0,2.0*np.pi,200.0)) ,
        #                        plot_args={'linewidth':5,'color':'black','linestyle':'--'})

        if j == 2:
            prsl.annotate_grids()
            extraStr = '_grid'

        if not scale_field:
            if field == 'temperature':
                if heating:
                     prsl.set_zlim(field, 5800.0,8E6)
 #                   prsl.set_zlim(field,8000.0, 3E6)   
                else:
                    prsl.set_zlim(field, 5800.0, 8E6)
            elif field == 'density':
                if heating:
                    prsl.set_zlim(field, 1.0E-28, 1.0E-23)
            #        prsl.set_zlim(field, 1.0E-28, 8.0E-25)
                else:
                    prsl.set_zlim(field,1.0E-28,1.0E-24)
            elif field == 'pressure':
                prsl.set_zlim(field,5E-15,5.E-12)
 
        if not os.path.isdir(outfolder + "/" ):
            os.makedirs(outfolder + "/")
        if not os.path.isdir(outfolder + "/" + field + extraStr):
            os.makedirs(outfolder + "/" + field + extraStr)
 

        if not field == 'g_r':
            prsl.annotate_title('t = %0.1f Myr'%(time))
         
        if field == 'velx' or field =='vely' or field == 'velz' or\
           field == 'v_r':
            data = prsl.data_source
            vmin,vmax = np.min(data[field]), np.max(data[field])
            vmin = 1.0E4
            vmax = 250E5
            prsl.set_zlim(field,vmin,vmax)
            
#            vlimit = np.max([np.abs(vmin),np.abs(vmax)])
 #           prsl.set_log(field, False)
  #          prsl.set_zlim(field, -vlimit,vlimit)
            #prsl.set_cmap(field, 'seismic')
 
        j = j + 1
        prsl.set_buff_size((4000,2000))

        if (field == 'Temperature' or field == 'temperature' or field == 'temp'):
            prsl.set_cmap(field, cmap=temperature_cmap)
        else:
            prsl.set_cmap(field, cmap=cmap)
        prsl.save(outfolder + "/" + field + extraStr)
 #           slc.save('slice/velocity/') 
#    i = i + 1

