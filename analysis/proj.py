import yt
import numpy as np
import os

name = 'blobtest_hdf5_chk_'
fields_list = ['density','temperature','pressure','velx']
#fields_list = ['pressure']
#field = 'pressure'
#fields_list = ['density']
imin, imax = 0,500
outfolder = 'slicez'
project = False
velocity = False
scale_field = True# if true, let yt decide how to scale color bar
i = imin
#for field in fields_list:

while i <= imax:
    for field in fields_list:

        ds = yt.load(name + "%0004i"%(i))
    
        if project:
            prsl = yt.ProjectionPlot(ds, 2, field, weight_field='density')
#            prsl.annotate_velocity()          
            
#            prsl.save('projections/')
        else:
            prsl = yt.SlicePlot(ds, 'z', field)
 #           slc.annotate_velocity()
 #           slc.annotate_velocity()
#        slc.save('slice/' + field)

        if velocity:
            prsl.annotate_velocity()

        if not scale_field:
            if field == 'temperature':
                prsl.set_zlim(field,1000, 1E7)
            elif field == 'density':
                prsl.set_zlim(field,1.0E-28,1.0E-24)
            elif field == 'pressure':
                prsl.set_zlim(field,1E-14,1E-12)
 
        if not os.path.isdir(outfolder + "/" ):
            os.makedirs(outfolder + "/")
        if not os.path.isdir(outfolder + "/" + field):
            os.makedirs(outfolder + "/" + field)
 
        prsl.save(outfolder + "/" + field)
 #           slc.save('slice/velocity/') 
    i = i + 1

