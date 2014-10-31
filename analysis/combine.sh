#!/bin/bash

name='blobtest_hdf5_chk_'

for i in $(seq -f "%04g" 0 120)
do

    fname1="./density/""$name$i""_Slice_z_density.png"
    fname2="./pressure/""$name$i""_Slice_z_pressure.png"
    fname3="./temperature/""$name$i""_Slice_z_temperature.png"
    fname4="./gpot/""$name$i""_Slice_z_gpot.png"

    outname="./panel/""$name$i""_panel_plot.png"

    montage $fname1 $fname2 $fname3 $fname4 -tile 2x2 -geometry 1085x520 $outname

#    j=expr "$i" % 10
#    if [ "$j" -eq 0 ]
#      then
#        echo "completed picture number $i"
#    fi

done
