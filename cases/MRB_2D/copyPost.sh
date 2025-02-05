#!/bin/bash
if [[ -d "./postProcessing" ]]
then
    echo "postProcessing exists, copying to postProcessing/Test"
else
    echo "postProcessing not found, exiting!"
    exit 0
fi
rm -rf postProcessing/Test
mkdir postProcessing/Test
#for f in postProcessing/new_plane/surface/*
for f in postProcessing/*/*
do
    tName=${f##*/} # strip time
    tName=${tName/./,} # change dot to comma
    
    tDec=$(printf "_%.8f" $tName) #make decimal number
    for k in $f/*.vtk
    do
        y=${k%.vtk}
        fName=${y##*/}
      # echo "$fName$tName.vtk"
        cpn="$fName$tDec.vtk" #new name
        cp $k postProcessing/Test/$cpn #copy
    done
done
