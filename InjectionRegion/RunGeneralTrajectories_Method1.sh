#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method1
#outputDirectoryPrefix=Method1
magneticFieldDirectory=MagneticFieldFiles   
magneticFieldFilePrefix=magneticField
magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
beamLatticeFileDirectory=BeamLatticeFiles
chicaneScaleDirectory=Method1
suffixInjection=InjectBeam_Method1
suffixWaste=WasteBeam_Method1
suffixClosed=ClosedBeam_Method1

suffix=($suffixInjection $suffixWaste $suffixClosed)

#echo "$outputDirectoryPrefix"_"$suffixInjection"
echo "Method1"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectoryPrefix"_"$suffixInjection" --magneticFieldFile "$magneticFieldFile" --beamLatticeFile "$beamLatticeFileDirectory"/"$suffixInjection" --chicaneScaleDirectory "$chicaneScaleDirectory"
for j in ${magneticFieldFiles[@]}
do
	for i in ${suffix[@]}
	do
		pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectory/$i"_"$j" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$j.txt" --beamLatticeFile "$beamLatticeFileDirectory"/"$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"_"$j"
	done
done