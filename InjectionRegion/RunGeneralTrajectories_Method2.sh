#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method2_strippersNotClosed
#outputDirectoryPrefix=Method1
magneticFieldDirectory=MagneticFieldFiles   
magneticFieldFilePrefix=magneticField
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
magneticFieldFiles=(LeftUp RightUp)
beamLatticeFileDirectory=BeamLatticeFiles
chicaneScaleDirectory=Method2_Part2_strippersNotClosed
suffixInjection=InjectBeam_Method2
suffixWaste=WasteBeam_Method2
suffixClosed=ClosedBeam_Method2_strippersNotClosed

suffix=($suffixInjection $suffixWaste $suffixClosed)

#echo "$outputDirectoryPrefix"_"$suffixInjection"
echo "Method2"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectoryPrefix"_"$suffixInjection" --magneticFieldFile "$magneticFieldFile" --beamLatticeFile "$beamLatticeFileDirectory"/"$suffixInjection" --chicaneScaleDirectory "$chicaneScaleDirectory"
for j in ${magneticFieldFiles[@]}
do
	for i in ${suffix[@]}
	do
		pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectory/$i"_"$j" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$j.txt" --beamLatticeFile "$beamLatticeFileDirectory"/"$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$j"
	done
done