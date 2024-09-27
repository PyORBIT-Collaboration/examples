#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method1_strippersNotClosed_ArrayConfig_Separate_Late
#outputDirectoryPrefix=Method1
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
magneticFieldFiles=(UpUp DownDown UpDown DownUp)
#magneticFieldFiles=(LeftUp RightUp LeftDown RightDown)
#magneticFieldFiles=(LeftUp)
beamLatticeFileDirectory=BeamLatticeFiles
beamLatticeFileSuffix=(LRL RLR)
chicaneScaleDirectory=Method1_Part2_strippersNotClosed_FloatLength
suffixInjection=InjectBeam_Method1_ArrayConfig_FloatLength_Strong_Late
suffixWaste=WasteBeam_Method1_ArrayConfig_FloatLength_Strong_Late
suffixClosed=ClosedBeam_Method1_strippersNotClosed_ArrayConfig

suffix=($suffixInjection $suffixWaste $suffixClosed)

#echo "$outputDirectoryPrefix"_"$suffixInjection"
echo "Method2"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectoryPrefix"_"$suffixInjection" --magneticFieldFile "$magneticFieldFile" --beamLatticeFile "$beamLatticeFileDirectory"/"$suffixInjection" --chicaneScaleDirectory "$chicaneScaleDirectory"
for j in ${beamLatticeFileSuffix[@]}
do
	for i in ${suffix[@]}
	do
	pyORBIT InjectionRegion_NoSpaceChargeStripNode_General_ArrayConfig.py --outputDirectory "$outputDirectory/$i"_"$j" --beamLatticeFile "$beamLatticeFileDirectory"/"$i"_""$j".txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$j"
	done
done