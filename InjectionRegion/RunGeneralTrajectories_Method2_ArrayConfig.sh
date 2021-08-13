#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method2_strippersNotClosed_ArrayConfig_Weak_Separate
#outputDirectoryPrefix=Method1
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
magneticFieldFiles=(UpUp DownDown UpDown DownUp)
#magneticFieldFiles=(LeftUp RightUp LeftDown RightDown)
#magneticFieldFiles=(LeftUp)
beamLatticeFileDirectory=BeamLatticeFiles
beamLatticeFileSuffix=(RUR LUL)
chicaneScaleDirectory=Method2_Part3_strippersNotClosed_FloatLength
suffixInjection=InjectBeam_Method2_ArrayConfig_FloatLength_Weak
suffixWaste=WasteBeam_Method2_ArrayConfig_FloatLength_Weak
suffixClosed=ClosedBeam_Method2_strippersNotClosed_ArrayConfig

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