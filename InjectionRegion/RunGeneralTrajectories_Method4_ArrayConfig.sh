#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method4_strippersNotClosed_ArrayConfig_PPU_SecondDoneRight
#outputDirectory=Method4_strippersNotClosed_ArrayConfig
#outputDirectoryPrefix=Method1
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
magneticFieldFiles=(UpUp DownDown UpDown DownUp)
#magneticFieldFiles=(LeftUp RightUp LeftDown RightDown)
#magneticFieldFiles=(LeftUp)
beamLatticeFileDirectory=BeamLatticeFiles
#beamLatticeFileSuffix=(RUR LUL)
#beamLatticeFileSuffix=(p20_ULRD p30_ULRD p40_ULRD p50_ULRD p60_ULRD p70_ULRD p80_ULRD)
beamLatticeFileSuffix=(p80_ULRD)
chicaneScaleDirectory=Method4_Part4_strippersNotClosed_FloatLength
#suffixInjection=PPU_InjectBeam_Method4_ArrayConfig_FloatLength
suffixInjection=PPU_InjectBeam_Method4_ArrayConfig_FloatLength_Second000
suffixWaste=PPU_WasteBeam_Method4_ArrayConfig_FloatLength
#suffixWaste=WasteBeam_Method4_ArrayConfig_FloatLength
suffixClosed=PPU_ClosedBeam_Method4_strippersNotClosed_ArrayConfig
#suffixClosed=ClosedBeam_Method4_strippersNotClosed_ArrayConfig

#suffix=($suffixInjection $suffixWaste $suffixClosed)  
suffix=($suffixInjection)

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