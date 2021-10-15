#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=NoMethod_SecondDoneRight
#outputDirectory=Method4_strippersNotClosed_ArrayConfig
#outputDirectoryPrefix=Method1
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
#magneticFieldFiles=(LeftUp RightUp LeftDown RightDown)
#magneticFieldFiles=(LeftUp)
configFileDirectory=ConfigFilesCombineWeightBunches
configFileName=(101_p5_011_p5.txt 101_1p0.txt)
#configFileName=(101_p5_011_p5_LR.txt 101_1p0_LR.txt)
#beamLatticeFileSuffix=(RUR LUL)
#beamLatticeFileSuffix=(p20_ULRD p30_ULRD p40_ULRD p50_ULRD p60_ULRD p70_ULRD p80_ULRD)
inputDirectorySuffix=(p80_ULRD)
#inputDirectorySuffix=(p80_LR)
inputDirectoryPrefix=$outputDirectory/PPU_InjectBeam_NoMethod_ArrayConfig__Second

#suffix=($suffixInjection $suffixWaste $suffixClosed)  
suffix=(200 020 101 110 011)

#echo "$outputDirectoryPrefix"_"$suffixInjection"
echo "Method2"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectoryPrefix"_"$suffixInjection" --magneticFieldFile "$magneticFieldFile" --beamLatticeFile "$beamLatticeFileDirectory"/"$suffixInjection" --chicaneScaleDirectory "$chicaneScaleDirectory"

for i in ${configFileName[@]}
do
pyORBIT CombineWeightedBunches_CalcEmit.py --outputDirectory "$outputDirectory" --configFile "$configFileDirectory"/"$i" --inputDirectoryPrefix "$inputDirectoryPrefix" --inputDirectorySuffix "$inputDirectorySuffix"
done
