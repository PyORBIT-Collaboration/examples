#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method2_strippersNotClosed_ArrayConfig
outputDirectorySub=Method2_Part1_strippersNotClosed_FloatLength
magneticFieldDirectory=MagneticFieldFiles   
magneticFieldFilePrefix=magneticField
beamLatticeFileSuffix=(RUR LUL)
#beamLatticeFileSuffix=(RUR)
beamLatticeFilePrefix=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed_ArrayConfig_
optimizerConfigFile=OptimizerConfigFiles/Method2_Part1_Settings_ModifyStripper.txt


chicaneScaleDirectory=Method2_Part1_strippersNotClosed_FloatLength

echo "Method2 Part 1"
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"
done

echo "Method2 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method2_Part2_Settings_ArrayConfig_FloatLength.txt
outputDirectorySub=Method2_Part2_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done