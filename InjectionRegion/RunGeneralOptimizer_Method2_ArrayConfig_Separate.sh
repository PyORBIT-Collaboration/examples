#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method2_strippersNotClosed_ArrayConfig_Weak_Early_Separate
outputDirectorySub=Method2_Part1_strippersNotClosed_FloatLength
beamLatticeFileSuffix=(RUL LUR)
#beamLatticeFileSuffix=(RUR)
beamLatticeFilePrefix=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed_ArrayConfig_Weak_Early_
optimizerConfigFile=OptimizerConfigFiles/Method2_Part1_Settings_ModifyStripper.txt


chicaneScaleDirectory=Method2_Part1_strippersNotClosed_FloatLength

echo "Method2 Part 1"
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"
done
chicaneScaleDirectory=Method2_Part1_strippersNotClosed_FloatLength
echo "Method2 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method2_Part2_Settings_ArrayConfig_FloatLength_JustY.txt
outputDirectorySub=Method2_Part2_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done
chicaneScaleDirectory=Method2_Part2_strippersNotClosed_FloatLength
optimizerConfigFile=OptimizerConfigFiles/Method2_Part3_Settings_ArrayConfig_FloatLength_JustX.txt
outputDirectorySub=Method2_Part3_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done
echo "Method2 Part 3"