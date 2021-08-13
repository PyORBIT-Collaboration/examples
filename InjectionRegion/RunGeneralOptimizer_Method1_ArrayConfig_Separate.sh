#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method1_strippersNotClosed_ArrayConfig_Separate_Late
outputDirectorySub=Method1_Part1_strippersNotClosed_FloatLength
beamLatticeFileSuffix=(LRL RLR)
#beamLatticeFileSuffix=(RUR)
beamLatticeFilePrefix=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed_ArrayConfig_Strong_Late_
optimizerConfigFile=OptimizerConfigFiles/Method1_Part1_Settings_ArrayConfig.txt


chicaneScaleDirectory=Method1_Part1_strippersNotClosed_FloatLength

echo "Method2 Part 1"
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"
done
chicaneScaleDirectory=Method1_Part1_strippersNotClosed_FloatLength
echo "Method2 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method1_Part2_Settings_ArrayConfig_FloatLength_JustY.txt
outputDirectorySub=Method1_Part2_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done