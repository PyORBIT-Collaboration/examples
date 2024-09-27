#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=PPU_Method2_strippersNotClosed_ArrayConfig_p20_Early_Separate_Horiz
outputDirectorySub=Method2_Part1_strippersNotClosed_FloatLength
#beamLatticeFileSuffix=(RUL LUR)
#beamLatticeFileSuffix=(ULLL URRR)
beamLatticeFileSuffix=(ULRD)
#beamLatticeFileSuffix=(RUR)
#beamLatticeFilePrefix=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed_ArrayConfig_p20_Early_
beamLatticeFilePrefix=OptimizerConfigFiles/PPU_DefaultBeamLattice_strippersNotClosed_ArrayConfig_p20_Early_
optimizerConfigFile=OptimizerConfigFiles/Method2_Part1_Settings_ModifyStripper.txt


chicaneScaleDirectory=Method2_Part1_strippersNotClosed_FloatLength

echo "Method2 Part 1"
for i in ${beamLatticeFileSuffix[@]}
do
	#continue
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"
done
chicaneScaleDirectory=Method2_Part1_strippersNotClosed_FloatLength
echo "Method2 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method2_Part2_Settings_ArrayConfig_FloatLength_JustY.txt
outputDirectorySub=Method2_Part2_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	#continue
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done
echo "Method2 Part 3"
chicaneScaleDirectory=Method2_Part2_strippersNotClosed_FloatLength
optimizerConfigFile=OptimizerConfigFiles/Method2_Part3_Settings_ArrayConfig_FloatLength_JustX_1st.txt
outputDirectorySub=Method2_Part3_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done
