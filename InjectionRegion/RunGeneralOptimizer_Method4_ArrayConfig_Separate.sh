#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method4_strippersNotClosed_ArrayConfig_PPU
#outputDirectory=Method4_strippersNotClosed_ArrayConfig
outputDirectorySub=Method4_Part1_strippersNotClosed_FloatLength
#beamLatticeFileSuffix=(RUL LUR)
#beamLatticeFileSuffix=(ULRD)
#beamLatticeFileSuffix=(p20_ULRD p30_ULRD p40_ULRD p50_ULRD p60_ULRD p70_ULRD p80_ULRD p90_ULRD)
beamLatticeFileSuffix=(p80_ULRD)
#beamLatticeFileSuffix=(p80_ULRU)
#beamLatticeFileSuffix=(RUR)
beamLatticeFilePrefix=OptimizerConfigFiles/PPU_DefaultBeamLattice_strippersNotClosed_ArrayConfig_strong_Early_first_
#beamLatticeFilePrefix=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed_ArrayConfig_strong_Early_first_
optimizerConfigFile=OptimizerConfigFiles/Method4_Part1_Settings_ArrayConfig_FloatLength.txt


chicaneScaleDirectory=Method4_Part1_strippersNotClosed_FloatLength

echo "Method4 Part 1"
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"
done
chicaneScaleDirectory=Method4_Part1_strippersNotClosed_FloatLength
echo "Method4 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method4_Part2_Settings_ArrayConfig_FloatLength.txt
outputDirectorySub=Method4_Part2_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done
echo "Method4 Part 3"
beamLatticeFilePrefix=OptimizerConfigFiles/PPU_DefaultBeamLattice_strippersNotClosed_ArrayConfig_strong_Early_second_
#beamLatticeFilePrefix=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed_ArrayConfig_strong_Early_second_
chicaneScaleDirectory=Method4_Part2_strippersNotClosed_FloatLength
optimizerConfigFile=OptimizerConfigFiles/Method4_Part3_Settings_ArrayConfig_FloatLength.txt
outputDirectorySub=Method4_Part3_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done
echo "Method4 Part 4"
chicaneScaleDirectory=Method4_Part3_strippersNotClosed_FloatLength
optimizerConfigFile=OptimizerConfigFiles/Method4_Part4_Settings_ArrayConfig_FloatLength.txt
outputDirectorySub=Method4_Part4_strippersNotClosed_FloatLength
for i in ${beamLatticeFileSuffix[@]}
do
	pyORBIT OptimizerGeneral_ArrayConfig.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFilePrefix""$i.txt" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done