#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method2_strippersNotClosed_ModifyStripper
outputDirectorySub=Method2_Part1_strippersNotClosed_FloatLength
magneticFieldDirectory=MagneticFieldFiles   
magneticFieldFilePrefix=magneticField
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
#magneticFieldFiles=(LeftUp RightUp LeftDown RightDown)
magneticFieldFiles=(UpUp DownDown UpDown DownUp)
#magneticFieldFiles=(LeftUp)
optimizerConfigFile=OptimizerConfigFiles/Method2_Part1_Settings_ModifyStripper.txt
#beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice.txt
beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed.txt
#beamLatticeFile=OptimizerConfigFiles/BeamLattice_DIfferentClosed.txt
chicaneScaleDirectory=Method2_Part1_strippersNotClosed_FloatLength

echo "Method2 Part 1"
for i in ${magneticFieldFiles[@]}
do
	pyORBIT OptimizerGeneral_ModifyStripper.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$i.txt" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFile" --chicaneScaleDirectory "$chicaneScaleDirectory"
done

echo "Method2 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method2_Part2_Settings_ModifyStripper_FloatLength.txt
outputDirectorySub=Method2_Part2_strippersNotClosed_FloatLength
for i in ${magneticFieldFiles[@]}
do
	pyORBIT OptimizerGeneral_ModifyStripper.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$i.txt" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFile" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done