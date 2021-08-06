#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method3_strippersNotClosed
outputDirectorySub=Method3_Part1_strippersNotClosed
magneticFieldDirectory=MagneticFieldFiles   
magneticFieldFilePrefix=magneticField
magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight)
optimizerConfigFile=OptimizerConfigFiles/Method3_Part1_Settings.txt
#beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice.txt
beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed.txt
chicaneScaleDirectory=Method3_Part1_strippersNotClosed

echo "Method3 Part 1"
for i in ${magneticFieldFiles[@]}
do
	pyORBIT OptimizerGeneral.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$i.txt" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFile" --chicaneScaleDirectory "$chicaneScaleDirectory"
done

echo "Method3 Part 2"
optimizerConfigFile=OptimizerConfigFiles/Method3_Part2_Settings.txt
outputDirectorySub=Method3_Part2_strippersNotClosed
for i in ${magneticFieldFiles[@]}
do
	pyORBIT OptimizerGeneral.py --outputDirectory "$outputDirectory/$outputDirectorySub"_"$i" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$i.txt" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFile" --chicaneScaleDirectory "$outputDirectory/$chicaneScaleDirectory"_"$i"
done