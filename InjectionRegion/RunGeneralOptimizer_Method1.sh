#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method1_strippersNotClosed
magneticFieldDirectory=MagneticFieldFiles   
magneticFieldFilePrefix=magneticField
#magneticFieldFiles=(UpUp DownDown UpDown DownUp LeftLeft LeftRight RightLeft RightRight)
#magneticFieldFiles=(LeftLeft LeftRight RightLeft RightRight)
magneticFieldFiles=(LeftUp RightUp)
optimizerConfigFile=OptimizerConfigFiles/Method1_Settings.txt
#beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice.txt
beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice_strippersNotClosed.txt
chicaneScaleDirectory=blank

echo "Method1"
for i in ${magneticFieldFiles[@]}
do
	pyORBIT OptimizerGeneral.py --outputDirectory "$outputDirectory/$outputDirectory"_"$i" --magneticFieldFile "$magneticFieldDirectory/$magneticFieldFilePrefix""$i.txt" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFile" --chicaneScaleDirectory "$chicaneScaleDirectory"
done

