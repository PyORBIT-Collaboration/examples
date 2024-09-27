#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectory=Method1_UpUp
magneticFieldFile=MagneticFieldFiles/magneticFieldUpUp.txt
optimizerConfigFile=OptimizerConfigFiles/DefaultSettings.txt
beamLatticeFile=OptimizerConfigFiles/DefaultBeamLattice.txt
#beamLatticeFile=OptimizerConfigFiles/BeamLattice_DIfferentClosed.txt
chicaneScaleDirectory=blank

echo "Method1"
pyORBIT OptimizerGeneral.py --outputDirectory "$outputDirectory" --magneticFieldFile "$magneticFieldFile" --optimizerConfigFile "$optimizerConfigFile" --beamLatticeFile "$beamLatticeFile" --chicaneScaleDirectory "$chicaneScaleDirectory"


