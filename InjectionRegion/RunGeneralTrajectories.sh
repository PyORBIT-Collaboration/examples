#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
outputDirectoryPrefix=Method1_UpUp
magneticFieldFile=MagneticFieldFiles/magneticFieldUpUp.txt
beamLatticeFileDirectory=BeamLatticeFiles
chicaneScaleDirectory=Method1_UpUp
suffixInjection=InjectBeam
suffixWaste=WasteBeam
suffixClosed=ClosedBeam

suffix=(InjectBeam WasteBeam ClosedBeam)

#echo "$outputDirectoryPrefix"_"$suffixInjection"
echo "Method1"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectoryPrefix"_"$suffixInjection" --magneticFieldFile "$magneticFieldFile" --beamLatticeFile "$beamLatticeFileDirectory"/"$suffixInjection" --chicaneScaleDirectory "$chicaneScaleDirectory"

for i in ${suffix[@]}
do
	pyORBIT InjectionRegion_NoSpaceChargeStripNode_General.py --outputDirectory "$outputDirectoryPrefix"_"$i" --magneticFieldFile "$magneticFieldFile" --beamLatticeFile "$beamLatticeFileDirectory"/"$i.txt" --chicaneScaleDirectory "$chicaneScaleDirectory"
done