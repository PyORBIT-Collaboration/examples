#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
directoryChicaneScales="ClosureTest"

echo "HERE0"
pyORBIT injectionRegionOnly_OptimizerStripperSplitChicaneInjectBeamClass.py 
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_ClosedBeam_SplitChicaneGeneral.py --outputDirectory "$directoryChicaneScales" --chicaneScaleDirectory "$directoryChicaneScales"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/NotPencilDipolesReplaceChicanes" --bunchFromFile "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.06 --stripperLength2 0.06

#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_AddChicaneToStripperIncludeChicaneFieldForStripPropabilityNoStrip_tester.py --outputDirectory "$directory/PencilDipolesReplaceChicanes" --pencilBeam "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.872201172553 --stripperLength2 0.990369023898

