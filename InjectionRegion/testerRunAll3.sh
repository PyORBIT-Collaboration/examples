#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam wastebeam
#directory2 is for closed beam
#directory3 is injectbeam ring
directory=WasteBeamSplitGeneralNewStripperChicaneFieldAddedCleanNewTestY_Inject3c_ReverseBoth
directory2=WasteBeamSplitGeneralClosedBeamNewY_Inject3c_ReverseBoth
directory3=InjectBeam3c_ReverseBoth
directoryChicaneScales=InjectBeam3c_ReverseBoth

echo "HERE0"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_AddChicaneToStripperIncludeChicaneFieldForStripPropability.py --outputDirectory "$directory" --chicaneScaleDirectory "$directoryChicaneScales"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_ClosedBeam_SplitChicaneGeneral.py --outputDirectory "$directory2" --chicaneScaleDirectory "$directoryChicaneScales"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_InjectBeam_SplitChicaneGeneral_AddChicaneToStripperIncludeChicaneFieldForStripPropability.py --outputDirectory "$directory3" --chicaneScaleDirectory "$directoryChicaneScales"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/NotPencilDipolesReplaceChicanes" --bunchFromFile "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.06 --stripperLength2 0.06

#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_AddChicaneToStripperIncludeChicaneFieldForStripPropabilityNoStrip_tester.py --outputDirectory "$directory/PencilDipolesReplaceChicanes" --pencilBeam "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.872201172553 --stripperLength2 0.990369023898

