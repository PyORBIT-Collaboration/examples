#!/bin/bash
echo "Hello World"
#filename=outputAddMagnetNoFoilNode35_Parts_3_XDipole_10
#filename=outputAddMagnetNoFoilNode35_Parts_3_Dipole_10
#directory is for injectbeam
#directory2 is for closed beam
directory=tester



echo "HERE0"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_OnlyChicane_tester.py --outputDirectory "$directory/OnlyChicanePencilNoDipoles" --pencilBeam "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_OnlyChicane_tester.py --outputDirectory "$directory/OnlyChicanePencilDipolesReplaceChicanes" --pencilBeam "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperLength1 0.87220117255

echo "HERE0b"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_OnlyChicane_tester.py --outputDirectory "$directory/OnlyChicaneNewBunchEntireNotPencilNoDipoles" --usePrintNode "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_OnlyChicane_tester.py --outputDirectory "$directory/OnlyChicaneNewBunchEntireNotPencilDipolesReplaceChicanes" --bunchFromFile "True" --bunchFromFileName "$directory/OnlyChicaneNewBunchEntireNotPencilNoDipoles/print_beg_0.txt" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperLength1 0.87220117255

pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_OnlyChicane_tester.py --outputDirectory "$directory/OnlyChicanePencilEntireNotPencilDipolesReplaceChicanes_CancelTest" --doDipoleKickers "True" --pencilBeam "True" --cancelFieldTest "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_OnlyChicane_tester.py --outputDirectory "$directory/OnlyChicaneNewBunchEntireNotPencilDipolesReplaceChicanes_CancelTest" --bunchFromFile "True" --bunchFromFileName "$directory/OnlyChicaneNewBunchEntireNotPencilNoDipoles/print_beg_0.txt" --doDipoleKickers "True" --cancelFieldTest "True"

echo "HERE1"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_tester.py --outputDirectory "$directory/PencilNoDipoles" --pencilBeam "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_tester.py --outputDirectory "$directory/PencilDipolesReplaceChicanes" --pencilBeam "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.87220117255 --stripperLength2 0.99036902389

echo "HERE2"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_tester.py --outputDirectory "$directory/NotPencilNoDipoles" --bunchFromFile "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_tester.py --outputDirectory "$directory/NotPencilDipolesReplaceChicanes" --bunchFromFile "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.87220117255 --stripperLength2 0.99036902389

echo "HERE3"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_tester.py --outputDirectory "$directory/NewBunchNotPencilNoDipoles" --usePrintNode "True" 
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplaceEntireChicaneWithStripper_tester.py --outputDirectory "$directory/NewBunchNotPencilDipolesReplaceChicanes" --usePrintNode "True" --bunchFromFile "True" --bunchFromFileName "$directory/NewBunchNotPencilNoDipoles/print_beg_0.txt" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.87220117255 --stripperLength2 0.99036902389

echo "HERE4"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/PartPencilNoDipoles" --pencilBeam "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/PartPencilDipolesReplaceChicanes" --pencilBeam "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.06 --stripperLength2 0.06

echo "HERE5"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/PartNotPencilNoDipoles" --bunchFromFile "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/PartNotPencilDipolesReplaceChicanes" --bunchFromFile "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.06 --stripperLength2 0.06

echo "HERE6"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/NewBunchPartNotPencilNoDipoles" --usePrintNode "True"
pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/NewBunchPartNotPencilDipolesReplaceChicanes" --bunchFromFile "True" --bunchFromFileName "$directory/NewBunchPartNotPencilNoDipoles/print_beg_0.txt" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.06 --stripperLength2 0.06
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/NotPencilNoDipoles" --bunchFromFile "True"
#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_ReplacePartChicaneWithStripper_tester.py --outputDirectory "$directory/NotPencilDipolesReplaceChicanes" --bunchFromFile "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.06 --stripperLength2 0.06

#pyORBIT InjectionRegion_NoSpaceChargeStripNode_WasteBeam_SplitChicaneGeneral_AddChicaneToStripperIncludeChicaneFieldForStripPropabilityNoStrip_tester.py --outputDirectory "$directory/PencilDipolesReplaceChicanes" --pencilBeam "True" --doDipoleKickers "True" --stripperStrengthMax1 0 --stripperStrengthMin1 0 --stripperStrengthMax2 0 --stripperStrengthMin2 0 --stripperLength1 0.872201172553 --stripperLength2 0.990369023898

