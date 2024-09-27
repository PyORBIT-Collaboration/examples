Author: Kevin Hildebrand
Questions/Comments:krthilde@gmail.com

The get plotable results script:
InjectionRegion_NoSpaceChargeStripNode_General_ArrayConfig.py

An example script that runs script:
RunGeneralTrajectories_Method4_ArrayConfig.sh

This main python code is run separately for each type of beam trajectory (inject Beam, Waste Beam, and closed beam). The example steering script (RunGeneralTrajectories_Method4_ArrayConfig.sh) runs all three, one after the other. For each beam a separate config file is needed.

For our example the three different config files are:
 InjectBeam:PPU_InjectBeam_Method4_ArrayConfig_FloatLength_p80_ULRD.txt
 WasteBeam:PPU_WasteBeam_Method4_ArrayConfig_FloatLength_p80_ULRD.txt
 ClosedBeam:PPU_ClosedBeam_Method4_strippersNotClosed_ArrayConfig_p80_ULRD.txt

A description of what can/needs to be in config file (there is some overlap from the optimizer beam lattice files):
	the MAD lattice file describing the lattice must be defined (it is a comma separated list of MAD files):
		lattice=MAD_Injection_Region_Lattice/PPU_InjectionRegionOnly_Chicane_Replaced_With_Kickers_DB12Laser_To_DB34.LAT,MAD_Injection_Region_Lattice/PPU_InjectionRegionOnly_Chicane_Replaced_With_Kickers_DHA13_To_DB_Waste.LAT
	the nominal chicane kick strength (these values will/can be scaled by the results from the optimizer):
		chicaneKickStrength=.042,-.05011,-.03752,.04563
	whether or not to include custom dipoles (in our example they are excluded from the closed Beam lattice):
		doDipoleStrippers=True
	the number of custom dipoles being added:
		numberOfStripperDipoles=4
	The lines describing each custom dipole (brief description given in Readme_Optimizer.txt):
		stripper1=MagneticFieldFiles/magneticFieldRight.txt,Dipole_DH_A11,True,DH_A11,Before,-1
		stripper2=MagneticFieldFiles/magneticFieldUniformLeftWeak.txt,Dipole_DH_A11B,False,Dipole_DH_A11,Before,-1
		stripper3=MagneticFieldFiles/magneticFieldUniformUpStrong.txt,Dipole_DH_A11C,False,DB12_Laser,.8,-1
		stripper4=MagneticFieldFiles/magneticFieldDown.txt,Dipole_DH_A12,True,DH_A12,Before,0.02	
	Whether or not to use secondary foil (for our example this is only used in waste beam lattice):
		useFoil=False
	The foil is not defined in the MAD file so if using it it must be added separately. You can define which lattice it is added on to the end of in the list of MAD lattice files defined by "lattice" in the config file. For our example the two lattices defined in "lattice" are chosen so they are split around where the secondary foil is so for the waste beam lattice we add it at the end of the first lattice file:
		latticeToAddFoilTo=0
	If you have run the optimzer and have values to read from the resulting chicaneScaleFile then you need to tell the script to read the chicaneScaleFile:
		useChicaneScaleFile=True
	If you are reading the chicaneScaleFile you need to tell the script which values from it are being using. Options will be mentioned but first the format of the chicaneScaleFile. It contains 1 line of comma separated values. The first 4 (0-3) are the corresponding chicane kick strength scales. The next 2 (4-5) are the inject beam initial offset (old and not used in any recent examples). The next 2 (6-7) are the closed beam initial offset (old and not used in any recent examples).The remaining numbers are the length, angle and offset of the custom dipoles in strange but consistent order.
The first beam lattice config file used in the example script will be described (the optimizer and beam lattice config file go hand in hand). 
	To use the chicane kick strength scale from the chicaneScaleFile:
		chicaneScaleFile_UseScales=True
	To use the scales from the chicaneScaleFile for the custom dipoles the format is:
		chicaneScaleFile_useStripperLength3=True
		chicaneScaleFile_useStripperAngle3=False
		chicaneScaleFile_useStripperOffset3=True
	For the closed beam the charge of the bunch never changes, for the inject beam the stripper dipoles handle the charge changing, but for the waste beam the first stripper handles changing the charge from -1 to 0 however the foil does not change the charge. So to properly handle this there is a config option to define at the end of which lattice to change the charge of the bunch to +1. For the waste beam we want this as the same place we placed the secondary foil so:
		latticeIndexToMakeBunchChargePlusOneAtEndOf=0
	For the closed beam and inject beam config we don't want/need to change charge to +1 manually so we set this value to:
		latticeIndexToMakeBunchChargePlusOneAtEndOf=-1
	Initial beam properties are set the same as in the optimizer so for injection and Waste Beam:
		e_kin_ini=1.3
		mass=0.93827231
		initial_charge=-1
		nParts=10000

		#units [m/rad]
		xOffset=0.25671
		pxOffset=-.042
		yOffset=0.046
		pyOffset=0

		alphaZ=0.0196
		betaZ=0.5844
		emittZ=0.24153

		alphaX=.224
		betaX=10.5
		emittX=1.445

		alphaY=.224
		betaY=10.5
		emittY=1.445	
	and for closed beam:
		e_kin_ini=1.3
		mass=0.93827231
		initial_charge=1
		nParts=10000

		#units [m/rad]
		xOffset=0
		pxOffset=0
		yOffset=0
		pyOffset=0

		alphaZ=0.0196
		betaZ=0.5844
		emittZ=0.24153

		alphaX=.224
		betaX=10.5
		emittX=1.445

		alphaY=.224
		betaY=10.5
		emittY=1.445	
	If you want to have the macroparticle bunch parameters to be printed at each node you need to set this in config:
		usePrintNode=True
	The overall bunch info created by calculateEmit.py node is always created by default and can't be turned off.
	How to plot the results is left to the user. I wrote scripts that use ROOT to do the plotting. For the curious the scripts are in a separate project called PythonPlotting. They require ROOT to be installed and configured. The two most interesting plotting scripts are plot_GeneralTrajectories_ArrayConfig.py and plot_phaseSpace.py which are run using RunPlotGeneralTraj_Method2_ArrayConfig.sh and run_plotPhaseSpace.sh, respectively.
	
	
