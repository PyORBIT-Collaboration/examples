Author: Kevin Hildebrand
Questions/Comments:krthilde@gmail.com

The optimizer script:
OptimizerGeneral_ArrayConfig.py

An example script that runs in this script (probably most complicated one):
RunGeneralOptimizer_Method4_ArrayConfig_Separate.sh

Two config files are used to define settings. They are located in the directory "OptimizerConfigFiles". One defines what is being floated/fixed in the optimizer and what variables (if any) are being read in from the chicane scale file (ie variables which were determined in previous steps of the optimization and need to be read in). The other config file defines the lattices that are being used.

The first beam lattice config file used in the example script will be described (the optimizer and beam lattice config file go hand in hand) (PPU_DefaultBeamLattice_strippersNotClosed_ArrayConfig_strong_Early_first_p80_ULRD.txt):
	Four different MAD lattice files are defined:
		1)latticeInjection
		2)latticeInjectionFull
		3)latticeClosedCompareToInjection
		4)latticeClosed
	How these are used will be described when talking about optimizer config file.
	The chicane kick strengths are defined in config file with:
		chicaneKickStrength=.042,-.05011,-.03752,.04563
		These are for PPU.
	The custom dipoles are not defined in the MAD files these are added separated. To include them in the injection lattice use:
		"doDipoleStrippersInjection=True"
	To include them in the Closed lattice:
		"doDipoleStrippersClosed=True"
	In the example config they are excluded form the closed lattice:
		"doDipoleStrippersClosed=False"
	If they are included in either lattice the number of custom dipole magnets being added must be defined in the config file with the following line:
		"numberOfStripperDipoles=4"
	Then there must be a separate line for each of these custom dipoles that defines its magnetic field,  where it is located,whether or not it can strip and how it strips (the optimizer code has not yet had added the ability to strip using excited lifetimes). The general definition for each stripper is as follows:
		stripperX=["<magnetic field file>","<Dipole Name>",isStripper<True/False>,<Name of node to place it in reference too>,<position with respect to reference node>,<fixed stripper length>]
		<position with respect to reference node> can be "before"(places it right before reference), "after"(places it right after reference), or float [0-1](places it this fraction into reference),<fixed stripper length> sets stripping length same for all particles if >0, if <0 then does nothing
	So in our example config file we have:
		stripper1=MagneticFieldFiles/magneticFieldRight.txt,Dipole_DH_A11,True,DH_A11,Before,-1
		stripper2=MagneticFieldFiles/magneticFieldUniformLeftWeak.txt,Dipole_DH_A11B,False,Dipole_DH_A11,Before,-1
		stripper3=MagneticFieldFiles/magneticFieldUniformUpStrong.txt,Dipole_DH_A11C,False,DB12_Laser,.8,-1
		stripper4=MagneticFieldFiles/magneticFieldDown.txt,Dipole_DH_A12,True,DH_A12,Before,0.02	
	Also included in this file is the initial beam configuration. So for example:
		e_kin_ini=1.3
		mass=0.93827231
		initial_chargeInjection=-1
		nPartsInjection=10000	
	As well as the initial offsets for the Inject beam and other beam properties:
		xOffsetInjection=0.25671
		pxOffsetInjection=-.042
		yOffsetInjection=0.046
		pyOffsetInjection=0

		alphaZInjection=0.0196
		betaZInjection=0.5844
		emittZInjection=0.24153

		alphaXInjection=.224
		betaXInjection=10.5
		emittXInjection=1.445

		alphaYInjection=.224
		betaYInjection=10.5
		emittYInjection=1.445
	These settings can also be set for the closed beam:
		initial_chargeClosed=1
		#currently this does nothing. its always pencil for closed
		pencilBeamClosed=True
		#if pencil beam closed is true then nPartsClosed is irrelevant
		nPartsClosed=10000		
		xOffsetClosed=0
		pxOffsetClosed=0
		yOffsetClosed=0
		pyOffsetClosed=0

		#these values are irrelevent when using pencil beam
		alphaZClosed=0.0196
		betaZClosed=0.5844
		emittZClosed=0.24153

		alphaXClosed=.224
		betaXClosed=10.5
		emittXClosed=1.445

		alphaYClosed=.224
		betaYClosed=10.5
		emittYClosed=1.445	


The first optimizer config file used in the example script will be described (Method4_Part1_Settings_ArrayConfig_FloatLength.txt):		
	1) and 3) are compared for the options "useParallelScore","useParallelScoreY","useParallelScoreX","usePartOffsetDifferenceX"
	4) is used to determine if the closed orbit is closed if the value of "usedClosedScore" is set to true in config file, ie it enters the lattice defined by 4) at (x,py,y,py)=(0,0,0,0)  and exits the lattice with (x,py,y,py)=(0,0,0,0). The target value it exits with can be changed from 0 by using the following variables in the config file: target_xOffsetClosed=0,target_pxOffsetClosed=0,target_yOffsetClosed=0,target_pyOffsetClosed=0
	2) and 4) are compared for the option "useFullOffsetDifferenceX=True" in which case then this option must be set "targetFullOffsetDifferenceX=0.05"
	
	these options allow the scaling of the chicane kicks to be floated or fixed in the optimation:
		fixChicaneScale10=True
		fixChicaneScale11=True
		fixChicaneScale12=True
		fixChicaneScale13=True
	These options tell the program whether or not to read the chicane scales from the chicaneScaleFile.
		readChicaneScale10FromFile=False
		readChicaneScale11FromFile=False
		readChicaneScale12FromFile=False
		readChicaneScale13FromFile=False
	There are also options to fix or float the initial injection and closed beams initial angles in the fit:
		fixInitialPXInjection=True
		fixInitialPYInjection=True
		fixInitialPXClosed=True
		fixInitialPYClosed=True	
	And corresponding options to read them from the chicaneScaleFile:
		readInitialPXInjectionFromFile=False
		readInitialPYInjectionFromFile=False
		readInitialPXClosedFromFile=False
		readInitialPYClosedFromFile=False	
	There are also options	for floating or fixing properties of the stripping magnets. If they have been modified in previous steps or being floated these two flags must be set to True (ugly was ment to be improved):
		readChicaneFile=True	
		modifyAStripper=True
	Then the options to float or fix the custom dipoles:
		floatStripperLength1=False
		floatStripperAngle1=False
		floatStripperOffset1=False

		floatStripperLength2=False
		floatStripperAngle2=False
		floatStripperOffset2=False

		floatStripperLength3=False
		floatStripperAngle3=False	
		floatStripperOffset3=False
		
		floatStripperLength4=False
		floatStripperAngle4=False	
		floatStripperOffset4=False	
	If these have been floated in previous optimization steps then to read in the values from the chicaneScaleFile the options used are as follows:
		chicaneScaleFile_useStripperLength2=True
		chicaneScaleFile_useStripperAngle2=True		
		chicaneScaleFile_useStripperOffset2=True
	
	
Going through the full example "RunGeneralOptimizer_Method4_ArrayConfig_Separate.sh":
	There are 4 custom dipoles in the full injection lattice. Prior to chicane 2, is where 3 are located with there fields in the direction Up, Left and Right. The first two are non stripping dipoles. And the Right one is a stripper dipole. The 4th custom dipole is located right before chicane 3 and is a stripper with direction Down.
	It runs the optimizer script (OptimizerGeneral_ArrayConfig.py) 4 times.
	
	During the 1st time (Method4_Part1_Settings_ArrayConfig_FloatLength.txt):
		The score is used to get the closed beam and injection beam parallel in Y by floating just the length of the second custom dipole.
	During the 2nd time (Method4_Part2_Settings_ArrayConfig_FloatLength.txt):
		The new length scale of 2nd custom dipole from first part is read in from chicane scale file.
		The first and second chicane kick strength are floated and the length of the 3rd custom dipole is floated to get the inject beam and closed beam separated by .05 m after chicane 2 as well as parallel in X after chicane 4.
	During the 3rd time (Method4_Part3_Settings_ArrayConfig_FloatLength.txt):
		The values floated during part 1 and part 2 are read in from the chicaneScaleFile and Fixed.
		The 3rd and 4th chicane kick strengths are floated to get the closed beam closed (ie enters with (x,px,y,py)=0 and exits at 0) just before the quad after chicane 4.
	During the 4th time (Method4_Part4_Settings_ArrayConfig_FloatLength.txt):
		The values floated during part 1, part 2 and part 3 are read in from the chicaneScaleFile and Fixed.
		The offset of the 3rd custom dipole is floated to get the inject and closed beam separated by .05 m after chicane 4.
		
	The net result of all 4 steps results in the closed beam being closed and the inject beam being parallel to the closed beam after chicane 4 in both X and Y and offset from the closed beam by .05 m in X.
	
