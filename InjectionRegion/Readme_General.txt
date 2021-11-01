Author: Kevin Hildebrand
Questions/Comments:krthilde@gmail.com

For these scripts to work they require my (khilde on github) fork of the PyOrbit code (specifically my changes to Function, Scorer(optimizer), and teapot, and possibly other changes not thought of).
Important steering scripts:
InjectionRegion_NoSpaceChargeStripNode_General_ArrayConfig.py
OptimizerGeneral_ArrayConfig.py

Important python files used by above scripts:
ConfigureFileClass.py
MagneticFieldClass.py
KevinPython/calculateEmit.py
KevinPython/printNode.py

KevinPython/printNode.py:
	This is a node that prints the macroparticle info to a txt file for many uses. Example: phase space plot
	Format of output file is:
		<macroparticle index> (x,px,y,py,z,pz,s)= (<x>,<px>,<y>,<py>,<z>,<dE>,<s>)
		Example:
			0 (x,px,y,py,z,pz,s)= (0.109542,-0.041435,0.050838,0.000692,-0.001529,0.001916,6.072170)	
	If only 1 turn is taken through lattice then the macroparticle index will be unique for each line in the output file.
	The output file is always appended too. So when more than 1 turn is taken in the lattice the file will be appended and you can see the macroparticle evolve each turn. 
	
KevinPython/calculateEmit.py:
	This is a node that prints bunch info to a txt file for many uses. Examples: the emmitance before or after a node, the X RMS, etc.
	
ConfigureFileClass.py:
	This is just a class used to read in configuration files
	
MagneticFieldClass.py:
	Class used for creating stripper dipole magnets or non stripper dipole magnets.
