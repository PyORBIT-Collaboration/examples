import sys
import os

folder1="WasteBeamSplitGeneralNewStripperChicaneFieldAddedOG"
folder2="WasteBeamSplitGeneralNewStripperChicaneFieldAddedClean"

#check directories exist
if not os.path.isdir(folder1):
	print "folder1= ",folder1," does not exist"
if not os.path.isdir(folder2):
	print "folder2= ",folder2," does not exist"	
	
#loop over files in folder1
for currentFile in os.listdir(folder1):
	if "emmit" == currentFile[0:5]:
		openedFile1=open("%s/%s"%(folder1,currentFile))
		#check file exists in folder2
		openedFile2=None
		if (not os.path.isfile("%s/%s"%(folder2,currentFile))):
			print "currrentFile= ", currentFile," does not exist in folder2= ",folder2
			sys.exit(0)
		else:
			openedFile2=open("%s/%s"%(folder2,currentFile))
		#check that values in files are the same
		lines1=openedFile1.readlines()
		lines2=openedFile2.readlines()
		numberOfLinesToCheck=len(lines1)
		if len(lines1)!=len(lines2):
			print "currrentFile= ", currentFile, " has len(lines1)= ",len(lines1), " and len(lines2)= ",len(lines2)
			if len(lines2)<len(lines1):
				numberOfLinesToCheck=len(lines2)
		for currentLineN in range(numberOfLinesToCheck):
			if lines1[currentLineN].strip()!=lines2[currentLineN].strip():
				print "currrentFile= ", currentFile, " lines1[",currentLineN,"]= ", lines1[currentLineN], " != lines2[",currentLineN,"]= ",lines2[currentLineN]
		