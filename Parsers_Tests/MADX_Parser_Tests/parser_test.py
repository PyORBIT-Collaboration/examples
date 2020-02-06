import matplotlib.pyplot as plt
import numpy as np

from bunch import Bunch
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
#------------------------------------------------------------------------------------------------------

def load_lattice(file_name,lattice_name,new_parser=False):

	lattice = TEAPOT_Lattice("lattice")

	if new_parser:
		print "Generate Lattice using MADX (new) parser from .SEQ file"
		lattice.readMADX(file_name,lattice_name)
	else:
		print "Generate Lattice using MAD parser from .LAT file"
		lattice.readMAD(file_name,lattice_name)
	lattice.setUseRealCharge(useCharge = 1)
	return lattice

def calculate_betas(teapot_lattice, bunch):
	matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_lattice,bunch)
	TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
	return np.transpose(TwissDataX[-1]), np.transpose(TwissDataY[-1])

#---------------------------------------------Bunch init---------------------------------------------
print "Start."

b = Bunch()
energy = 1
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)

#---------------------------------------------Make a Teapot Lattices from MAD-like file---------------

sis18_mad = load_lattice("sis18.lat","SIS18",False)
sector_mad = load_lattice("YR_sector.lat","yring",False)

#---------------------------------------------Make a Teapot Lattices from MADX-like file-------------- new parser

sis18_madx = load_lattice("sis18.seq","sis18",True)
sector_madx = load_lattice("YR_sector.seq","yring",True)

#------------------------------------------------------------------------------------------------------

sis_betx0,sis_bety0 = calculate_betas(sis18_mad,b)
sis_betx,sis_bety = calculate_betas(sis18_madx,b)

yr_betx0,yr_bety0 = calculate_betas(sector_mad,b)
yr_betx,yr_bety = calculate_betas(sector_madx,b)

#-------------------------------plotting---------------------------------------------

plt.figure()
plt.plot(sis_betx0[0],sis_betx0[1],label = r"$\beta_x$, mad")
plt.plot(sis_bety0[0],sis_bety0[1],label = r"$\beta_y$, mad")

plt.plot(sis_betx[0],sis_betx[1], ls="--",label = r"$\beta_x$, madx")
plt.plot(sis_bety[0],sis_bety[1], ls="--",label = r"$\beta_y$, madx")

plt.title("SIS18 lattice")
plt.xlabel("s [m]")
plt.ylabel(r"$\beta$ [m]")
plt.legend()
#-----
plt.figure()
plt.plot(yr_betx0[0],yr_betx0[1],label = r"$\beta_x$, mad")
plt.plot(yr_bety0[0],yr_bety0[1],label = r"$\beta_y$, mad")

plt.plot(yr_betx[0],yr_betx[1], ls="--",label = r"$\beta_x$, madx")
plt.plot(yr_bety[0],yr_bety[1], ls="--",label = r"$\beta_y$, madx")

plt.title("YRING lattice")
plt.xlabel("s [m]")
plt.ylabel(r"$\beta$ [m]")
plt.legend()
plt.show()


