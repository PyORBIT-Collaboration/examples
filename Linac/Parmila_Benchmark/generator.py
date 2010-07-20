#!/usr/bin/env python

#--------------------------------------------------------
# The script will generates bunches for Parmila and Impact
# codes the by using gerenrators from bunch_generators.py
#--------------------------------------------------------

import math
import sys
import os
import random

from distribution_generators import TwissContainer
from distribution_generators import KVDist2D, KVDist3D
from distribution_generators import GaussDist2D, GaussDist3D
from distribution_generators import WaterBagDist2D, WaterBagDist3D
from distribution_generators import TwissAnalysis

random.seed(1)

n_particles = 20000   # number of particles
e_kin = 2.5          # MeV - kinetic energy
mass = 939.3014      # MeV - mass of H-
c = 2.99792458e+10   # speed of light in cm/sec
freq = 402.5e+6      # cavity frequency in Hz
beam_current = 0.01  # beam current in mA , design = 38 mA 
lambd = c/freq       # lambda in cm
gamma = 1.0 + e_kin/mass
beta = math.sqrt(1.0 - 1.0/gamma**2)

#=======INITIAL PHASE in DEG
phase_init_deg = -45.0

#-----------x plane ----------
# for x and y axises: alpha - dimensionless, beta - [cm/radian], emitt - [cm*radian]
alpha_x = -1.9619
beta_x = 18.3140
emitt_x = 0.00002/(beta*gamma)

alpha_x = -0.1149*25
beta_x = 55.8710*4
emitt_x = (0.0000160/(beta*gamma))/1

alpha_x = 4.
beta_x = 200.
emitt_x = (0.000016/(beta*gamma))/1


#-----------y plane ----------
# for x and y axises: alpha - dimensionless, beta - [cm/radian], emitt - [cm*radian]
alpha_y = -1.7681
beta_y = 16.1030
emitt_y = 0.00002/(beta*gamma)

alpha_y = -0.0933*25
beta_y = 13.4717*4
emitt_y = (0.0000162/(beta*gamma))/1.

alpha_y = -4.
beta_y = 200.
emitt_y = (0.000016/(beta*gamma))/1

#Distribution contains 20000 particles, current=    0.0100 mA, freq= 402.500 MHz
#Beam distribution parameters: STRANGE
#        rms(n)  100% ellipse  alfa    beta(u)
#        (cm-mr)    (cm-mr)   (cm/rad),(deg/MeV)
#    1    0.0225    0.4853    4.9469  209.1178
#    2    0.0200    0.1574   -0.0504   10.8433
#    3    0.1312    1.2519    0.1884  207.6669

#----------z plane ----------
# for z axis: alpha - dimensionless, beta - [degree/MeV], emitt - [MeV*degree]
# The numbers are from the Parmila file: they are related to the output
# lambda = c*T = c/freq
# beta_z[Parmila_out in degrees/Mev] = beta_z[Parmila_input in cm/radian]*360/(m*c^2*(beta*gamma)^3*lambda)
# emitt_z[Parmila_out, normalized, deg*MeV] = emitt_z[Parmila_in, cm*radian]*360*m*c^2*beta*gamma^3/lambda
alpha_z = 0.0196
beta_z = 58.2172*360/(mass*(beta*gamma)**3*lambd)
emitt_z = 0.1281

alpha_z = 0.0665 
beta_z = 388.0340
emitt_z = 0.1004

alpha_z = 0.1884
beta_z = 207.6669
emitt_z = 0.10

#Beam distribution parameters:
#        rms(n)  100% ellipse  alfa    beta(u)
#        (cm-mr)    (cm-mr)   (cm/rad),(deg/MeV)
#    1    0.0160    0.1269   -0.1149   55.8710
#    2    0.0162    0.1240   -0.0933   13.4717
#    3    0.1004    0.7741    0.0665  388.0340


"""WaterBagDist3D
Parmila input
;input -2 100000  -1.9619 18.3140    0.0021824 ! Toutatis 38 mA emittance to 7mm into MEBT
;                  1.7681 16.1030    0.0021856
;                  0.0196 58.2172    0.003088

;input -2 100000  -1.9899  19.6360   0.0017573 ! exact match to rms properties of ReadDist distribution
;                  1.92893 17.778    0.0017572
;                  0.015682 67.0939  0.002420
"""

#---------------------------------------------
# Set up Twiss for X,Y,Z
#---------------------------------------------
twissX = TwissContainer(alpha = alpha_x, beta = beta_x, emittance = emitt_x)
twissY = TwissContainer(alpha = alpha_y, beta = beta_y, emittance = emitt_y)
twissZ = TwissContainer(alpha = alpha_z, beta = beta_z, emittance = emitt_z)
print "-------------------input parameters----------------------------------------"
print "X alpha= %12.5g    beta [cm/rad] =%12.5g    emitt[cm*rad] = %12.5g "%(alpha_x,beta_x,emitt_x)
print "Y alpha= %12.5g    beta [cm/rad] =%12.5g    emitt[cm*rad] = %12.5g "%(alpha_y,beta_y,emitt_y)
print "Z alpha= %12.5g    beta [deg/MeV] =%12.5g   emitt[deg*MeV] = %12.5g "%(alpha_z,beta_z,emitt_z)

#distributor = GaussDist3D(twissX,twissY,twissZ,cut_off = 4.0)
distributor = WaterBagDist3D(twissX,twissY,twissZ)
#distributor = KVDist3D(twissX,twissY,twissZ)

twiss_analysis = TwissAnalysis(3)

#----------------------------------
#open file for Parmila
#----------------------------------
parmila_out = open("parmila_bunch_tmp.txt","w")
parmila_out.write("Parmila data from *****  Generated ex= %5.4f  ey= %5.4f ez = %6.5f           \n"%(emitt_x,emitt_y,emitt_z))
parmila_out.write("Structure number       =          1 \n")
parmila_out.write("Cell or element number =          0 \n")
parmila_out.write("Design particle energy =%11.6g     MeV \n"%e_kin)
parmila_out.write("Number of particles    =%11d           \n"%n_particles)
parmila_out.write("Beam current           =%11.7f         \n"%beam_current)
parmila_out.write("RF Frequency           =   402.5000     MHz \n")
parmila_out.write("Bunch Freq             =   402.5000     MHz \n")
parmila_out.write("Chopper fraction       =   0.680000  \n")   
parmila_out.write("The input file particle coordinates were written in double precision. \n")
parmila_out.write("   x(cm)             xpr(=dx/ds)       y(cm)             ypr(=dy/ds)       phi(radian)        W(MeV) \n")
	 
	 
#----------------------------------
#open file for Impact
#----------------------------------	 
impact_out = open("impact.dat","w")
impact_coeff_x = 1/(c/(2*math.pi*freq))
impact_coeff_xp = gamma*beta
impact_coeff_phi = math.pi/180
impact_coeff_dE = 1./mass
impact_q_m = -1.0/(mass*1.0e+6)
impact_macrosize = (beam_current*1.0e-3/freq)/n_particles
	 
i = 0
while(i < n_particles):
	i = i + 1
	if(i % 10000 == 0): print "i=",i
	#results are ([m],[rad],[m],[rad],[deg],[MeV])
	(x,xp,y,yp,phi,dE) = distributor.getCoordinates()
	twiss_analysis.account((x,xp,y,yp,phi + phase_init_deg,dE))
	parmila_out.write("%18.11g%18.11g%18.11g%18.11g%18.11g%18.11g \n"%(x,xp,y,yp,(phi + phase_init_deg)*math.pi/180.,dE+e_kin))
	#---write impact line --------
	s = ""
	s = s + " %14.7e"%(x*impact_coeff_x)
	s = s + " %14.7e"%(xp*impact_coeff_xp)
	s = s + " %14.7e"%(y*impact_coeff_x)
	s = s + " %14.7e"%(yp*impact_coeff_xp)
	s = s + " %14.7e"%(phi*impact_coeff_phi)
	s = s + " %14.7e"%(dE*impact_coeff_dE)
	s = s + " %14.7e"%(impact_q_m)
	s = s + " %14.7e"%(impact_macrosize)
	s = s + " %14.7e"%(1.0*i)
	impact_out.write(s+"\n")
	
parmila_out.close()
impact_out.close()

print "n total =",n_particles

#------------------------------
# print the parameters of the generated distribution
#------------------------------
bg = gamma*beta
(alpha_x,beta_x,gamma_x,emitt_x) = twiss_analysis.getTwiss(0)
(alpha_y,beta_y,gamma_y,emitt_y) = twiss_analysis.getTwiss(1)
(alpha_z,beta_z,gamma_z,emitt_z) = twiss_analysis.getTwiss(2)
emitt_x = 1000.*emitt_x*bg
emitt_y = 1000.*emitt_y*bg
print "-------------------bunch's twiss parameters----------------------------------------"
print "X alpha= %12.5g    beta [cm/rad]  =%12.5g    gamma = %12.5g  norm. emitt[cm*mrad]  = %12.5g "%(alpha_x,beta_x,gamma_x,emitt_x)
print "Y alpha= %12.5g    beta [cm/rad]  =%12.5g    gamma = %12.5g  norm. emitt[cm*mrad]  = %12.5g "%(alpha_y,beta_y,gamma_y,emitt_y)
print "Z alpha= %12.5g    beta [deg/MeV] =%12.5g    gamma = %12.5g  emitt[deg*MeV] = %12.5g "%(alpha_z,beta_z,gamma_z,emitt_z)

print "-------------------centroid params--------------------------------------"
print "X x_avg [cm]    =%12.5g    xp_avg [rad]  =%12.5g "%twiss_analysis.getAvgU_UP(0)
print "Y y_avg [cm]    =%12.5g    yp_avg [rad]  =%12.5g "%twiss_analysis.getAvgU_UP(1)
print "Z phi_avg [deg] =%12.5g    Ek [Mev]      =%12.5g "%twiss_analysis.getAvgU_UP(2)

(x_rms,xp_rms) =  twiss_analysis.getRmsU_UP(0)
(y_rms,yp_rms) =  twiss_analysis.getRmsU_UP(1)
(z_rms,zp_rms) =  twiss_analysis.getRmsU_UP(2)
print "-------------------Rms--------------------------------------"
print "X x_rms [cm]    =%12.5g  xp_rms [deg] =%12.5g "%(x_rms,xp_rms*180./math.pi)
print "Y y_rms [cm]    =%12.5g  yp_rms [deg] =%12.5g "%(y_rms,yp_rms*180./math.pi)
print "Z phi_rms [deg] =%12.5g  Ek_rms [MeV] =%12.5g "%(z_rms,zp_rms)

(x_max,xp_max) =  twiss_analysis.getMaxU_UP(0)
(x_min,xp_min) =  twiss_analysis.getMinU_UP(0)
(y_max,yp_max) =  twiss_analysis.getMaxU_UP(1)
(y_min,yp_min) =  twiss_analysis.getMinU_UP(1)
(z_max,zp_max) =  twiss_analysis.getMaxU_UP(2)
(z_min,zp_min) =  twiss_analysis.getMinU_UP(2)
print "-------------------Min Max--------------------------------------"
print "X x_min_max [cm]    =%12.5g %12.5g  xp_min_max [rad] =%12.5g %12.5g "%(x_min,x_max,xp_min,xp_max)
print "Y y_min_max [cm]    =%12.5g %12.5g  yp_min_max [rad] =%12.5g %12.5g "%(y_min,y_max,yp_min,yp_max)
print "Z phi_min_max [deg] =%12.5g %12.5g  Ek_min_max [MeV] =%12.5g %12.5g "%(z_min,z_max,zp_min,zp_max)

print "-------------------generate parmila binary dst file----------------------------------------"

os.system("del part_rfq.dst")
os.system("del part_rfq.txt")
os.system("PARMILA_TXT_2_DST.exe parmila_bunch_tmp.txt")
os.system("unfseq_lf95_to_lf90.exe part_rfq_lf95.dst part_rfq.dst")
os.system("readdst.exe part_rfq.dst part_rfq.txt")
os.system("del part_rfq_lf95.dst")
os.system("del parmila_bunch_tmp.txt")
print "Done! Stop."

