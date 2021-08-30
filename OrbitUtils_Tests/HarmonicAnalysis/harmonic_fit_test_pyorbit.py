#-------------------------------------------------------------------------------
# This a test for harmonic data fitting by using PyORBIT optimization package
# HarmonicData includes includes f(phase) points and order=4 approximation 
# representation with 9 parameters (a0, a1,delta-phase1 a2,delta-phase2 ...)
# In the fitting procedure we find these parameters to reproduce data.
#-------------------------------------------------------------------------------

import sys
import math

from orbit_utils import Function
from orbit_utils import HarmonicData

from orbit.utils.fitting import Solver
from orbit.utils.fitting import Scorer
from orbit.utils.fitting import SolveStopperFactory
from orbit.utils.fitting import VariableProxy
from orbit.utils.fitting import TrialPoint

from orbit.utils.fitting import SimplexSearchAlgorithm

#---- data generation. They are inside the Function f
f = Function()

for ind in range(0,360,30):
	x = 1.0*ind
	y = 2.0*math.cos((math.pi/180.)*(x+ 25.))+ 4.0*math.cos((math.pi/180.)*(4*x+ 35.)) + 0.5
	y_err = 0.01*abs(y)
	f.add(x,y,y_err)
#----------------------------------------------------

#------ harmonic data with assumed order = 4
order = 4

harmonic_data = HarmonicData(order,f)

x_arr = [0.8,2.1,25.2,0.,0.,0.,0.,4.3,35.4]
for x_ind in range(len(x_arr)):
	harmonic_data.parameter(x_ind,x_arr[x_ind])

harmonic_data.parameter(0,0.5)
harmonic_data.parameter(1,2.1)
harmonic_data.parameter(2,25.2)
harmonic_data.parameter(7,4.3)
harmonic_data.parameter(8,35.4)

#---- class to provide difference between data and fit function
class HarmonicFitFunction(Scorer):
	"""
	The implementation of the abstract Score class 
	as harmonic function score
	"""
	def __init__(self,harmonic_data):
		self.harmonic_data = harmonic_data
	
	def getScore(self,trialPoint):
		x_arr = trialPoint.getVariableProxyValuesArr()
		for x_ind in range(len(x_arr)):
			self.harmonic_data.parameter(x_ind,x_arr[x_ind])
		return harmonic_data.sumDiff2()		

#---- Initial point initilization 
variableProxy_arr = []
variableProxy_arr.append(VariableProxy("amp0", 0.8,0.5))
variableProxy_arr.append(VariableProxy("amp1", 2.1,0.5))
variableProxy_arr.append(VariableProxy("pha1",22.0,0.5))
variableProxy_arr.append(VariableProxy("amp2", 0.0,0.5))
variableProxy_arr.append(VariableProxy("pha2", 0.0,0.5))
variableProxy_arr.append(VariableProxy("amp3", 0.0,0.5))
variableProxy_arr.append(VariableProxy("pha3", 0.0,0.5))
variableProxy_arr.append(VariableProxy("amp4", 3.0,0.5))
variableProxy_arr.append(VariableProxy("pha4",30.0,0.5))

trialPoint = TrialPoint()
for variableProxy in variableProxy_arr:
	trialPoint.addVariableProxy(variableProxy)

#---- Search algorithm from PyORBIT native package
searchAlgorithm = SimplexSearchAlgorithm()

#maxIter = 1000
#solverStopper = SolveStopperFactory.maxIterationStopper(maxIter)

max_time = 0.04
solverStopper = SolveStopperFactory.maxTimeStopper(max_time)

solver = Solver()
solver.setAlgorithm(searchAlgorithm)
solver.setStopper(solverStopper)

scorer = HarmonicFitFunction(harmonic_data)

solver.solve(scorer,trialPoint)	

#---- the fitting process ended, now about results
print "=============================================================="
print "??????????????????????????????????????????????????????????????"
solver.getScoreboard().printScoreBoard()
print "===== best score ========== fitting time =",solver.getScoreboard().getRunTime()
bestScore = solver.getScoreboard().getBestScore()	
print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
trialPoint = solver.getScoreboard().getBestTrialPoint()
print trialPoint.textDesciption()	

#----- this will set the trial point for best score to the harmonic_data
best_score = scorer.getScore(trialPoint)

print "======= fitting quality ================"
for ind in range(harmonic_data.dataSize()):
	x = harmonic_data.valueX(ind)
	y = harmonic_data.valueY(ind)
	y_err = harmonic_data.valueErr(ind)
	y_fit = harmonic_data.fitValueY(x)
	print "ind= %3d "%ind, " (x,y+-y_err,y_fit) = ( %8.1f , %+8.5f +- %8.5f ,%8.5f)"%(x,y,y_err,y_fit)
	
	