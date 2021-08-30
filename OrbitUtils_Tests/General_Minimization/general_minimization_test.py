#-------------------------------------------------------------------------------
# This a test of a general minimization package
#-------------------------------------------------------------------------------

from orbit.utils.fitting import Solver
from orbit.utils.fitting import Scorer
from orbit.utils.fitting import SolveStopperFactory
from orbit.utils.fitting import VariableProxy
from orbit.utils.fitting import TrialPoint

from orbit.utils.fitting import SimplexSearchAlgorithm
from orbit.utils.fitting import RandomSearchAlgorithm

class MyScorer(Scorer):
	""" The implementation of the abstract Score class """
	def __init__(self):
		Scorer.__init__(self)

	def getScore(self,trialPoint):
		x0 = trialPoint.getVariableProxyArr()[0].getValue()
		x1 = trialPoint.getVariableProxyArr()[1].getValue()
		x2 = trialPoint.getVariableProxyArr()[2].getValue()	
		score = (x0-1.0)**2 + (x1-2.0)**2 + (x2-3.0)**2 
		return score	

scorer = MyScorer()

searchAlgorithm = RandomSearchAlgorithm()
#searchAlgorithm = SimplexSearchAlgorithm()

max_time = 0.05
solverStopper = SolveStopperFactory.maxTimeStopper(max_time)

solver = Solver()
solver.setAlgorithm(searchAlgorithm)
solver.setStopper(solverStopper)

trialPoint = TrialPoint()
trialPoint.addVariableProxy(VariableProxy(name = "x0", value = 0., step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x1", value = 0., step = 0.1))
trialPoint.addVariableProxy(VariableProxy(name = "x2", value = 0., step = 0.1))

solver.solve(scorer,trialPoint)

print "===== best score ========== fitting time =",solver.getScoreboard().getRunTime()
bestScore = solver.getScoreboard().getBestScore()	
print "best score=",bestScore," iteration=",solver.getScoreboard().getIteration()
trialPoint = solver.getScoreboard().getBestTrialPoint()
print trialPoint.textDesciption()

print "=================================="
print "Done."