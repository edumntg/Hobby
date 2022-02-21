from pyomo.environ import *
import numpy as np

from OPF_Helpers import BuildYbus
from OPF_AC import *
from OPF_DC import*
from OPF_LPAC_Gurobi import*
from math import pi
import random
import logging
from PF_IPOPT import BuildPFModel
import copy
import time

logging.getLogger('pyomo.core').setLevel(logging.ERROR)


Sb = 100 # MVAbase 10MW

#		ID	Type Volt	Angle	Pg	Qg	Pl	Ql

buses = { 1: [1, 0, 1.06, 0.0, 0.0, 0, 0.0, 0.0],
		  2: [2, 1, 1.043, 0.0, 0.0, 0, 21.7/Sb, 12.7/Sb],
		  3: [3, 2, 1.021, 0.0, 0.0, 0, 2.4/Sb, 1.2/Sb],
		  4: [4, 2, 1.012, 0.0, 0.0, 0, 7.6/Sb, 1.6/Sb],
		  5: [5, 2, 1.01, 0.0, 0.0, 0, 94.2/Sb, 19.0/Sb],
		  6: [6, 2, 1.01, 0.0, 0.0, 0, 0.0, 0.0],
		  7: [7, 2, 1.002, 0.0, 0.0, 0, 22.8/Sb, 10.9/Sb],
		  8: [8, 2, 1.01, 0.0, 0.0, 0, 30.0/Sb, 30.0/Sb],
		  9: [9, 2, 1.051, 0.0, 0.0, 0, 0.0, 0.0],
		  10: [10, 2, 1.045, 0.0, 0.0, 0, 5.8/Sb, 2.0/Sb],
		  11: [11, 2, 1.082, 0.0, 0.0, 0, 0.0, 0.0],
		  12: [12, 2, 1.057, 0.0, 0.0, 0, 11.2/Sb, 7.5/Sb],
		  13: [13, 2, 1.071, 0.0, 0.0, 0, 0.0, 0.0],
		  14: [14, 2, 1.042, 0.0, 0.0, 0, 6.2/Sb, 1.6/Sb],
		  15: [15, 2, 1.038, 0.0, 0.0, 0, 8.2/Sb, 2.5/Sb],
		  16: [16, 2, 1.045, 0.0, 0.0, 0, 3.5/Sb, 1.8/Sb],
		  17: [17, 2, 1.040, 0.0, 0.0, 0, 9.0/Sb, 5.8/Sb],
		  18: [18, 2, 1.028, 0.0, 0.0, 0, 3.2/Sb, 0.9/Sb],
		  19: [19, 2, 1.026, 0.0, 0.0, 0, 9.5/Sb, 3.4/Sb],
		  20: [10, 2, 1.03, 0.0, 0.0, 0, 2.2/Sb, 0.7/Sb],
		  21: [21, 2, 1.033, 0.0, 0.0, 0, 17.5/Sb, 11.2/Sb],
		  22: [22, 2, 1.033, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  23: [23, 2, 1.027, 0.0, 0.0, 0, 3.2/Sb, 1.6/Sb],
		  24: [24, 2, 1.021, 0.0, 0.0, 0, 8.7/Sb, 6.7/Sb],
		  25: [25, 2, 1.017, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  26: [26, 2, 1.00, 0.0, 0.0, 0, 3.5/Sb, 2.3/Sb],
		  27: [27, 2, 1.023, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  28: [28, 2, 1.007, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  29: [29, 2, 1.003, 0.0, 0.0, 0, 2.4/Sb, 0.9/Sb],
		  30: [30, 2, 0.992, 0.0, 0.0, 0, 10.6/Sb, 1.9/Sb]}

lineas = { 1: [1, 2, 0.0192, 0.0575, 0.0264, 1, 130/Sb],
		   2: [1, 3, 0.0452, 0.1852, 0.0204, 1, 130/Sb],
		   3: [2, 4, 0.0570, 0.1737, 0.0184, 1, 65/Sb],
		   4: [3, 4, 0.0132, 0.0379, 0.0042, 1, 130/Sb],
		   5: [2, 5, 0.0472, 0.1983, 0.0209, 1, 130/Sb],
		   6: [2, 6, 0.0581, 0.1763, 0.0187, 1, 65/Sb],
		   7: [4, 6, 0.0119, 0.0414, 0.0045, 1, 90/Sb],
		   8: [5, 7, 0.0460, 0.1160, 0.0102, 1, 70/Sb],
		   9: [6, 7, 0.0267, 0.0820, 0.0085, 1, 130/Sb],
		   10: [6, 8, 0.0120, 0.0420, 0.0045, 1, 32/Sb],
		   11: [6, 9, 0.0000, 0.2080, 0.0000, 1.0155, 65/Sb],
		   12: [6, 10, 0.0000, 0.5560, 0.0000, 0.9629, 32/Sb],
		   13: [9, 11, 0.0000, 0.2080, 0.0000, 1, 65/Sb],
		   14: [9, 10, 0.0000, 0.1100, 0.0000, 1, 65/Sb],
		   15: [4, 12, 0.0000, 0.2560, 0.0000, 1.0129, 65/Sb],
		   16: [12, 13, 0.0000, 0.1400, 0.0000, 1, 65/Sb],
		   17: [12, 14, 0.1231, 0.2559, 0.0000, 1, 32/Sb],
		   18: [12, 15, 0.0662, 0.1304, 0.0000, 1, 32/Sb],
		   19: [12, 16, 0.0945, 0.1987, 0.0000, 1, 32/Sb],
		   20: [14, 15, 0.2210, 0.1997, 0.0000, 1, 16/Sb],
		   21: [16, 17, 0.0824, 0.1932, 0.0000, 1, 16/Sb],
		   22: [15, 18, 0.1070, 0.2185, 0.0000, 1, 16/Sb],
		   23: [18, 19, 0.0639, 0.1292, 0.0000, 1, 16/Sb],
		   24: [19, 20, 0.0340, 0.0680, 0.0000, 1, 32/Sb],
		   25: [10, 20, 0.0936, 0.2090, 0.0000, 1, 32/Sb],
		   26: [10, 17, 0.0324, 0.0845, 0.0000, 1, 32/Sb],
		   27: [10, 21, 0.0348, 0.0749, 0.0000, 1, 32/Sb],
		   28: [10, 22, 0.0727, 0.1499, 0.0000, 1, 32/Sb],
		   29: [21, 22, 0.0116, 0.0236, 0.0000, 1, 32/Sb],
		   30: [15, 23, 0.1000, 0.2020, 0.0000, 1, 16/Sb],
		   31: [22, 24, 0.1150, 0.1790, 0.0000, 1, 16/Sb],
		   32: [23, 24, 0.1320, 0.2700, 0.0000, 1, 16/Sb],
		   33: [24, 25, 0.1885, 0.3292, 0.0000, 1, 16/Sb],
		   34: [25, 26, 0.2544, 0.3800, 0.0000, 1, 16/Sb],
		   35: [25, 27, 0.1093, 0.2087, 0.0000, 1, 16/Sb],
		   36: [28, 27, 0.0000, 0.3690, 0.0000, 0.9581, 65/Sb],
		   37: [27, 29, 0.2198, 0.4153, 0.0000, 1, 16/Sb],
		   38: [27, 30, 0.3202, 0.6027, 0.0000, 1, 16/Sb],
		   39: [29, 30, 0.2399, 0.4533, 0.0000, 1, 16/Sb],
		   40: [8, 28, 0.0636, 0.2000, 0.0214, 1, 32/Sb],
		   41: [6, 28, 0.0169, 0.0599, 0.0065, 1, 32/Sb]}
		   
shunts = {1: [10, -0.19],
		  2: [24, -0.043]}
		  
# shunts = {}
		   
gens = {1: [1, 50/Sb, 200/Sb, -0/Sb, 0/Sb, 0.00375, 2.0, 0],
		2: [2, 20/Sb, 80/Sb, -20/Sb, 100/Sb, 0.0175, 1.75, 0],
		3: [5, 15/Sb, 50/Sb, -15/Sb, 80/Sb, 0.0625, 1.0, 0],
		4: [8, 10/Sb, 35/Sb, -15/Sb, 60/Sb, 0.00834, 3.25, 0],
		5: [11, 10/Sb, 30/Sb, -10/Sb, 50/Sb, 0.0250, 3.0, 0],
		6: [13, 12/Sb, 40/Sb, -15/Sb, 60/Sb, 0.0250, 3.0, 0]}


# modelDC = ModelOPF_DC(buses, lineas, gens)
# modelAC = ModelOPF_AC(buses, lineas, gens, shunts)


# solver = SolverFactory('gurobi', solver_io='python')
solver = SolverFactory('gurobi', solver_io='python')
solverPF = SolverFactory('ipopt')

# results = solver.solve(modelDC)
# PrintOPFDCResults(modelDC, buses, lineas, gens)

# results = solver.solve(modelAC)
# PrintOPFACResults(modelAC, buses, lineas, gens, shunts)

Pl_total = sum(buses[i][6] for i in buses)
print(Pl_total)
cases = 100
ndmg = 3

Plsup_mean = 0
successOPF = 0
successPF = 0
for k in range(1, cases+1):

	lineas_ = copy.deepcopy(lineas)
	
	damagedLines = random.sample(lineas_.items(), ndmg)
	for (i, v) in damagedLines:
		lineas_.pop(i)
		
	nl = len(lineas_)	
	model = ModelOPF_LPAC(buses, lineas_, gens, shunts)
	# model = ModelOPF_DC(buses, lineas_, gens)
	
	try:
		result = solver.solve(model)
		print(result.solver.termination_condition)
		if result.solver.termination_condition == TerminationCondition.optimal or result.solver.termination_condition == TerminationCondition.globallyOptimal or result.solver.termination_condition == TerminationCondition.locallyOptimal or result.solver.termination_condition == TerminationCondition.minFunctionValue:
			successOPF += 1
			
			busesPF = copy.deepcopy(buses)
			lineasPF = copy.deepcopy(lineas_)
			shuntsPF = copy.deepcopy(shunts)
			
			for i in buses:
				if busesPF[i][1] == 0: # slack
					busesPF[i][4] = model.Pg[i]()
					busesPF[i][5] = model.Qg[i]()
				elif busesPF[i][1] == 1: # PV
					busesPF[i][3] = model.th[i]()
					busesPF[i][5] = model.Qg[i]()
				elif busesPF[i][1] == 2: # PQ
					busesPF[i][2] += model.fi[i]()
					busesPF[i][3] = model.th[i]()
					
				busesPF[i][6] *= model.l[i]()
				busesPF[i][7] *= model.l[i]()
				
			modelPF = BuildPFModel(busesPF, lineasPF, shuntsPF)
			resultPF = solverPF.solve(modelPF)
			
			print(resultPF.solver.termination_condition)
			if resultPF.solver.termination_condition == TerminationCondition.optimal or resultPF.solver.termination_condition == TerminationCondition.globallyOptimal or resultPF.solver.termination_condition == TerminationCondition.locallyOptimal or resultPF.solver.termination_condition == TerminationCondition.minFunctionValue:
				successPF += 1
				Pl_sup = sum(busesPF[i][6] for i in busesPF)
				perc = Pl_sup/Pl_total
				Plsup_mean += perc
			
		if successOPF > 0 and successPF > 0:
			print('[' + str(k) + '/' + str(cases) + ']: -> OPF: ' + str(successOPF) + ', PF: ' + str(successPF) + ' - Mean: {0:.4f}.'.format(Plsup_mean/successPF *100))
		
	except ValueError:
		print("Exception in case {0:.0f}".format(k))

	# print('[' + str(k) + '/' + str(cases) + ']: Damaged Power Lines: ' + str([k for (k,v) in damagedLines]) + '.')

# PrintOPFLPACResults(modelLPAC, buses, lineas, gens, shunts)
# print(results)