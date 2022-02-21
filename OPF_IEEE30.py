from pyomo.environ import *
import numpy as np

from OPF_DC import*
from OPF_LPAC_Gurobi import*
from math import pi
import random
import logging
# import scipy as sp
# import csv
# import matplotlib.pyplot as plt
import math
import copy

from PF_IPOPT import BuildPFModel

logging.getLogger('pyomo.core').setLevel(logging.ERROR)

FILENAME = 'cases-3dmg.csv'

solver = SolverFactory('gurobi', solver_io='python')

Sb = 100 # MVAbase 10MW

#		ID	Type Volt	Angle	Pg		Qg		Pl		Ql

buses = { 1: [1, 0, 1.00, 0.0, 0.0, 0, 0.0, 0.0],				# slack
		  2: [2, 1, 1.01, 0.0, 0.0, 0, 21.7/Sb, 12.7/Sb],		# PV
		  3: [3, 2, 1.00, 0.0, 0.0, 0, 2.4/Sb, 1.2/Sb],
		  4: [4, 2, 1.00, 0.0, 0.0, 0, 7.6/Sb, 1.6/Sb],
		  5: [5, 1, 1.00, 0.0, 0.0, 0, 94.2/Sb, 19.0/Sb],		# PV
		  6: [6, 2, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
		  7: [7, 2, 1.00, 0.0, 0.0, 0, 22.8/Sb, 10.9/Sb],
		  8: [8, 1, 1.00, 0.0, 0.0, 0, 30.0/Sb, 30.0/Sb],		# PV
		  9: [9, 2, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
		  10: [10, 2, 1.00, 0.0, 0.0, 0, 5.8/Sb, 2.0/Sb],
		  11: [11, 1, 1.00, 0.0, 0.0, 0, 0.0, 0.0],				# PV
		  12: [12, 2, 1.00, 0.0, 0.0, 0, 11.2/Sb, 7.5/Sb],
		  13: [13, 1, 1.00, 0.0, 0.0, 0, 0.0, 0.0],				# PV
		  14: [14, 2, 1.00, 0.0, 0.0, 0, 6.2/Sb, 1.6/Sb],
		  15: [15, 2, 1.00, 0.0, 0.0, 0, 8.2/Sb, 2.5/Sb],
		  16: [16, 2, 1.00, 0.0, 0.0, 0, 3.5/Sb, 1.8/Sb],
		  17: [17, 2, 1.00, 0.0, 0.0, 0, 9.0/Sb, 5.8/Sb],
		  18: [18, 2, 1.00, 0.0, 0.0, 0, 3.2/Sb, 0.9/Sb],
		  19: [19, 2, 1.00, 0.0, 0.0, 0, 9.5/Sb, 3.4/Sb],
		  20: [10, 2, 1.00, 0.0, 0.0, 0, 2.2/Sb, 0.7/Sb],
		  21: [21, 2, 1.00, 0.0, 0.0, 0, 17.5/Sb, 11.2/Sb],
		  22: [22, 2, 1.00, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  23: [23, 2, 1.00, 0.0, 0.0, 0, 3.2/Sb, 1.6/Sb],
		  24: [24, 2, 1.00, 0.0, 0.0, 0, 8.7/Sb, 6.7/Sb],
		  25: [25, 2, 1.00, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  26: [26, 2, 1.00, 0.0, 0.0, 0, 3.5/Sb, 2.3/Sb],
		  27: [27, 2, 1.00, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  28: [28, 2, 1.00, 0.0, 0.0, 0, 0.0/Sb, 0.0],
		  29: [29, 2, 1.00, 0.0, 0.0, 0, 2.4/Sb, 0.9/Sb],
		  30: [30, 2, 1.00, 0.0, 0.0, 0, 10.6/Sb, 1.9/Sb]}

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



# Variables to store the results to later write them in .csv file
nb = len(buses)
nl = len(lineas)

TOTAL_SYSTEM_LOAD = sum(buses[i][6] for i in buses)

DCOPF_SUPLOAD_MEAN = 0
LPACOPF_SUPLOAD_MEAN = 0
N_SUCCESS_DCOPF = 0
N_SUCCESS_LPACOPF = 0
N_SUCCESS_DCPF = 0
N_SUCCESS_LPACPF = 0

# Evaluating the cases
CASES = 10000 # Number of cases to evaluate
N_DAMAGED = 19 # Number of damaged lines

solver = SolverFactory('gurobi', solver_io = 'python')
solverPF = SolverFactory('ipopt')

print("Starting evaluation for {0:.0f} cases".format(CASES))
for k in range(1, CASES+1):
	busesDC = copy.deepcopy(buses)
	busesLPAC = copy.deepcopy(buses)
	lineas_ = copy.deepcopy(lineas)
	gensDC = copy.deepcopy(gens)
	gensLPAC = copy.deepcopy(gens)
	shuntsDC = copy.deepcopy(shunts)
	shuntsLPAC = copy.deepcopy(shunts)
	
	DCOPF_SOLVED = 0
	LPACOPF_SOLVED = 0
	DCPF_SOLVED = 0
	LPACPF_SOLVED = 0
	
	damagedLines = random.sample(lineas_.items(), N_DAMAGED)
	for (i, v) in damagedLines:
		lineas_.pop(i)
		
	nl = len(lineas_)
	
	lineasDC = copy.deepcopy(lineas_)
	lineasLPAC = copy.deepcopy(lineas_)
	
	modelDC = ModelOPF_DC(busesDC, lineasDC, gensDC)
	try:
		result = solver.solve(modelDC)
		if result.solver.termination_condition == TerminationCondition.optimal:
			
			DCOPF_SOLVED = 1		
			N_SUCCESS_DCOPF += 1
			
			LOAD_SUPPLIED_THIS = sum(busesDC[i][6]*modelDC.l[i]() for i in busesDC)
			DCOPF_SUPLOAD_MEAN += LOAD_SUPPLIED_THIS/TOTAL_SYSTEM_LOAD

			# Now we will use this solution as setpoint for AC Power Flow
			busesPF = copy.deepcopy(busesDC)
			lineasPF = copy.deepcopy(lineasDC)
			shuntsPF = copy.deepcopy(shuntsDC)
			
			for i in modelDC.buses:
				if busesPF[i][1] == 0: # slack
					busesPF[i][4] = modelDC.Pg[i]()
				elif busesPF[i][1] == 1: # PV
					busesPF[i][3] = modelDC.th[i]()
				elif busesPF[i][1] == 2: # PQ
					busesPF[i][3] = modelDC.th[i]()
					
				busesPF[i][6] *= modelDC.l[i]()
				busesPF[i][7] *= modelDC.l[i]()
				
			modelPF = BuildPFModel(busesPF, lineasPF, shuntsPF)
			resultPF = solverPF.solve(modelPF)
			if resultPF.solver.termination_condition == TerminationCondition.optimal:
				DCPF_SOLVED = 1
				N_SUCCESS_DCPF += 1
			
			# Vpf, thetapf, Pgpf, Qgpf, Pflowpf, Qflowpf, x, exitflag = ExecutePF(busesPF, lineasPF, shuntsPF)
			# if exitflag == 1: # converged
				# DCPF_SOLVED = 1
				# N_SUCCESS_DCPF += 1
				
				# LOAD_SUPPLIED_THIS = sum(busesPF[i][6] for i in busesPF)
				# DCOPF_SUPLOAD_MEAN += LOAD_SUPPLIED_THIS/TOTAL_SYSTEM_LOAD
				
	except ValueError:
		print('Exception DCOPF in case ' + str(k) + ' for lines ' + str([k for (k,v) in damagedLines]) + '.')
		
		
	modelLPAC = ModelOPF_LPAC(busesLPAC, lineasLPAC, gensLPAC, shuntsLPAC)
	try:
		result = solver.solve(modelLPAC)
		if result.solver.termination_condition == TerminationCondition.optimal:
			
			LPACOPF_SOLVED = 1
			N_SUCCESS_LPACOPF += 1
			
			LOAD_SUPPLIED_THIS = sum(busesLPAC[i][6]*modelLPAC.l[i]() for i in busesLPAC)
			LPACOPF_SUPLOAD_MEAN += LOAD_SUPPLIED_THIS/TOTAL_SYSTEM_LOAD
			
			# Now we will use this solution as setpoint for AC Power Flow
			busesPF = copy.deepcopy(busesLPAC)
			lineasPF = copy.deepcopy(lineasLPAC)
			shuntsPF = copy.deepcopy(shuntsLPAC)
			
			for i in busesLPAC:
				if busesPF[i][1] == 0: # slack
					busesPF[i][2] += modelLPAC.fi[i]()
					busesPF[i][3] = modelLPAC.th[i]()
					busesPF[i][4] = modelLPAC.Pg[i]()
					busesPF[i][5] = modelLPAC.Qg[i]()
				elif busesPF[i][1] == 1: # PV
					busesPF[i][2] += modelLPAC.fi[i]()
					busesPF[i][3] = modelLPAC.th[i]()
					busesPF[i][4] = modelLPAC.Pg[i]()
					busesPF[i][5] = modelLPAC.Qg[i]()
				elif busesPF[i][1] == 2: # PQ
					busesPF[i][2] += modelLPAC.fi[i]()
					busesPF[i][3] = modelLPAC.th[i]()
					busesPF[i][4] = modelLPAC.Pg[i]()
					busesPF[i][5] = modelLPAC.Qg[i]()
					
				busesPF[i][6] *= modelLPAC.l[i]()
				busesPF[i][7] *= modelLPAC.l[i]()
				
			modelPF = BuildPFModel(busesPF, lineasPF, shuntsPF)
			resultPF = solverPF.solve(modelPF)
			if resultPF.solver.termination_condition == TerminationCondition.optimal:
				LPACPF_SOLVED = 1
				N_SUCCESS_LPACPF += 1
				
			# Vpf, thetapf, Pgpf, Qgpf, Pflowpf, Qflowpf, x, exitflag = ExecutePF(busesPF, lineasPF, shuntsPF)
			# if exitflag == 1: # converged
				# LPACPF_SOLVED = 1
				# N_SUCCESS_LPACPF += 1
				# LOAD_SUPPLIED_THIS = sum(busesPF[i][6] for i in busesPF)
				# LPACOPF_SUPLOAD_MEAN += LOAD_SUPPLIED_THIS/TOTAL_SYSTEM_LOAD
			
				
	except ValueError:
		print('Exception LPACOPF in case ' + str(k) + ' for lines ' + str([k for (k,v) in damagedLines]) + '.')
		
	if N_SUCCESS_DCPF != 0 and N_SUCCESS_LPACPF != 0:
		print('[' + str(k) + '/' + str(CASES) + ']: DPL: ' + str([k for (k,v) in damagedLines]) + ' - [OPFDC: ' + str(DCOPF_SOLVED) + ', PFDC: ' + str(DCPF_SOLVED) + ', OPFLPAC: ' + str(LPACOPF_SOLVED) + ', PFLPAC: ' + str(LPACPF_SOLVED) + '] - [DC: {0:.2f}% LPAC: {1:.2f}%].'.format(DCOPF_SUPLOAD_MEAN/N_SUCCESS_DCOPF *100, LPACOPF_SUPLOAD_MEAN/N_SUCCESS_LPACOPF *100))
	


print("			E-DC-LPP					E-LPAC-LPP")
print("		Solved	SolvedPF	% Delivered		Solved	SolvedPF	% Delivered")
print("		{0:.0f}	{1:.0f}		{2:.2f}%			{3:.0f} 	{4:.0f}		{5:.2f}%".format(N_SUCCESS_DCOPF, N_SUCCESS_DCPF, DCOPF_SUPLOAD_MEAN/N_SUCCESS_DCOPF *100, N_SUCCESS_LPACOPF, N_SUCCESS_LPACPF, LPACOPF_SUPLOAD_MEAN/N_SUCCESS_LPACOPF *100))
# print("LPAC:			{0:.2f}%		{1:.2f}%		{2:.2f}%".format(percLPAC, ps_LPAC, pf_LPAC))
# print("AC:			{0:.2f}%		{1:.2f}%".format(percAC, ps_AC))

