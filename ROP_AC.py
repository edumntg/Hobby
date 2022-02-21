from pyomo.environ import *
import numpy as np
from math import pi
import random

# solver = SolverFactory("ipopt")

Sb = 100 # MVAbase 10MW

#		ID	Type Volt	Angle	Pg		Qg		Pl		Ql

buses = { 1: [1, 0, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
		  2: [2, 1, 1.01, 0.0, 0.0, 0, 21.7/Sb, 12.7/Sb],
		  3: [3, 2, 1.00, 0.0, 0.0, 0, 2.4/Sb, 1.2/Sb],
		  4: [4, 2, 1.00, 0.0, 0.0, 0, 7.6/Sb, 1.6/Sb],
		  5: [5, 2, 1.00, 0.0, 0.0, 0, 94.2/Sb, 19.0/Sb],
		  6: [6, 2, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
		  7: [7, 2, 1.00, 0.0, 0.0, 0, 22.8/Sb, 10.9/Sb],
		  8: [8, 2, 1.00, 0.0, 0.0, 0, 30.0/Sb, 30.0/Sb],
		  9: [9, 2, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
		  10: [10, 2, 1.00, 0.0, 0.0, 0, 5.8/Sb, 2.0/Sb],
		  11: [11, 2, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
		  12: [12, 2, 1.00, 0.0, 0.0, 0, 11.2/Sb, 7.5/Sb],
		  13: [13, 2, 1.00, 0.0, 0.0, 0, 0.0, 0.0],
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

lineas = { 1: [1, 2, 0.0192, 0.0575, 0.0264, 1, 130/Sb, 1],
		   2: [1, 3, 0.0452, 0.1852, 0.0204, 1, 130/Sb, 1],
		   3: [2, 4, 0.0570, 0.1737, 0.0184, 1, 65/Sb, 1],
		   4: [3, 4, 0.0132, 0.0379, 0.0042, 1, 130/Sb, 1],
		   5: [2, 5, 0.0472, 0.1983, 0.0209, 1, 130/Sb, 1],
		   6: [2, 6, 0.0581, 0.1763, 0.0187, 1, 65/Sb, 1],
		   7: [4, 6, 0.0119, 0.0414, 0.0045, 1, 90/Sb, 1],
		   8: [5, 7, 0.0460, 0.1160, 0.0102, 1, 70/Sb, 1],
		   9: [6, 7, 0.0267, 0.0820, 0.0085, 1, 130/Sb, 1],
		   10: [6, 8, 0.0120, 0.0420, 0.0045, 1, 32/Sb, 1],
		   11: [6, 9, 0.0000, 0.2080, 0.0000, 1.0155, 65/Sb, 1],
		   12: [6, 10, 0.0000, 0.5560, 0.0000, 0.9629, 32/Sb, 1],
		   13: [9, 11, 0.0000, 0.2080, 0.0000, 1, 65/Sb, 1],
		   14: [9, 10, 0.0000, 0.1100, 0.0000, 1, 65/Sb, 1],
		   15: [4, 12, 0.0000, 0.2560, 0.0000, 1.0129, 65/Sb, 1],
		   16: [12, 13, 0.0000, 0.1400, 0.0000, 1, 65/Sb, 1],
		   17: [12, 14, 0.1231, 0.2559, 0.0000, 1, 32/Sb, 1],
		   18: [12, 15, 0.0662, 0.1304, 0.0000, 1, 32/Sb, 1],
		   19: [12, 16, 0.0945, 0.1987, 0.0000, 1, 32/Sb, 1],
		   20: [14, 15, 0.2210, 0.1997, 0.0000, 1, 16/Sb, 1],
		   21: [16, 17, 0.0824, 0.1932, 0.0000, 1, 16/Sb, 1],
		   22: [15, 18, 0.1070, 0.2185, 0.0000, 1, 16/Sb, 1],
		   23: [18, 19, 0.0639, 0.1292, 0.0000, 1, 16/Sb, 1],
		   24: [19, 20, 0.0340, 0.0680, 0.0000, 1, 32/Sb, 1],
		   25: [10, 20, 0.0936, 0.2090, 0.0000, 1, 32/Sb, 1],
		   26: [10, 17, 0.0324, 0.0845, 0.0000, 1, 32/Sb, 1],
		   27: [10, 21, 0.0348, 0.0749, 0.0000, 1, 32/Sb, 1],
		   28: [10, 22, 0.0727, 0.1499, 0.0000, 1, 32/Sb, 1],
		   29: [21, 22, 0.0116, 0.0236, 0.0000, 1, 32/Sb, 1],
		   30: [15, 23, 0.1000, 0.2020, 0.0000, 1, 16/Sb, 1],
		   31: [22, 24, 0.1150, 0.1790, 0.0000, 1, 16/Sb, 1],
		   32: [23, 24, 0.1320, 0.2700, 0.0000, 1, 16/Sb, 1],
		   33: [24, 25, 0.1885, 0.3292, 0.0000, 1, 16/Sb, 1],
		   34: [25, 26, 0.2544, 0.3800, 0.0000, 1, 16/Sb, 1],
		   35: [25, 27, 0.1093, 0.2087, 0.0000, 1, 16/Sb, 1],
		   36: [28, 27, 0.0000, 0.3690, 0.0000, 0.9581, 65/Sb, 1],
		   37: [27, 29, 0.2198, 0.4153, 0.0000, 1, 16/Sb, 1],
		   38: [27, 30, 0.3202, 0.6027, 0.0000, 1, 16/Sb, 1],
		   39: [29, 30, 0.2399, 0.4533, 0.0000, 1, 16/Sb, 1],
		   40: [8, 28, 0.0636, 0.2000, 0.0214, 1, 32/Sb, 1],
		   41: [6, 28, 0.0169, 0.0599, 0.0065, 1, 32/Sb, 1]}
	
R = {} # Lineas a considerar dañadas
	
shunts = {}
		   
gens = {1: [50/Sb, 200/Sb, -0/Sb, 0/Sb, 0.00375, 2.0, 0],
		2: [20/Sb, 80/Sb, -20/Sb, 100/Sb, 0.0175, 1.75, 0],
		3: [15/Sb, 50/Sb, -15/Sb, 80/Sb, 0.0625, 1.0, 0],
		4: [10/Sb, 35/Sb, -15/Sb, 60/Sb, 0.00834, 3.25, 0],
		5: [10/Sb, 30/Sb, -10/Sb, 50/Sb, 0.0250, 3.0, 0],
		6: [12/Sb, 40/Sb, -15/Sb, 60/Sb, 0.0250, 3.0, 0]}
		


def ROPModel(buses, lineas, shunts, gens):

	nb = len(buses)
	nl = len(lineas)
	ng = len(gens)
	ns = len(shunts)

	# Definition of model
	solver = SolverFactory('ipopt') # Solver for nin-linear equations
	model = ConcreteModel()

	model.buses = Set(initialize=buses.keys())
	model.lineas = Set(initialize=lineas.keys())
	model.gens = Set(initialize=gens.keys())
	model.shunts = Set(initialize=shunts.keys())
	model.z = Var(model.lineas, within = Binary)
	model.o = Var(model.lineas, within = Binary)

	model.Pg = Var(model.gens, within = NonNegativeReals)
	model.Qg = Var(model.gens)
	model.v = Var(model.buses, bounds=(0.9,1.1),within = NonNegativeReals)
	model.th = Var(model.buses, bounds=(-pi/4, pi/4))
	model.l = Var(model.buses, bounds=(0,1))

	model.Pl = Param(model.buses)
	model.Ql = Param(model.buses)

	model.Pflow = Var(model.buses, model.buses)
	model.Qflow = Var(model.buses, model.buses)
	model.Qshunt = Var(model.shunts)


	fromDmg = 3
	toDmg = 20
	# Variables will be stored as Var[k][n] where k is the step and n is the var_id
	Vk = np.zeros((toDmg-fromDmg,nb))
	thk = np.zeros((toDmg-fromDmg, nb))
	Pgk = np.zeros((toDmg-fromDmg,ng))
	Qgk = np.zeros((toDmg-fromDmg,ng))
	if ns > 0:
		Qshunt = np.zeros((toDmg-fromDmg,ns))

	R = {}
	for ndmg in range(fromDmg, toDmg): # from 3 to 20 damaged lines
		# LET'S SELECT THE RANDOM ndmg DAMAGED LINES
		R_ = random.sample(set(range(0, nl-1)), ndmg) # get the random dmg lines
		for k in range(0, ndmg-1):
			R[k] = R_
			#Now set these lines as damaged
			for i in model.lines
				lineas[i][7] = 0 # 1 for active, 0 for inactive

			# Lets build the Ybus
			Ybus, g, b = BuildYbus(lineas)
			G = Ybus.real;
			B = Ybus.imag;
		
		
		
	return V, th, Pg, Qg, l, o, z, Pik, Pki, Qik, Qki, 


# Ybus builder
def BuildYbus(lineas):
	Ybus = np.zeros((nb,nb), dtype=np.complex128)
	g = np.zeros((nb,nb))
	b = np.zeros((nb,nb))
	# print(Ybus)
	for l in model.lineas:
		if lineas[l][7] == 1: # line is active
			i = int(lineas[l][0])-1
			j = int(lineas[l][1])-1
			R = lineas[l][2]
			X = lineas[l][3]
			B = 1j*lineas[l][4]
			a = lineas[l][5] #tap
			Z = R + 1j*X
			Ybus[i][i] += (1/Z)/(a**2)
			if (i != j):
				Ybus[j][j] += (1/Z)/(a**2)
				Ybus[i][j] -= (1/Z)/a
				Ybus[j][i] -= (1/Z)/a
				
			Ybus[i][i] += B
			if (i != j):
				Ybus[j][j] += B
				b[i][j] = B.imag
				b[j][i] = B.imag

	# Lets add shunts
	for i in model.shunts:
		bus = shunts[i][0]
		B = 1j*shunts[i][1]
		Ybus[bus][bus] += B

	return Ybus, g, b
	
def ObjectiveFunc(model):
	return sum(buses[i][6]*model.l[i] for i in model.buses)

model.obj = Objective(rule = ObjectiveFunc, sense=maximize)

def SlackBusAngle(model): # define the slack bus angle set to zero
	return model.th[1] == 0

def KirchoffBusesP(model, bus):
	Pik = 0
	Pgbus = 0
	if bus > ng: # this bus doesnt have a generator
		Pgbus = 0
	else:
		Pgbus = model.Pg[bus]
		
	for linea in model.lineas:
		i = lineas[linea][0]
		if i == bus: # gen id is the same as bus id
			j = lineas[linea][1]
			Pik += model.Pflow[i,j]
			
	for linea in model.lineas:
		i = lineas[linea][1]
		if i == bus: # gen id is the same as bus id
			j = lineas[linea][0]
			Pik += model.Pflow[i,j]
		
	return Pgbus == buses[bus][6]*model.l[bus] + Pik

def KirchoffBusesQ(model, bus):
	Qik = 0
	Qgbus = 0
	Qshunt = 0
	if bus > ng: # this bus doesnt have a generator
		Qgbus = 0
	else:
		Qgbus = model.Qg[bus]
		
	for linea in model.lineas:
		i = lineas[linea][0]
		if i == bus: # gen id is the same as bus id
			j = lineas[linea][1]
			Qik += model.Qflow[i,j]
			
	for linea in model.lineas:
		i = lineas[linea][1]
		if i == bus: # gen id is the same as bus id
			j = lineas[linea][0]
			Qik += model.Qflow[i,j]
			
	key = [key for (key, v) in shunts.items() if v[0] == bus]
	if key:
		Qshunt = model.Qshunt[key[0]]
		
	return Qgbus == buses[bus][7]*model.l[bus] + Qik + Qshunt
										
def MinGen_P(model, gen):
	return model.Pg[gen] >= gens[gen][0]
	
def MaxGen_P(model, gen):
	return model.Pg[gen] <= gens[gen][1]

def MinGen_Q(model, gen):
	return model.Qg[gen] >= gens[gen][2]
	
def MaxGen_Q(model, gen):
	return model.Qg[gen] <= gens[gen][3]

def PflowEq1(model, linea):
	i = lineas[linea][0]
	j = lineas[linea][1]
	return model.Pflow[i, j] == (-G[i-1][j-1] + g[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(G[i-1][j-1]*cos(model.th[i]-model.th[j]) + B[i-1][j-1]*sin(model.th[i]-model.th[j]))

def PflowEq2(model, linea):
	i = lineas[linea][1]
	j = lineas[linea][0]
	return model.Pflow[i, j] == (-G[i-1][j-1] + g[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(G[i-1][j-1]*cos(model.th[i]-model.th[j]) + B[i-1][j-1]*sin(model.th[i]-model.th[j]))

def QflowEq1(model, linea):
	i = lineas[linea][0]
	j = lineas[linea][1]
	return model.Qflow[i, j] == (B[i-1][j-1] - b[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(-B[i-1][j-1]*cos(model.th[i]-model.th[j]) + G[i-1][j-1]*sin(model.th[i]-model.th[j]))

def QflowEq2(model, linea):
	i = lineas[linea][1]
	j = lineas[linea][0]
	return model.Qflow[i, j] == (B[i-1][j-1] - b[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(-B[i-1][j-1]*cos(model.th[i]-model.th[j]) + G[i-1][j-1]*sin(model.th[i]-model.th[j]))

def QShunt(model, shunt):
	bus = shunts[shunt][0]
	Z = 1j*shunts[shunt][1]
	Z = np.conj(Z)
	B = 1/Z;
	B = B.imag;
	return model.Qshunt[shunt] == (model.v[bus]**2)*B #V^2 /conj(Z)

def MaxMVAline(model, linea):
	i = lineas[linea][0]
	j = lineas[linea][1]
	return model.Pflow[i,j]**2 + model.Qflow[i,j]**2 <= lineas[linea][6]**2

model.c0 = Constraint(rule = SlackBusAngle)
model.c1 = Constraint(model.buses, rule=KirchoffBusesP)
model.c2 = Constraint(model.buses, rule=KirchoffBusesQ)

model.c3 = Constraint(model.gens, rule = MinGen_P)
model.c4 = Constraint(model.gens, rule = MaxGen_P)
model.c5 = Constraint(model.gens, rule = MinGen_Q)
model.c6 = Constraint(model.gens, rule = MaxGen_Q)
	
model.c7 = Constraint(model.lineas, rule = PflowEq1)
model.c8 = Constraint(model.lineas, rule = PflowEq2)
model.c9 = Constraint(model.lineas, rule = QflowEq1)
model.c10 = Constraint(model.lineas, rule = QflowEq2)

model.c11 = Constraint(model.lineas, rule = MaxMVAline)

# Shunts constraints
if ns > 0:
	model.c12 = Constraint(model.shunts, rule = QShunt)

results = solver.solve(model)

print('BusID	V	th	Pg	Qg	l	Pl	Ql	Qshunt\n')
l = {}
for i in model.buses:
	Pg = 0
	Qg = 0
	Qshunt = 0
	if i < ng:
		Pg = model.Pg[i]()
		Qg = model.Qg[i]()
	
	l[i] = model.l[i]()
	if l[i] is None:
		l[i] = 0
	
	keys = [key for (key, v) in shunts.items() if v[0] == i]
	if keys:
		Qshunt = model.Qshunt[keys[0]]()
		
	print("{0:.0f}	{1:.4f}	{2:.4f}	{3:.4f}	{4:.4f}	{5:.4f}	{6:.4f}	{7:.4f}	{8:.4f}".format(i,model.v[i](),model.th[i](),Pg,Qg, l[i], buses[i][6]*l[i], buses[i][7]*l[i], Qshunt))

Pgtotal = sum(model.Pg[i]() for i in model.gens)
Qgtotal = sum(model.Qg[i]() for i in model.gens)

Ploadtotal = sum(buses[i][6]*l[i] for i in model.buses)
Qloadtotal = sum(buses[i][7]*l[i] for i in model.buses)

print("\n")
print("TOTAL			{0:.4f}	{1:.4f}		{2:.4f}	{3:.4f}".format(Pgtotal, Qgtotal, Ploadtotal, Qloadtotal))
print("\n\n")

print("Busi	Busk	Pik	Pki	Qik	Qki")
Pik = np.zeros((nb,nb))
Pki = np.zeros((nb,nb))
Qik = np.zeros((nb,nb))
Qki = np.zeros((nb,nb))

for l in model.lineas:
	i = lineas[l][0]
	j = lineas[l][1]
	print("{0:.0f}	{1:.0f}	{2:.4f}	{3:.4f}	{4:.4f}	{5:.4f}".format(i,j,model.Pflow[i,j](),model.Pflow[j,i](),model.Qflow[i,j](),model.Qflow[j,i]()))
	
Ploss = 0
Qloss = 0
for l in model.lineas:
	i = lineas[l][0]
	j = lineas[l][1]
	Ploss += model.Pflow[i,j]() + model.Pflow[j,i]()
	Qloss += model.Qflow[i,j]() + model.Qflow[j,i]()

print("\n")
print("Total Ploss: {0:.4f}\nTotal Qloss: {1:.4f}".format(Ploss,Qloss))

print(results)