from pyomo.environ import *
import numpy as np
from math import pi, inf

def ModelOPF_DC(buses, lineas, gens):

	nb = len(buses)
	nl = len(lineas)
	ng = len(gens)
	
	model = ConcreteModel()

	model.buses = Set(initialize=buses.keys())
	model.lineas = Set(initialize=lineas.keys())
	model.gens = Set(initialize=gens.keys())

	model.Pg = Var(model.buses, initialize = 0)
	model.th = Var(model.buses, bounds = (-pi, pi), initialize = 0)
	model.l = Var(model.buses, bounds=(0,1), initialize=1)

	# Lets build the Ybus
	Ybus = np.zeros((nb,nb), dtype=np.complex128)
	b = np.zeros((nb,nb))
	# print(Ybus)
	for l in model.lineas:
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
			
	B = Ybus.imag
	
	def ObjectiveFunc(model):
		return sum(buses[i][6]*model.l[i] for i in model.buses)

	def KirchoffBusesP(model, bus):
		Pik = 0
		Pgbus = model.Pg[bus]
			
		for linea in model.lineas:
			i = lineas[linea][0]
			if i == bus: # gen id is the same as bus id
				j = lineas[linea][1]
				Pik += B[i-1][j-1]*(model.th[i]-model.th[j])
				
		for linea in model.lineas:
			i = lineas[linea][1]
			if i == bus: # gen id is the same as bus id
				j = lineas[linea][0]
				Pik += B[i-1][j-1]*(model.th[i]-model.th[j])
			
		return Pgbus == buses[bus][6]*model.l[bus] + Pik
										
	def MinGen_P(model, bus):
		keys = [key for (key, v) in gens.items() if v[0] == bus]
		lb = 0
		if keys:
			lb = gens[keys[0]][1]
			
		return model.Pg[bus] >= 0
		
	def MaxGen_P(model, bus):
		keys = [key for (key, v) in gens.items() if v[0] == bus]
		ub = 0
		if keys:
			ub = gens[keys[0]][2]
			
		return model.Pg[bus] <= ub

	def PflowEq1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return model.Pflow[i, j] == B[i-1][j-1]*(model.th[i]-model.th[j])

	def PflowEq2(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		return model.Pflow[i, j] == B[i-1][j-1]*(model.th[i]-model.th[j])

	def MaxMVAline_ik1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		Pik = B[i-1][j-1]*(model.th[i]-model.th[j])
		return B[i-1][j-1]*(model.th[i]-model.th[j]) <= lineas[linea][6]
		
	def MinMVAline_ik1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return B[i-1][j-1]*(model.th[i]-model.th[j]) >= -lineas[linea][6]
		
	def MaxMVAline_ik2(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		return B[i-1][j-1]*(model.th[i]-model.th[j]) <= lineas[linea][6]
		
	def MinMVAline_ik2(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		return B[i-1][j-1]*(model.th[i]-model.th[j]) >= -lineas[linea][6]
		
	model.obj = Objective(rule = ObjectiveFunc, sense = maximize)
	
	model.slack = Constraint(expr=model.th[1] == 0)
	
	model.c0 = Constraint(model.buses, rule=KirchoffBusesP)

	model.c1 = Constraint(model.buses, rule = MinGen_P)
	model.c2 = Constraint(model.buses, rule = MaxGen_P)
		


	model.c5 = Constraint(model.lineas, rule = MaxMVAline_ik1)
	model.c6 = Constraint(model.lineas, rule = MaxMVAline_ik2)
	model.c7 = Constraint(model.lineas, rule = MinMVAline_ik1)
	model.c8 = Constraint(model.lineas, rule = MinMVAline_ik2)
	
	return model

# print(Ybus)

def PrintOPFDCResults(model, buses, lineas, gens):
	print('BusID	th	Pg	l	Pl\n')
	nb = len(buses)
	nl = len(lineas)
	ng = len(gens)
	
	l = {}
	th = {}
	Pg = {}
	for i in model.buses:
		Pg[i] = model.Pg[i]()
		
		l[i] = model.l[i]()
		if l[i] is None:
			l[i] = 0
			
		th[i] = model.th[i]()
		if th[i] is None:
			th[i] = 0
		
		# print(i,v[i],th[i],Pg[i],Qg[i], l[i])
		print("{0:.0f}	{1:.4f}	{2:.4f}	{3:.4f} {4:.4f}".format(i,th[i],Pg[i], l[i], buses[i][6]*l[i]))

	Pgtotal = sum(model.Pg[i]() for i in model.gens)

	Ploadtotal = sum(buses[i][6]*l[i] for i in model.buses)

	print("\n")
	print("TOTAL			{0:.4f}	{1:.4f}".format(Pgtotal, Ploadtotal))
	print("\n\n")

	print("Busi	Busk	Pik	Pki")
	Pik = np.zeros((nb,nb))
	Pki = np.zeros((nb,nb))

	for l in model.lineas:
		i = lineas[l][0]
		j = lineas[l][1]
		print("{0:.0f}	{1:.0f}	{2:.4f}	{3:.4f}".format(i,j,model.Pflow[i,j](),model.Pflow[j,i]()))
		
	Ploss = 0
	for l in model.lineas:
		i = lineas[l][0]
		j = lineas[l][1]
		Ploss += model.Pflow[i,j]() + model.Pflow[j,i]()


	Pl_supplied = sum(buses[i][6]*model.l[i]() for i in model.buses)
	Pl_total = sum(buses[i][6] for i in model.buses)
	perc_supplied = (Pl_supplied/Pl_total)*100
	
	print("\n")
	print("Total Ploss: {0:.4f}".format(Ploss))
	print("Total Load Supplied: {0:.4f}%".format(perc_supplied))