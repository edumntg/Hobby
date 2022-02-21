from pyomo.environ import *
import numpy as np
from math import pi, inf

def ModelOPF_LPAC(buses, lineas, gens, shunts):

	nb = len(buses)
	nl = len(lineas)
	ng = len(gens)
	ns = len(shunts)
	
	lh = pi/3
	h = 5
	d = 2*lh/(h+1)

	model = ConcreteModel()

	model.buses = Set(initialize=buses.keys())
	model.lineas = Set(initialize=lineas.keys())
	model.gens = Set(initialize=gens.keys())
	model.shunts = Set(initialize=shunts.keys())

	model.Pg = Var(model.buses, initialize = 0)
	model.Qg = Var(model.buses, initialize = 0)
	model.th = Var(model.buses, bounds=(-pi, pi), initialize = 0)
	model.l = Var(model.buses, bounds=(0,1), initialize = 1)
	model.cosnm = Var(model.buses, model.buses, bounds = (0,1), initialize = 0.1)

	model.Pflow = Var(model.buses, model.buses, initialize = 0)
	model.Qflow = Var(model.buses, model.buses, initialize = 0)
	if ns > 0:
		model.Qshunt = Var(model.shunts, initialize = 0)

	model.v = Var(model.buses, bounds = (1,1), initialize = 1)
	model.fi = Var(model.buses, initialize = 0, bounds = (-0.1, 0.1)) # if voltage limits are 0.95 and 1.05

	# Lets build the Ybus

	Ybus = np.zeros((nb,nb), dtype=np.complex128)
	g = np.zeros((nb,nb))
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
		
	# Lets add shunts
	for i in model.shunts:
		bus = shunts[i][0]-1
		B = 1j*shunts[i][1]
		Ybus[bus][bus] += B
			

	G = Ybus.real;
	B = Ybus.imag;
	
	
	def ObjectiveFunc(model):
		return sum(buses[i][6]*model.l[i] for i in model.buses)

	def KirchoffBusesP(model, bus):
		Pik = 0
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
		Qgbus = model.Qg[bus]
		Qshunt = 0
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
				
		if ns > 0:
			key = [key for (key, v) in shunts.items() if v[0] == bus]
			if key:
				Qshunt = model.Qshunt[key[0]]
				
		return Qgbus == buses[bus][7]*model.l[bus] + Qik + Qshunt
										
	def MinGen_P(model, bus):
		keys = [key for (key, v) in gens.items() if v[0] == bus]
		lb = 0
		if keys:
			lb = gens[keys[0]][1]
			
		return model.Pg[bus] >= lb
		
	def MaxGen_P(model, bus):
		keys = [key for (key, v) in gens.items() if v[0] == bus]
		ub = 0
		if keys:
			ub = gens[keys[0]][2]
			
		return model.Pg[bus] <= ub

	def MinGen_Q(model, bus):
		keys = [key for (key, v) in gens.items() if v[0] == bus]
		lb = 0
		if keys:
			lb = gens[keys[0]][3]
			
		return model.Qg[bus] >= lb
		
	def MaxGen_Q(model, bus):
		keys = [key for (key, v) in gens.items() if v[0] == bus]
		ub = 0
		if keys:
			ub = gens[keys[0]][4]
			
		return model.Qg[bus] <= ub

	def PflowEq1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return model.Pflow[i, j] == -G[i-1][j-1] + (G[i-1][j-1]*model.cosnm[i,j] + B[i-1][j-1]*(model.th[i]-model.th[j]))
		# return model.Pflow[i, j] == (-G[i-1][j-1] + g[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(G[i-1][j-1]*cos(model.th[i]-model.th[j]) + B[i-1][j-1]*sin(model.th[i]-model.th[j]))

	def PflowEq2(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		return model.Pflow[i, j] == -G[i-1][j-1] + (G[i-1][j-1]*model.cosnm[i,j] + B[i-1][j-1]*(model.th[i]-model.th[j]))
		# return model.Pflow[i, j] == (-G[i-1][j-1] + g[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(G[i-1][j-1]*cos(model.th[i]-model.th[j]) + B[i-1][j-1]*sin(model.th[i]-model.th[j]))

	def QflowEq1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return model.Qflow[i, j] == B[i-1][j-1] + (-B[i-1][j-1]*model.cosnm[i,j] + G[i-1][j-1]*(model.th[i]-model.th[j])) + B[i-1][j-1]*(model.fi[i]-model.fi[j])
		# return model.Qflow[i, j] == (B[i-1][j-1] - b[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(-B[i-1][j-1]*cos(model.th[i]-model.th[j]) + G[i-1][j-1]*sin(model.th[i]-model.th[j]))

	def QflowEq2(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		return model.Qflow[i, j] == B[i-1][j-1] + (-B[i-1][j-1]*model.cosnm[i,j] + G[i-1][j-1]*(model.th[i]-model.th[j])) + B[i-1][j-1]*(model.fi[i]-model.fi[j])
		# return model.Qflow[i, j] == (B[i-1][j-1] - b[i-1][j-1])*model.v[i]**2 + model.v[i]*model.v[j]*(-B[i-1][j-1]*cos(model.th[i]-model.th[j]) + G[i-1][j-1]*sin(model.th[i]-model.th[j]))

	def QShunt(model, shunt):
		bus = shunts[shunt][0]
		Bs = shunts[shunt][1]
		# Z = np.conj(Z)
		# B = 1/Z;
		# B = B.imag;
		return model.Qshunt[shunt] == Bs*(model.v[bus] + model.fi[bus])*(model.v[bus] + model.fi[bus]) #V^2 /conj(Z)

	def MaxMVAline1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return model.Pflow[i,j]**2 + model.Qflow[i,j]**2 <= lineas[linea][6]**2
		
	def MaxMVAline2(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return model.Pflow[i,j]**2 + model.Qflow[i,j]**2 <= lineas[linea][6]**2
		
	def Cosnm1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		return model.cosnm[i,j] >= cos(lh)
		
	def Cosnm2(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		return model.cosnm[i,j] >= cos(lh)
		
	def Cosnm_tang_1(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		t = 1
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_2(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		t = 2
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_3(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		t = 3
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_4(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		t = 4
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_5(model, linea):
		i = lineas[linea][0]
		j = lineas[linea][1]
		t = 5
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_6(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		t = 1
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_7(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		t = 2
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_8(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		t = 3
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
		
	def Cosnm_tang_9(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		t = 4
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - lh - d) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def Cosnm_tang_10(model, linea):
		i = lineas[linea][1]
		j = lineas[linea][0]
		t = 5
		# return model.cosnm[i,j] <= -sin(lh + d)*(model.th[i] - model.th[j] - d - lh) + cos(d + lh)
		return model.cosnm[i,j] <= -sin(t*d-lh)*(model.th[i] - model.th[j] - t*d + lh) + cos(t*d - lh)
		
	def FiConst(model, bus):
		return model.v[bus] == buses[bus][2] + model.fi[bus] # operational voltage defined as setpoint + fi

	
	model.obj = Objective(rule = ObjectiveFunc, sense=maximize)

	model.c0 = Constraint(expr=model.th[1] == 0) # angle reference, in this case is bus 1
	model.c1 = Constraint(model.buses, rule=KirchoffBusesP)
	model.c2 = Constraint(model.buses, rule=KirchoffBusesQ)

	model.c3 = Constraint(model.buses, rule = MinGen_P)
	model.c4 = Constraint(model.buses, rule = MaxGen_P)
	model.c5 = Constraint(model.buses, rule = MinGen_Q)
	model.c6 = Constraint(model.buses, rule = MaxGen_Q)
		
	model.c7 = Constraint(model.lineas, rule = PflowEq1)
	model.c8 = Constraint(model.lineas, rule = PflowEq2)
	model.c9 = Constraint(model.lineas, rule = QflowEq1)
	model.c10 = Constraint(model.lineas, rule = QflowEq2)

	model.c11 = Constraint(model.lineas, rule = MaxMVAline1)
	model.c12 = Constraint(model.lineas, rule = MaxMVAline2)
	# model.c13 = Constraint(model.lineas, rule = Cosnm1)
	# model.c14 = Constraint(model.lineas, rule = Cosnm2)
	# model.c15 = Constraint(model.lineas, rule = Cosnm3)
	# model.c16 = Constraint(model.lineas, rule = Cosnm4)

	# model.c17 = Constraint(model.buses, rule = FiConst)
	
	model.c19 = Constraint(model.lineas, rule = Cosnm_tang_1)
	model.c20 = Constraint(model.lineas, rule = Cosnm_tang_2)
	model.c21 = Constraint(model.lineas, rule = Cosnm_tang_3)
	model.c22 = Constraint(model.lineas, rule = Cosnm_tang_4)
	model.c23 = Constraint(model.lineas, rule = Cosnm_tang_5)
	model.c24 = Constraint(model.lineas, rule = Cosnm_tang_6)
	model.c25 = Constraint(model.lineas, rule = Cosnm_tang_7)
	model.c26 = Constraint(model.lineas, rule = Cosnm_tang_8)
	model.c27 = Constraint(model.lineas, rule = Cosnm_tang_9)
	model.c28 = Constraint(model.lineas, rule = Cosnm_tang_10)
	
	# Shunts constraints
	if ns > 0:
		model.c18 = Constraint(model.shunts, rule = QShunt)
	
	return model

def PrintOPFLPACResults(model, buses, lineas, gens, shunts):
	print('BusID	V	th	Pg	Qg	l	Pl	Ql	Qshunt\n')
	nb = len(buses)
	nl = len(lineas)
	ng = len(gens)
	ns = len(shunts)
	
	l = {}
	v = {}
	th = {}
	Pg = {}
	Qg = {}
	Qshunt = {}
	for i in model.buses:
		Pg[i] = abs(model.Pg[i]())
		Qg[i] = model.Qg[i]()
		
		l[i] = model.l[i]()
			
		v[i] = model.v[i]()
			
		th[i] = model.th[i]()
		
		keys = [key for (key, v) in shunts.items() if v[0] == i]
		if keys:
			Qshunt[i] = model.Qshunt[keys[0]]()
		else:
			Qshunt[i] = 0
		
		# print(i,v[i],th[i],Pg[i],Qg[i], l[i])
		print("{0:.0f}	{1:.4f}	{2:.4f}	{3:.4f}	{4:.4f}	{5:.4f}	{6:.4f}	{7:.4f} {8:.4f}".format(i,v[i],th[i],Pg[i],Qg[i], l[i], buses[i][6]*l[i], buses[i][7]*l[i], Qshunt[i]))

	Pgtotal = sum(model.Pg[i]() for i in model.buses)
	Qgtotal = sum(model.Qg[i]() for i in model.buses)

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

	Pl_supplied = sum(buses[i][6]*model.l[i]() for i in model.buses)
	Pl_total = sum(buses[i][6] for i in model.buses)
	perc_supplied = (Pl_supplied/Pl_total)*100
	
	print("\n")
	print("Total Ploss: {0:.4f}\nTotal Qloss: {1:.4f}".format(Ploss,Qloss))
	print("Total Load Supplied: {0:.4f}%".format(perc_supplied))
	print("\n")
	print("i	j	fi	cosnm")
	for linea in model.lineas:
		i = lineas[linea][0]
		j = lineas[linea][1]
		print("{0:.0f}	{1:.0f}	{2:.4f} {3:.4f}".format(i, j, model.fi[i](), model.cosnm[i,j]()))

