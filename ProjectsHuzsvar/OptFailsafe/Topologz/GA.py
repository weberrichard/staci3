import os
import numpy as np
import time
import math
import csv
from geneticalgorithm import geneticalgorithm as ga
from os import listdir

#-----------#
aa = 0
bb = 0
Network_out = ""
#-----------#

def func(X):
	#print(X)
	a = Parameters_aa()
	b = Parameters_bb()
	f = open("input.csv", "w")
	Network = Parameters_Network()
	for i in range(len(X)):
		#print(i)
		if(i%2 == 0 and i <len(X)-1):
			f.write(str(X[i])+","+str(X[i+1])+"\n")
			#print(str(X[i])+","+str(X[i+1])+"\n")
	f.close()
	cmd = './isoFF.out ' + Network + ' Gamma input.csv'
	#print(cmd)
	so = os.popen(cmd).read()
	ff = open("results.txt", "r")
	so = ff.read()
	#print(so)
	return float(so)

def Parameters_I(Network, a, b):
	global aa
	global bb
	global Network_out
	aa = a
	bb = b
	Network_out = Network

def Parameters_aa():
	global aa
	a_ret = aa
	return a_ret

def Parameters_bb():
	global bb
	b_ret = bb
	return b_ret

def Parameters_Network():
	global Network_out
	Net_re = Network_out
	return Net_re

multipliers = ["0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1","1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2"]
Network_List = []
with open('Results.csv', newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		if(row[0] == "Network"):
			a_per_b_list = row[1:]
		else:
			Network_List.append(row[0])

for i in range(32,len(Network_List)):
	ActualNetwork = Network_List[i]
	for j in range(len(multipliers)):
		a = 1
		b = 1

		callTopology = './isoFF.out '+str(ActualNetwork)+' Min'
		Limit_Min = os.popen(callTopology).read()
		f = open("results.txt", "r")
		Limit_Min = f.read()
		f.close()

		callTopology = './isoFF.out '+str(ActualNetwork)+' Max'
		Limit_Max = os.popen(callTopology).read()
		f = open("results.txt", "r")
		Limit_Max = f.read()
		f.close()
		
		callTopology = './isoFF.out '+str(ActualNetwork)+' Valve'
		Limit_Valve = os.popen(callTopology).read()
		f = open("results.txt", "r")
		Limit_Valve = f.read()
		f.close()

		print(Limit_Min," ", Limit_Max, " ", Limit_Valve)
		#print(Limit_Min,Limit_Max,Limit_Valve)
		varbound=np.array([[int(Limit_Min),int(Limit_Max)-1],[0,1]]*(round(int(Limit_Valve)*float(multipliers[j]))))
		vartype=np.array([['int'],['int']]*round(int(Limit_Valve)*float(multipliers[j])))

		Parameters_I(ActualNetwork, a, b)
		algorithm_param = {'max_num_iteration': None,\
		                   'population_size':50,\
		                   'mutation_probability':0.01,\
		                   'elit_ratio': 0.05,\
		                   'crossover_probability': 0.8,\
		                   'parents_portion': 0.1,\
		                   'crossover_type':'uniform',\
		                   'max_iteration_without_improv':100}

		model=ga(function=func,\
		            dimension=round(int(Limit_Valve)*float(multipliers[j]))*2,\
		            variable_type_mixed=vartype,\
		            variable_boundaries=varbound,\
		            algorithm_parameters=algorithm_param)

		model.run()
		convergence=model.report
		solution=model.output_dict
		
		print(convergence)
		print(solution)
		g = open("results_"+str(ActualNetwork)+"_"+str(multipliers[j])+"_Convergence.txt", "w")
		g.write(str(convergence))
		g.close
		e = open("results_"+str(ActualNetwork)+"_"+str(multipliers[j])+"_Results.txt", "w")
		e.write(str(solution))
		e.close
