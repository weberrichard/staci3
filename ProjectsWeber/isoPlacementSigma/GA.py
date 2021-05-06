import os
import numpy as np
import time
import math
import csv
from geneticalgorithm import geneticalgorithm as ga
from os import listdir

case = ''

def func(X):
	network = case
	f = open("input.csv", "w")
	for i in range(len(X)):
		#print(i)
		if(i%2 == 0 and i <len(X)-1):
			f.write(str(X[i])+","+str(X[i+1])+"\n")
			#print(str(X[i])+","+str(X[i+1])+"\n")
	f.close()
	cmd = './isoFF.out ' + network + ' Sigma_gamma_approx input.csv'
	os.system(cmd)
	#print(cmd)
	ff = open("results.txt", "r")
	so = ff.read()
	#print(so)
	return float(so)

network_list = []
'''network_list.append("ujhermes")
network_list.append("ohermes")
network_list.append("alomhegy")
network_list.append("gloriette")
network_list.append("jerevan")
network_list.append("csillahegy")
network_list.append("ebergoc")
network_list.append("tozeggyarmajor")
network_list.append("gorbehalom")
network_list.append("fertoujlak")
network_list.append("meszlen")
network_list.append("nyarliget")
network_list.append("kutyahegy")
network_list.append("pusztacsalad")
network_list.append("brennberg")
network_list.append("rojtokmuzsaj")
network_list.append("und")
network_list.append("csapod")
network_list.append("balf")
network_list.append("nagylozs")
network_list.append("csaford")
network_list.append("acsad")
network_list.append("simasag")'''
network_list.append("kofejto")
network_list.append("agyagosszergeny")
network_list.append("ivan")
network_list.append("dudlesz")
network_list.append("sopronkovesd")
network_list.append("pozsonyiut")
network_list.append("harka")
network_list.append("kohegy")
network_list.append("szakov")
network_list.append("tomalom")
network_list.append("becsidomb")
network_list.append("varis")
network_list.append("vashegy")
network_list.append("nagycenk")
network_list.append("lovo")
network_list.append("sanchegy")
network_list.append("ferto")

for c in network_list:
	case = c
	print('Network: ' + case)
	call = './isoFF.out ' + case + ' Min'
	os.system(call)
	f = open("results.txt", "r")
	limit_min = f.read()
	f.close()

	call = './isoFF.out ' + case + ' Max'
	os.system(call)
	f = open("results.txt", "r")
	limit_max = f.read()
	f.close()

	call = './isoFF.out ' + case + ' Valve'
	os.system(call)
	f = open("results.txt", "r")
	limit_valve = f.read()
	f.close()

	varbound=np.array([[int(limit_min),int(limit_max)-int(limit_min)],[0,1]]*int(limit_valve))
	vartype=np.array([['int'],['int']]*int(limit_valve))
	#varbound=np.array([[int(limit_min),int(limit_max)-int(limit_min)]]*int(limit_valve))
	#vartype=np.array(['int']*int(limit_valve))

	algorithm_param = {'max_num_iteration': None,\
	                   'population_size':50,\
	                   'mutation_probability':0.01,\
	                   'elit_ratio': 0.05,\
	                   'crossover_probability': 0.8,\
	                   'parents_portion': 0.1,\
	                   'crossover_type':'uniform',\
	                   'max_iteration_without_improv':500}

	model=ga(function=func,\
            dimension=int(limit_valve)*2,\
            variable_type_mixed=vartype,\
            variable_boundaries=varbound,\
            algorithm_parameters=algorithm_param)

	model.run()
	convergence=model.report
	solution=model.output_dict
	#print(convergence)
	#print(solution)
	g = open("results_" + case + "_sigmaGammaApprox_Convergence.txt", "w")
	g.write(str(convergence))
	g.close
	e = open("results_" + case + "_sigmaGammaApprox_Results.txt", "w")
	e.write(str(solution))
	e.close
