import math
import os
import numpy as np
from geneticalgorithm import geneticalgorithm as ga

length = ' ./LeakageCalculator Balf SEL'
readed = os.popen(length).read()
readresults = open('results.txt', 'r')
SEL = readresults.read()
print("A hálózat összes csövének hossza: " + SEL + "m")
#print(SEL)
#print(type(SEL))

#Fajlagos veszteség a cikk szerinti számítása a soproni hálózattal
qv1 = 0.48 * float(SEL) / 3600000
print("A Somos Éva-féle cikk alapján a hálózat várható szivárgása: " + str(qv1) + "m3/s")
#print(qv1)
#print(type(qv1))

#A c értékei
varbound = np.array([[0.00000000001,0.00099]]*1)

#Optimalizálandó fv bevitele, of kiadása
def f(X):
	writetosajt = open('sajt.txt','w')
	writetosajt.write(str(X[0]) + " ")
	writetosajt.close()
	results=' ./LeakageCalculator balf sajt.txt Calculate_Loss'
	getBack=os.popen(results).read()
	resultreader = open('results.txt','r')
	qv2 = resultreader.readline()
	print(qv2+"\n")
	if(qv2 == "inf"):
		qv2 = 99999
	of = abs(qv1-float(qv2))
	return of

algorithm_param = {'max_num_iteration': None,
					'population_size': 500,
					'mutation_probability': 0.9,
					'elit_ratio': 0.01,
					'crossover_probability': 0.9,
					'parents_portion': 0.4,
					'crossover_type':'uniform',
					'max_iteration_without_improv': 10000}


model=ga(function=f,dimension=1,variable_type='real',
	variable_boundaries=varbound,
	algorithm_parameters=algorithm_param)

model.run()


#Kiértékelés
convergence = model.report
solution = model.output_dict
g = open('c_Convergence.txt', 'a')
g.write(str(convergence) + "\n")
g.close()
e = open('c_Results.txt', 'a')
e.write(str(solution) + "\n" + str(algorithm_param) + "\n")
e.close()