import json
import numpy as np
import os

network_list = []
network_list.append("ujhermes")
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
network_list.append("simasag")
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
#network_list.append("sanchegy")
#network_list.append("ferto")

gamma_orig = []
gamma_orig_std = []
gamma_opt = []
gamma_opt_std = []

for case in network_list:
	print("Case: " + case)
	if(case=="tomalom"):
		file_name = "results_" + case + "_sigmaGammaApprox_Results2.txt"
	else:
		file_name = "results_" + case + "_sigmaGammaApprox_Results.txt"
	
	print(file_name)

	f = open(file_name,'r')
	data = f.read()
	s = data.find("array")+7
	e = data.find("function")-5
	data = data[s:e]
	v = np.fromstring(data, dtype=float, sep=',')
	f.close()

	f = open('input_vis.csv','w')
	for i in range(0,len(v),2):
		f.write(str(v[i]) + ', ' + str(v[i+1]) + '\n')
	f.close()

	call = './isoVis.out ' + case + ' input_vis.csv'
	os.system(call)
	f = open('results_vis.txt')
	data = f.read()
	v = np.fromstring(data, dtype=float, sep=',')
	gamma_orig.append(v[0])
	gamma_orig_std.append(v[1])
	gamma_opt.append(v[2])
	gamma_opt_std.append(v[3])

f = open('gamma_cases.txt','w')
for i in network_list:
	f.write(str(i)+'\n')
f.close()
f = open('gamma_orig.txt','w')
for i in gamma_orig:
	f.write(str(i)+'\n')
f.close()
f = open('gamma_orig_std.txt','w')
for i in gamma_orig_std:
	f.write(str(i)+'\n')
f.close()
f = open('gamma_opt.txt','w')
for i in gamma_opt:
	f.write(str(i)+'\n')
f.close()
f = open('gamma_opt_std.txt','w')
for i in gamma_opt_std:
	f.write(str(i)+'\n')
f.close()
