import os
import numpy as np
import time
import math
import csv
from geneticalgorithm import geneticalgorithm as ga
from os import listdir
import yaml
import ast

Network_List = []
stack = []
with open('Results.csv', newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		if(row[0] == "Network"):
			a_per_b_list = row[1:]
		else:
			Network_List.append(row[0])

multipliers = ["0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1","1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2"]
for i in range(0,len(multipliers)):
	ActualNetwork = "sopronkovesd"
	ActualMultiplier = multipliers[i]
	d = {}
	valvelist = []
	valves = []
	start_it = False
	stack = []
	stack2 = ""
	with open("results_"+str(ActualNetwork)+"_"+str(multipliers[i])+"_Results.txt") as f:
		lines = [line.strip().split('[') for line in f]
		for line in lines:
			d = str(line[:])
			for character in line:
				valvelist.append(character)
				for char in character:
					#print(char)
					#input("Press Enter to continue...")
					if char == "]" or char == ".":
						start_it = False
						try:
							stack.append(int(stack2))
							stack2 = ""
						except:
							print("kicked")
					if start_it == True:
						try:
							stack2 += str(int(char))
						except:
							print(char)
					if char == "[" or  char == " " or char == ",":
						start_it = True
		try:
			stack.pop()
		except:
			print(stack)

		ww = open("input3.txt", "w")
		stringW = ""
		stringR = ""
		i = 0
		for st in stack:
			if(i%2 == 0):
				stringR = str(st)
				i = i + 1
			else:
				stringW = stringR + "," + str(st) + " "+"\n"
				ww.write(stringW)
				i = i + 1
		ww.close()

		callTopology2 = './isoFF.out '+ ActualNetwork +' gamma_distribution input3.txt'
		out = os.popen(callTopology2).read()
		print(callTopology2)

		callTopology3 = './isoFF.out '+ ActualNetwork +' get_Topology input3.txt '+ str(ActualMultiplier)
		out = os.popen(callTopology3).read()
		print(callTopology3)

callTopology = './isoFF.out '+ ActualNetwork + ' GetOriginalGamma'
out = os.popen(callTopology).read()
print(callTopology)