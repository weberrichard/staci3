import os

case = 'balf.inp'

# getting minimum index of pipes
s = './DMAFF.out ' + case + ' Min'
os.system(s)
f = open("results.txt", "r")
min_index = int(f.read())
print(min_index)

# getting maximum index of pipes
s = './DMAFF.out ' + case + ' Max'
os.system(s)
f = open("results.txt", "r")
max_index = int(f.read())
print(max_index)

# writing indecies to input.csv
f = open("input.csv", "w")
f.write("34\n")
f.write("65\n")
f.write("156\n")
f.close()

s = './DMAFF.out ' + case + ' nDMA input.csv'
os.system(s)
f = open("results.txt", "r")
nDMA = int(f.read())
print(nDMA)

s = './DMAFF.out ' + case + ' demand input.csv'
os.system(s)
f = open("results.txt", "r")
dem = f.read()
dem = dem.split('\n')
print(dem)
print(dem[1])
print(float(dem[1]))
