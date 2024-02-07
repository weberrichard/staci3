import numpy as np
from geneticalgorithm import geneticalgorithm as ga
import os

of_best = 1e10
diam_best = ''

def of(x):
	np.savetxt('diameter.txt', x)
	os.system('./relative_shortfall.out ' + case_name + ' of')
	f = open("of.txt", "r")
	of = f.read()
	of = np.array(of)
	of = of.astype(float)

	global of_best

	if of_best>of:
		of_best
		of_best = of
		diam_best = x
		print(of_best)
		np.savetxt('diameter_best.txt', diam_best)


	return of

#np.savetxt('diameter.txt', X, delimiter=',')

case_name = "balf_mat_year_simp"

os.system('./relative_shortfall.out ' + case_name + ' pn')
f = open("pn.txt", "r")
pn = f.read()
pn = int(pn)

varbound=np.ones((pn,2))

varbound[:,0] = 50
varbound[:,1] = 125

model=ga(function=of,dimension=pn,variable_type='real',variable_boundaries=varbound)
model.run()


'''
varbound=np.array([[-10,10]]*3)
print(varbound)

model=ga(function=f,dimension=3,variable_type='real',variable_boundaries=varbound)

model.run()
'''

# possible networks
# "villasor_mat_year", "ferto_mat_year", "sanchegy_mat_year", "buk_mat_year", "lovo_mat_year", "nagycenk_mat_year", "vashegy_mat_year", "varis_mat_year", "becsidomb_mat_year", "tomalom_mat_year", "szakov_mat_year", "kohegy_mat_year", "harka_mat_year", "pozsonyiut_mat_year", "sopronkovesd_mat_year", "dudlesz_mat_year", "ivan_mat_year", "agyagosszergeny_mat_year", "kofejto_mat_year", "simasag_mat_year", "acsad_mat_year", "csaford_mat_year", "nagylozs_mat_year", "balf_mat_year", "csapod_mat_year", "und_mat_year", "rojtokmuzsaj_mat_year", "brennberg_mat_year", "pusztacsalad_mat_year", "kutyahegy_mat_year", "nyarliget_mat_year", "meszlen_mat_year", "fertoujlak_mat_year", "gorbehalom_mat_year", "tozeggyarmajor_mat_year", "ebergoc_mat_year", "csillahegy_mat_year", "jerevan_mat_year", "gloriette_mat_year", "alomhegy_mat_year", "ohermes_mat_year", "ujhermes_mat_year"
