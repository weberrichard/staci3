from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.problems import get_problem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
import numpy as np
from pymoo.core.problem import Problem
import os

case_name = "balf_mat_year_simp"

d_min = 50
d_max = 150

os.system('./relative_shortfall.out ' + case_name + ' pn')
f = open("pn.txt", "r")
pn = f.read()
pn = int(pn)

class OF1(Problem):

    def __init__(self):
        lower = d_min*np.ones(pn+1)
        lower[0] = 10
        upper = d_max*np.ones(pn+1)
        upper[0] = 50
        super().__init__(n_var=pn+1, n_obj=2, n_ieq_constr=0, xl=lower, xu=upper)

    def _evaluate(self, x, out, *args, **kwargs):
        n = len(x)
        of = np.empty(shape=(n,2))
        of.fill(0)
        for i in range(n):
            y = x[i,:]
            np.savetxt('diameter.txt', y[1:])
            run_str = './relative_shortfall.out ' + case_name + ' of_p1 ' + str(y[0])
            os.system(run_str)
            f = open("of_p1.txt", "r")
            f = f.read()
            of[i,:] = np.fromstring(f, dtype=float, sep='\n')

        print(of)
        out["F"] = of


algorithm = NSGA2(pop_size=100)

problem = OF1()

res = minimize(problem,algorithm,('n_gen', 200),seed=1,verbose=False)

print(res)

plot = Scatter()
plot.add(res.F, facecolor="none", edgecolor="red")
plot.show()