/*===================================================================*\
                                 chlorine
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef chlorine_H
#define chlorine_H
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

class chlorine
{
public:
  chlorine();
  ~chlorine();
  // determine water age, filling up Node chlorine variables, calling ODE45Solver
  vector<double> sourceTermChlorine(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX);

  void setParameters(std::string which, double value);

  double kb = 0.0006417; // Chlorine bulk decay coefficient, based on the work of L. Monteiro et al. (Modeling of chlorine decay in drinking water supply systems using 
                            // EPANET MSX)
  double n = 1; // The order of the decay, n-th order -> 1.2 (kb must be modified)

  int ModelDimension = 1;
protected:

private:  
};

#endif