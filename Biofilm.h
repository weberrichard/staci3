/*===================================================================*\
                                 Biofilm
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef biofilm_H
#define biofilm_H
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

class biofilm
{
public:
  biofilm();
  ~biofilm();
  // determine water age, filling up Node biofilm variables, calling ODE45Solver
  vector<double> sourceTermBiofilm(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX);

  void setParameters(std::string which, double value);

  double Diff = 0.05;
  double Dcs = 1;
  double W = 0.01;
  double mu_f = 2017/3600;
  double mu_b = 7117/3600;
  double Ycs = 0.35;
  double Kf = 30.82;
  double Kb = 265.7;
  double C_fm = 6.775;
  double C_bm = 5737/1000;
  double Yf = 2.607*pow(10,5);
  double Yb = 1.729*pow(10,8);

  int ModelDimension = 4;
protected:

private:  
};

#endif