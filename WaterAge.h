/*===================================================================*\
                                 WaterAge
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef WATERAGE_H
#define WATERAGE_H
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;



class WaterAge
{
public:
  WaterAge();
  ~WaterAge();
  int ModelDimension = 1;
  vector<string> listOfParameters = {"waterAge"};
  vector<double> sourceTermWaterAge(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX);
  vector<double> nodeEquationWaterAge(vector< vector<double> > nodeInputs, vector<double> VolFlowRates);
  vector<double> poolEquationWaterAge(double poolFlow, vector<double> poolActual, vector<double> poolUpstreamNode, double poolVolumeActual, double h);
protected:

private:  

};

#endif