/*===================================================================*\
                              QualitySolver
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef QUALITYSOLVER_H
#define QUALITYSOLVER_H

#include "SeriesHydraulics.h"

class QualitySolver : public SeriesHydraulics
{
public:

  // setting adaptive parameters for dt handling
  setAdaptiveParameters(double e_max, double e_min, double dt_inc, double dt_dec);
  // setting spatial parameters
  setSpatialParameters(double dx_min, double dx_max);

protected:

  // general ode45 solver with S at the right side, i.e. dx/dt = S
  ODE45Solver(vector<double> S);

private:

  // for adaptive timestep handling
  double error_max;
  double error_min;
  double dt_increase;
  double dt_decrease; 

  // actual timestep
  double dt;

  // spacial division
  double dx_min, dx_max;

  // general x in ode
  vector<double> x;

};

#endif