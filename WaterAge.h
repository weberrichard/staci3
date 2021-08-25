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

#include "QualitySolver.h"

class WaterAge : public QualitySolver
{
public:

  // determine water age, filling up Node waterAge variables, calling ODE45Solver
  calculateWaterAge();

protected:

private:  

};

#endif