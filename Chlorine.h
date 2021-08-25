/*===================================================================*\
                                 Chlorine
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef CHLORINE_H
#define CHLORINE_H

#include "WaterAge.h"

class Chlorine : public WaterAge
{
public:

  // determine chlorine concentration, fill up chlorineConcentration variables in Node
  calculateChlorine();

protected:

private:

  // parameter of chlorine model
  double chlorineDecay;

};

#endif