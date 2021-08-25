/*===================================================================*\
                                 Biofilm
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef BIOFILM_H
#define BIOFILM_H

#include "Chlorine.h"

class Biofilm : public Chlorine
{
public:

  // determine biofilm distribution
  calculateBiofilm();

protected:

private:

  // additional source terms for simplifying if necessary
  double S1();
  double S2();
  double S3();
  // ...


  // list of parameters
  double p1,p2,p3; // ....

};

#endif