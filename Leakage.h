/*===================================================================*\
                             Vulnerability
                            ---------------
	
	This class is capable of calculating the vulnerability.
	For def see article XY.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef LEAKAGE_H
#define LEAKAGE_H

#include "HydraulicSolver.h"

class Leakage : public HydraulicSolver
{
public:
	Leakage(string spr_filename);
	~Leakage();

	// it calculates the nodal leakage for the whole network
	void calculateLeakage();


	double getSummarizedLeakage()
	{
		return summarizedLeakage;
	}

	// relative demand loss in case of the loss of the segment
	double summarizedLeakage;
};

#endif