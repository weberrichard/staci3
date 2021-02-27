/*===================================================================*\
                              Sensitivity
                            ---------------

  Derived from HydraulicSolver class. This class is capable of
  calculating the sensitivity to a paramter on a per pipe base.
  Parameter can be Pipe diameter, roughness or Node demand.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/
#ifndef SENSITIVITY_H
#define SENSITIVITY_H

#include "HydraulicSolver.h"

class Sensitivity : public HydraulicSolver
{
public:
  MatrixXd massFlowRateSensitivity; // Sensitivity Matrix
  MatrixXd pressureSensitivity; // Sensitivity Matrix
  double massFlowRateSensitivityEdgeMax;
  double pressureSensitivityNodeMax;

  MatrixXd getPSensMatrix() // It makes reachable the nodal part of the sens. matrix.
  {
    return pressureSensitivity;
  }
  MatrixXd getMFRSensMatrix() // It makes reachable the link section of the sens. matrix.
  {
    return massFlowRateSensitivity;
  }
  double getLargestSensitivity() // It makes reachable the largest NODAL sensitivity
  {
    return SensitivityNodalMax;
  }
  void setLargestSensitivity(double a) // This function modifies the maximal NODAL sensitivity, it should be used only in the cases, when something modifies the sensitivity matrix.
  {
    SensitivityNodalMax = a;
  }

  double CalculateAverageSensitivity(); // It calculates the average value of the NODAL part of the sensitivity matrix, which means all of the rows are summarized and average row sum value is calculated.

	Sensitivity(string spr_filename);
	~Sensitivity();

  /*! Calculating the sensitivity on a per-pipe basis.
      Paramter can be either Pipe diameter, roughness or Node demand.
      diameter | friction_coeff | demand */
	bool calculateSensitivity(string parameter);
private:
  double SensitivityNodalMax = 0.0;  // The maximal value of the NODAL sensitivity.
};

#endif