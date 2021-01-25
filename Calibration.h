/*===================================================================*\
                              Calibration
                            ---------------

  Main class for calibration techniques, derived from Sensitivity.
  Contains the loading of measurement data and generating artificial
  measurment data for debugging and comparison.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef CALIBRATION_H
#define CALIBRATION_H

#include "Sensitivity.h"
#include "SeriesHydraulics.h"

class Calibration : public SeriesHydraulics
{
public:
	Calibration(string fileName);
	~Calibration();

  void loadMeasurement();
  void savePipes();
  
//protected:
	MatrixXd measuredPressure; // Pressure values at actual measure nodes
  vector<double> measuredTime; // Time appointment in which measurement occured
  double measuredTimeConvert;

  vector<string> measuredPressureID; // Node IDs for pressure
	vector<int> measuredPressureIndex; // Node indicies for pressure

	int numberPressure;  // Number of measured pressure nodes
	int numberOperationalPoint; // Number of measured operational points
};

#endif
