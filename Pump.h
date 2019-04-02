/*===================================================================*\
                                  Pump
                            ---------------

  Derived from Edge class and can be instantiated. Pump element.
  It gives pressure increment based on characteristic curve.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef PUMP_H
#define PUMP_H

#include "Edge.h"

class Pump: public Edge {

public:
	Pump(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_referenceCrossSection, vector<double> a_volumeFlowRate, vector<double> a_head, const double a_massFlowRate);
	~Pump();

  /// Provides basic informations
	string info();

  /// A line of F(x) = equation, rearranged to 0 in w.c.m.
	double function(vector<double>);

  /// Jacobian: df/dhe, df/dhv, df/dmp
	vector<double> functionDerivative(vector<double>);

  /// Initialization, mode: 0->automatic | 1-> using value
	void initialization(int mode, double value);

  //========================
  //GETSETGETSETGETSETGETSET
  //========================
  double getDoubleProperty(string prop);
  void setDoubleProperty(string prop, double value);

private:
	vector<double> volumeFlowRate, head;
	vector<double> coefficients; // coefficients of fitted polynomails
	int curveOrder;

	double characteristicCurve(double q);
};

#endif