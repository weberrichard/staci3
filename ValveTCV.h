/*===================================================================*\
                                ValveTCV
                            ----------------

  ValveTCV class derived from Valve class. Basic valve with loss curve.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef VALVETCV_H
#define VALVETCV_H

#include "Valve.h"

class ValveTCV : public Valve
{
public:
  ValveTCV(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_referenceCrossSection, vector<double> a_charX, vector<double> a_charY, double a_position, const double a_massFlowRate);
  ~ValveTCV();

  /// Provides basic informations
  string info();
  
  /// A line of F(x) = equation, rearranged to 0 in w.c.m.
  double function(vector<double> x);
  
  /// Jacobian: df/dhe, df/dhv, df/dmp
  vector<double> functionDerivative(vector<double>);

  /// Initialization, mode: 0->automatic | 1-> using value
  void initialization(int mode, double value);

  //========================
  //GETSETGETSETGETSETGETSET
  //========================
  void setDoubleProperty(string prop, double value);
  void setIntProperty(string property, int value);
  double getDoubleProperty(string prop);
  string getStringProperty(string prop);
  int getIntProperty(string property);

  vector<double> charX, charY; // X-Y coordinate of characteristic curve, that is position-dzeta
  double position, loss;
  double headLoss;
  
private:
  void updateLoss();
};

#endif