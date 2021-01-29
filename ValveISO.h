/*===================================================================*\
                                ValveISO
                            ---------------

  ValveISO class derived from Valve class. Isolation valves. 
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef VALVEISO_H
#define VALVEISO_H

#include "Valve.h"

class ValveISO : public Valve
{
public:
  ValveISO(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_referenceCrossSection, const double a_massFlowRate);
  ~ValveISO();

  /// Provides basic informations
  string info();
  
  /// A line of F(x) = equation, rearranged to 0 in w.c.m.
  double function(const VectorXd &ppq, VectorXd &fDer);
  
  /// Jacobian: df/dhe, df/dhv, df/dmp
  //vector<double> functionDerivative(vector<double>);

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

  int startSegment, endSegment;
};

#endif