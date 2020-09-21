/*===================================================================*\
                                  Pipe
                            ---------------

  Derived from Edge class and can be instantiated. Simply pipeline
  element. Equation is a Bernoulli equation with friction losses
  and geodetic height.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef PIPE_H
#define PIPE_H

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#include "Edge.h"

class Pipe : public Edge
{

public:
  Pipe(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_length, const double a_diameter, const double a_rougness, const double a_massFlowRate, bool a_isCheckValve, int a_frictionModel, double a_relativeViscosity);
  ~Pipe();

  /// Provides basic informations
  string info();

  /// A line of F(x) = equation, rearranged to 0 in w.c.m.
  double function(const VectorXd &ppq, VectorXd &fDer);

  /// Function ferivative w.r.t. parameter (diameter / roughness)
  double functionParameterDerivative(int parameter);

  /// Initialization, mode: 0->automatic | 1-> using value
  void initialization(int mode, double value);
  
  //========================
  //GETSETGETSETGETSETGETSET
  //========================
  double getDoubleProperty(string prop);
  int getIntProperty(string prop);
  void setDoubleProperty(string prop, double value);
  void setIntProperty(string prop, int value);
  
  // Basic pipe data
  double length, diameter, lambda;

  // Calculating a constant for pressure drop calulations to reduce math operations
  double pipeConst;

  // Pipe roughness (H-W: C factor, D-W: relative roughness, C-F: friction coefficient)
  double roughness; 

  // Relative viscosity of the fluid related to 20 degree water (1e-6 m^2/s)
  double relativeViscosity;

  // 0 - Hazen-Williams (H-W), 1 - Darcy-Weisbach (D-W), 2 - Constant friction coefficient (C-F)
  int frictionModel=-1; 

private:
  // setting the pipeConst value
  void setPipeConst();

  // darcy friction coefficient factor
  double getDarcy(double q);
};

#endif