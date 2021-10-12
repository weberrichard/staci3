/*===================================================================*\
                             PressurePoint
                            ---------------

  Derived from Edge class and can be instantiated. Simply pool
  element. Prescribes a certain pressure at a Node.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef PRESSUREPOINT_H
#define PRESSUREPOINT_H

#include "Edge.h"

class PressurePoint: public Edge
{
public:
    PressurePoint(const string a_name, const double a_referenceCrossSection, const string a_startNodeName, const double a_density, const double a_head, const double a_totalHead, const double a_massFlowRate);
    ~PressurePoint();

    /// Provides basic informations
    string info();

    /// prescribed head
    double head;
    double headNominal;
    /// head can have pattern
    string patternID="";
    int patternIndex=-1;

    /// Added ammount of chlorine (mg/L)
    double waterAge = 0.;
    double addedChlorine = 0.4;
    double Cf = 1., Cb = 0.1, Sf = 1., Sb = 0.;

    /// A line of F(x) = equation, rearranged to 0 in w.c.m.
    double function(const VectorXd &ppq, VectorXd &fDer);

    /// Initialization, mode: 0->automatic | 1-> using value
    void initialization(int mode, double value);

    //========================
    //GETSETGETSETGETSETGETSET
    //========================
    double getDoubleProperty(string prop);
    void setDoubleProperty(string prop, double value);
    string getStringProperty(string prop);
    void setStringProperty(string prop, string value);
    int getIntProperty(string prop);
    void setIntProperty(string prop, int value);
};

#endif