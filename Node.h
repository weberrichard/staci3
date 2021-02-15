/*===================================================================*\
                                  Node
                            ---------------

  Node class for organizing every Node property and function.
  Staci has a vector<Node*> type variable, called nodes.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <stdlib.h>

#include </usr/include/eigen3/Eigen/Eigen>

using namespace std;
using namespace Eigen;

class Node
{
public:
  /// Containing the INgoing and OUTgoing edges, important to create the equation system
  vector<int> edgeIn;
  vector<int> edgeOut;

  Node(const string a_name, const double a_xPosition, const double a_yPosition, const double a_geodeticHeight, const double a_demand, const double a_pressure, const double a_density);
  ~Node();

  /// Giving initial values for pressure
  void initialization(int mode, double value);

  // In case of pressure dependent demands, df/dx i.e. dd/dp is not zero
  double function(const VectorXd &pq, bool isPressureDemand, VectorXd &fDer);

  // Getting the consumption in case of pressure dependent consumption
  double getConsumption(double head);

  // Calculating the function derivative with respect to the parameter for sensitivity calculation
  double functionParameterDerivative(bool isPressureDemand);

  /// Setting a certain node double property
  /// demand|head|pressure|density|height|xPosition|yPosition|user1|user2
  void setProperty(string mit, double value);

  /// Getting a certain node double property
  /// demand|head|pressure|density|height|xPosition|yPosition|user1|user2
  double getProperty(string mit);

  /// Printing basic information about the node to console and log file
  string info(bool check_if_lonely);

  /// If every connecting edges are closed, then the node will be as well, basically closed if edgeIn.size() + edgeOut.size() is zero
  //bool isClosed = false;
  int status = 1;

  // [m] for series calculations
  vector<double> vectorHead;
  // [l/s] actual served water in time
  vector<double> vectorConsumption;
  // [-] status of the nodes, 1: open, 0: closed
  vector<int> vectorStatus;

  // [l/s] for series calculations, node can have multiple demand with different pattern
  vector<double> vectorDemand;
  // patterns for series calculations
  vector<string> vectorPatternID;
  vector<int> vectorPatternIndex;

  string name;
  double head;
  double density;
  double xPosition, yPosition;
  double geodeticHeight;
  double demand; // independent from pressure, however it can be varying in time
  int segment=-1; // the node takes place in which segment
  double pdExponent = 2.5, pdDesiredPressure = 25., pdMinPressure = 10.; // in case of pressure dependent consumptions
  double userOutput;

private:
  double consumption = 0.0; // in case of presure dependent demands consumption can be smaller than demand
  double consumptionPercent = 0.0; // consumption/demand*100 [%]
  
};
#endif