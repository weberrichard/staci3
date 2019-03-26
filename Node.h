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

using namespace std;

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

  /// Setting a certain node double property
  /// demand|head|pressure|density|height|xPosition|yPosition|user1|user2
  void setProperty(string mit, double value);

  /// Getting a certain node double property
  /// demand|head|pressure|density|height|xPosition|yPosition|user1|user2
  double getProperty(string mit);

  /// Printing basic information about the node to console and log file
  string info(bool check_if_lonely);

  //========================
  //GETSETGETSETGETSETGETSET
  //========================
  string getName()
  {
    return name;
  }

private:
  string name;
  double head;
  double density;
  double xPosition, yPosition;
  double geodeticHeight;
  double demand;
  double user1, user2;
};
#endif