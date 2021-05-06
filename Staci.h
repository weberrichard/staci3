/*===================================================================*\
                                 Staci
                            ---------------

  Main Staci class. Contains basic variables (e.g. vector for Node
  and Edge) and functions (e.g. building the system topology for
  solving the equations).
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef STACI_H
#define STACI_H

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#include "BasicFileIO.h"
#include "Edge.h"
#include "FlowMeter.h"
#include "Graph.h"
#include "IOxml.h"
#include "Node.h"
#include "Pipe.h"
#include "PressurePoint.h"
#include "Pump.h"
#include "Pool.h"
#include "Statistic.h"
#include "Valve.h"
#include "ValveISO.h"
#include "ValveFCV.h"
#include "ValvePRV.h"
#include "ValvePSV.h"
#include "ValveTCV.h"

#include </usr/include/eigen3/Eigen/Eigen>

#include <string>
#include <chrono>
#include <iomanip>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h> // mkdir

using namespace Eigen;

class Staci
{
public:
  Staci(string spr_filename);
  ~Staci();

  // Basic Node and Edge list
  vector<Node *> nodes;
  vector<Edge *> edges;
  // number of edges and nodes
  int numberEdges, numberNodes;

  // indicies vector for edge types, avoiding the for cycles
  // e.g. poolIndex contains the indicies of the pools in the edges list
  vector<int> poolIndex, presIndex, pumpIndex, valveIndex, valveISOIndex, pipeIndex, pipeCVIndex, flowMeterIndex;

  // Prints everything 
  void listSystem();

  // Constants for hydraulics, note: there are constants in Edge.h
  const double gravity = 9.81; // [m/s2]
  const double density = 1000.; // [kg/m3]
  double relativeViscosity = 1.0; // compared to 20 degree C, 1e-6; equals 1.0 by default

  // Converts node names (IDs) to indicies i.e. finds the node name in the list of the nodes
  //vector<int> ID2Index(const vector<string> &id);

  // Checking the IDs of the edges, if one has identical ones drops exit(-1)
  void checkSystem();

  /// Saving results to file
  void saveResult(string property, string element);

  /// Loading the system from INP | IOinp.cpp
  void loadSystem();
  /// Saving the system to INP | IOinp.cpp
  void saveSystem(string newFileName);

  // level of printing
  int printLevel = 0;

  int nodeIDtoIndex(string ID);
  int edgeIDtoIndex(string ID);
  
  vector<string> line2sv(string line); // cutting string line to pieces
  string frictionModel; // Darcy-Weisbach (DW) or Hazen-Williams (HW)
  
protected:

  // UNITS
  double demandUnit, headUnit; // LPS, GPM etc. to LPS
  string unit; // SI or US

  // name of the network without extension or folders
  string caseName;
  string definitionFile;

  // Creates the indicies for the nodes of edges and edges of nodes i.e. indexing of the sparse Jacobian
  void buildSystem();
  // Create the indexing, making the code more efficient
  void buildIndexing();

private:
  bool isInitialization = false;
};

#endif //STACI_H