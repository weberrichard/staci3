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

#include "BasicFileIO.h"
#include "Edge.h"
#include "GetNew.h"
#include "IOxml.h"
#include "Node.h"
#include "PressurePoint.h"
#include "Statistic.h"
#include "Graph.h"

#include "Eigen/Eigen/Eigen"

#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace Eigen;

class Staci
{
public:
    Staci(string spr_filename);
    ~Staci();

    // Basic Node and Edge list
    vector<Node *> nodes;
    vector<Edge *> edges;

    // Creates the indicies for the nodes of edges and edges of nodes i.e. indexing of the sparse Jacobian
    void buildSystem();

    // Prints everything 
    void listSystem();

    // Converts node names (IDs) to indicies i.e. finds the node name in the list of the nodes
    vector<int> ID2Index(const vector<string> &id);

    // Checking the IDs of the edges, if one has identical ones drops exit(-1)
    void checkSystem();

    // Closing/opening edges (close: state = false | open: state = true)
    void changeEdgeStatus(string ID, bool state);

    /// Containing indicies of open edges and nodes for HydSolver
    vector<int> openEdges, openNodes;

    void addNewEdge(Edge* edge);
    
    /// ************************************************ ///
    /// GETSET GETSET GETSET GETSET GETSET GETSET GETSET ///
    /// ************************************************ ///
    string getResultFile()
    {
        return resultFile;
    }
    void setResultFile(string xml_fnev)
    {
        resultFile = xml_fnev;
    }
    void setOutputFile(string fnev)
    {
        outputFile = fnev;
    }
    string getOutputFile()
    {
        return outputFile.c_str();
    }
    string getDefinitionFile()
    {
        return definitionFile;
    }
    bool getIsInitialization()
    {
        return isInitialization;
    }
    string getInitializationFile()
    {
        return initializationFile;
    }
    void setInitializationFile(string xml_fnev)
    {
        initializationFile = xml_fnev;
    }
    void setIsInitialization(const bool van_e)
    {
        isInitialization = van_e;
    }
    string getFrictionModel(){
        return frictionModel;
    }

protected:
    /// Finds the index of value in the v vector
    int getVectorIndex(const vector<int> &v, int value);

private:
    bool isInitialization;
    string definitionFile, outputFile, initializationFile, resultFile;
    string frictionModel;
};

#endif //STACI_H