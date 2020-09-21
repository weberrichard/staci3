/*===================================================================*\
                                  Graph
                            ---------------
	
	This file provides multiple algorithms corresponding to Graphs and
	their structural properties (e.g. diameter, clustering coeff.).
	
  Independently works from staci3.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>

using namespace std;

vector<vector<int> > segmenter(vector<int> ev); //WR
double globalClustering(vector<int> ev);
double avPathLength(vector<int> ev);
int graphDiameter(vector<int> ev);
int articulationPoints(vector<int> ev);

#endif
