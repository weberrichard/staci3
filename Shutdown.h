/*===================================================================*\
                                Shutdown
                            ---------------
	
	This class is capable of build the segment graph from the original
	network topology. Also creates shutdown plan.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef SHUTDOWN_H
#define SHUTDOWN_H

#include "HydraulicSolver.h"

class Shutdown : public HydraulicSolver
{
public:
	Shutdown(string spr_filename);
	~Shutdown();
	void buildSegmentGraph();
	vector<int> closeDisconnectedParts(); // output: additionally closed segment indecies
	vector<int> closeDisconnectedSegments(int segmentToClose); // output: additionally closed segment indecies
	//vector<string> shutdownPlan(string pipeID);	
	vector<int> closeSegment(int segmentToClose); // output: additionally closed segment indecies
	void openEverything(); // setting every status to 1
	vector<int> getEdgeVector()
	{
		return edgeVector;
	}
	vector<int> getSegmentEdgeVector()
	{
		return segmentEdgeVector;
	}
	vector<vector<int> > getSegmentVector()
	{
		return segmentVector;
	}
	int getNumberSegment()
	{
		return numberSegment;
	}

protected:
	vector<int> edgeVector; // edge vector of the original network
	vector<int> segmentEdgeVector; // edge vector of the segment graph
	vector<vector<int> > segmentVector; // contains the original node indicies of pipes
	vector<int> segmentRank; // rank of each segment
	int numberSegment;

	void updateEdgeVector();
	vector<int> findConnectionError(vector<int> connectingNodes);
	vector<int> calculateRank(const vector<int> &ev);
};

#endif