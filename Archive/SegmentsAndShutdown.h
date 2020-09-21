/*===================================================================*\
                                  Pool
                            ---------------

  Derived from Sensitivity class. This class is capable of
  calculating the segment graph, it's graph and hydarulical properties and 
  can create the shutdown plan of a water distribution system.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/


#include "Sensitivity.h"

class SegmentsAndShutdown : public Sensitivity
{//vector<int> edgeVector, vector<int> segmentEdgeVector, vector< vector<int> > segmentNest
private:
	vector<int> edgeVectorStack;
	vector<int> segmentEdgeVectorStack;
	vector< vector<int> > segmentNestStack;
public:
	SegmentsAndShutdown(string spr_filename);
	~SegmentsAndShutdown();
	class Graph;
	int GraphDiameter(vector<int> Segment_EV);
	double Spectral_Gap_Calculator(vector<int> Segment_EV, vector< vector<int> > Szegmensek);
	double GlobalSegmentClustering(vector<int> Segment_EV);
	double AvPathLength(vector<int> Segment_EV);
	vector<double> LocalSegmentClustering(vector<int> Segment_EV);
	vector<double> DegreeDistribution(vector<int> Segment_EV, bool selfloop);
	vector<double> Atriculation_point(vector<int> Segment_EV);
	vector<double> Get_Segment_Demand(vector< vector<int> > Szegmensek);
	vector<double> Get_Segment_Length(vector< vector<int> > Szegmensek);
	vector<double> Get_Segment_Volume(vector< vector<int> > Szegmensek);
	vector< vector<int> > segmenter(vector<int> Sledge_vector);
	vector< vector<int> > segmenter2(vector<int> Sledge_vector);
	vector< vector<double> > Segments_Flow_Matrix(vector< vector<int> > Szegmensek);
	void PlotWDS();
	void PlotSegmentGraph();
	void planner(string Kizarando_ID, vector< vector<int> > Segmensek);
	void Segmentation();
	void allGraphProperties(vector<int> Segment_EV, vector< vector<int> > Szegmensek);
	void Hunt_For_Critical_Segment(vector<int> Segment_EV, vector<double> Get_Segment_Demand, vector< vector<int> > Szegmensek, vector<int> Well, vector<double> DegreeDistribution, vector<int> Segment_CoordsX, vector<int> Segment_CoordsY);
	void positionOfSources(vector<int> Kutak);
	void Megzavart_szegmensek(vector<int> Segment_EV, vector<double> Get_Segment_Demand, vector< vector<int> > Szegmensek);
	
	void setEdgeVector(vector<int> edgeVector)
	{
		edgeVectorStack = edgeVector;
	}
	void setSegmentEdgeVector(vector<int> segmentEdgeVector)
	{
		segmentEdgeVectorStack = segmentEdgeVector;
	}
	void setSegmentNest(vector< vector<int> > segmentNest)
	{
		segmentNestStack = segmentNest;
	}
	vector<int> getEdgeVector(vector<int> edgeVector)
	{
		return edgeVectorStack;
	}
	vector<int> getSegmentEdgeVector(vector<int> segmentEdgeVector)
	{
		return segmentEdgeVectorStack;
	}
	vector< vector<int> > getSegmentNest(vector< vector<int> > segmentNest)
	{
		return segmentNestStack;
	}
};

class Graph
{
    int V;    // No. of vertices in graph
    list<int> *adj; // Pointer to an array containing adjacency lists
 
    // A recursive function used by printAllPaths()
    void printAllPathsUtil(int , int , bool [], int [], int &);
 
public:
    Graph(int V);   // Constructor
    void addEdge(int u, int v);
    void printAllPaths(int s, int d);
};