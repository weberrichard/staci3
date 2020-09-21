#include "SegmentsAndShutdown.h"
#include <boost/tuple/tuple.hpp>
#include </mnt/d/Linux/LinuxPrograms/Gnuplot/gnuplot-iostream/gnuplot-iostream.h>
#include </usr/include/eigen3/Eigen/Eigen>
#include </usr/include/eigen3/Eigen/Dense>
#include <igraph.h>

SegmentsAndShutdown::SegmentsAndShutdown(string spr_filename) : Sensitivity(spr_filename)
{

}
SegmentsAndShutdown::~SegmentsAndShutdown(){}

void SegmentsAndShutdown::Segmentation()
{
    cout << " eddig megvan " << endl;
    vector<int> edgeVector;
    vector<int> connectingEdges;
    vector<int> wells;
    vector<int> segmentsWithTanks;
    int numberOfValves = 0;
    int changeTo = 0, changeBack = 0;
    int counter = 0 , untilThis = 1;
    double sum = 0;
    for(int i = 0; i < edges.size(); i++)
    {	
    	if(edges[i]->getEdgeStringProperty("type") == "Pipe")
    	{
    		edgeVector.push_back(edges[i]->getEdgeIntProperty("startNodeIndex"));
    		counter++;
    		edgeVector.push_back(edges[i]->getEdgeIntProperty("endNodeIndex"));
    		counter++;
    	}
    	else if(edges[i]->getEdgeStringProperty("type") == "Valve")
    	{
    		connectingEdges.push_back(edges[i]->getEdgeIntProperty("startNodeIndex"));
    		connectingEdges.push_back(edges[i]->getEdgeIntProperty("endNodeIndex"));
            numberOfValves += 1;
    	}
    	else if(edges[i]->getEdgeStringProperty("type") == "Pool" || edges[i]->getEdgeStringProperty("type") == "Press")
    	{
            segmentsWithTanks.push_back(i);   
    		wells.push_back(edges[i]->getEdgeIntProperty("startNodeIndex"));
            wells.push_back(edges[i]->getEdgeIntProperty("endNodeIndex")); 	
        }
    }
    int count = edgeVector.size();
    vector< vector<int> > connections(connectingEdges.size());
    vector<int> connectionsOriginal(connectingEdges.size());
    vector<int> valveDegreeDistributionCalculator(connectingEdges.size());
    positionOfSources(segmentsWithTanks);
    for(int i = 0; i < connectingEdges.size(); i++)
    {
    	cout << i << " .iteracio" << endl;
    	int e = 0;
    	for(int j = 0; j < edgeVector.size(); j++)
    	{
    		if(connectingEdges[i] == edgeVector[j])
    		{
    			cout << j << " .iteracio" << endl;
    			connectionsOriginal[i] = edgeVector[j];
    			e++;
    			connections[i].push_back(count + untilThis);
    			e++;
    			cout << 60 << " .iteracio" << endl;
    			edgeVector[j] = (count + untilThis);
    			untilThis++;
    			changeTo += 1;
    			cout << 64 << " .iteracio" << endl;
                valveDegreeDistributionCalculator[i] += 1;
    		}
    	}
    }
    cout << "eddig elmegy 66" << endl;
    if(edgeVector.size() % 2 != 0)
    {
        cout << "INVALID EDGE VECTOR" << endl;
        cin.get();
    }
    /*for (int i = 0; i < edgeVector.size(); ++i)
    {
        cout << " Edge vector element: " << edgeVector[i] << endl;
    }*/
    vector< vector<int> > segmentNest = segmenter(edgeVector);
    vector< vector<int> > connectionsBTWSegments;
    vector<int> segmentEdgeVector;

    for (int i = 0; i < connectionsOriginal.size(); ++i)
    {
        connectionsBTWSegments.push_back(vector<int>());
        for (int j = 0; j < connections[i].size(); ++j)
        {
            for (int k = 0; k < segmentNest.size(); ++k)
            {
                for (int l = 0; l < segmentNest[k].size(); ++l)
                {
                    if (segmentNest[k][l] == connections[i][j])
                    {
                        segmentNest[k][l] = connectionsOriginal[i];
                        connectionsBTWSegments[i].push_back(k);
                        changeBack += 1;
                    }
                }
            }
        }
    }
    if( changeTo == changeBack )
    {
        cout << "segmentation completed" << endl;
    }
    int connectionStack, connectionCounter;
    for (int i = 0; i < connectionsBTWSegments.size(); ++i)
    {
        connectionCounter = 0;
        for (int j = 0; j < connectionsBTWSegments[i].size(); ++j)
        {
            connectionCounter += 1;
            if(connectionCounter == 1)
            {
                connectionStack = connectionsBTWSegments[i][j];
            }
            if (connectionCounter == 2)
            {
                cout << "EV1: " << connectionStack;
                segmentEdgeVector.push_back(connectionStack);
                cout << "EV2: " << connectionsBTWSegments[i][j] << " ";
                segmentEdgeVector.push_back(connectionsBTWSegments[i][j]);
            }
        }
        cout << endl;
    }
    for (int i = 0; i < segmentEdgeVector.size(); ++i)
    {
        cout << "Szegmens: " << segmentEdgeVector[i] << endl;
    }
    cin.get();
    vector<int> segmentsCoordinatesX;
    vector<int> segmentsCoordinatesY;
    for (int i = 0; i < segmentNest.size(); ++i)
    {
        segmentsCoordinatesX.push_back(0);
        segmentsCoordinatesY.push_back(0);
        for (int j = 0; j < segmentNest[i].size(); ++j)
        {
    //        cout << "  " << nodes.at(segmentNest[i][j])->Get_nev() << endl;//nodes.at(segmentNest[i][j])->getEdgeDoubleProperty("xcoord");
            segmentsCoordinatesX[i] += nodes.at(segmentNest[i][j])->getProperty("xPosition");
            segmentsCoordinatesY[i] += nodes.at(segmentNest[i][j])->getProperty("yPosition");
        }
        segmentsCoordinatesX[i] /= segmentNest[i].size();
        segmentsCoordinatesY[i] /= segmentNest[i].size();
        cout << i <<". Szegmens X: " << segmentsCoordinatesX[i] << endl;
        cout << i <<". Szegmens Y: " << segmentsCoordinatesY[i] << endl;
    }
    // planner("PIPE_6290143", segmentNest);//PIPE_1351966PIPE_777945,PIPE_1386965,PIPE_1388906, PIPE_1388906,PIPE_780652,PIPE_1351966
    if(segmentEdgeVector.size() % 2 != 0)
    {
        cout << "INVALID EDGE VECTOR" << endl;
        cin.get();
    }
    for (int i = 0; i < nodes.size(); ++i)
    {
        cout << i << ". segment's demand " << nodes.at(i)->getProperty("demand") << endl;
        if(nodes.at(i)->getProperty("demand") != 0)
        {
            sum += nodes.at(i)->getProperty("demand");
        }
        cout << "The summarized demand of the system: " << sum << endl;
    }
    allGraphProperties(segmentEdgeVector, segmentNest);
    vector<double> Segment_AllDemand = Get_Segment_Demand(segmentNest);
    //Hunt_For_Critical_Segment( segmentEdgeVector, Segment_AllDemand, segmentNest,  wells, DegreeDist, segmentsCoordinatesX, segmentsCoordinatesY);
    Megzavart_szegmensek(segmentEdgeVector, Segment_AllDemand, segmentNest);
    vector<int> valveDegreeDistribution(round(valveDegreeDistributionCalculator.size()/2));
    int Tulcsordulo_Tolozarak = 0;
    counter = 0;
    for (int i = 0; i < valveDegreeDistribution.size(); ++i)
    {
        valveDegreeDistribution[i] += valveDegreeDistributionCalculator[counter];
        counter += 1;
        valveDegreeDistribution[i] += valveDegreeDistributionCalculator[counter];
        counter += 1;
        if(valveDegreeDistribution[i] > 2)
        {
            Tulcsordulo_Tolozarak += 1;
        }
    }
}

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int u, int v)
{
    adj[u].push_back(v); // Add v to u’s list.
}
 
// Prints all paths from 's' to 'd'
void Graph::printAllPaths(int s, int d)
{
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
 
    // Create an array to store paths
    int *path = new int[V];
    int path_index = 0; // Initialize path[] as empty
 
    // Initialize all vertices as not visited
    for (int i = 0; i < V; i++)
        visited[i] = false;
 
    // Call the recursive helper function to print all paths
    printAllPathsUtil(s, d, visited, path, path_index);
}
 
// A recursive function to print all paths from 'u' to 'd'.
// visited[] keeps track of vertices in current path.
// path[] stores actual vertices and path_index is current
// index in path[]
void Graph::printAllPathsUtil(int u, int d, bool visited[],
                              int path[], int &path_index)
{
    // Mark the current node and store it in path[]
    visited[u] = true;
    path[path_index] = u;
    path_index++;
 
    // If current vertex is same as destination, then print
    // current path[]
    if (u == d)
    {
        for (int i = 0; i<path_index; i++)
            cout << path[i] << " ";
        cout << endl;
    }
    else // If current vertex is not destination
    {
        // Recur for all the vertices adjacent to current vertex
        list<int>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (!visited[*i])
                printAllPathsUtil(*i, d, visited, path, path_index);
    }
 
    // Remove current vertex from path[] and mark it as unvisited
    path_index--;
    visited[u] = false;
}
/*void SegmentsAndShutdown::PlotWDS()
{
	std::vector<std::vector<std::pair<double, double> > > all_segments;
	int counter = 0;
	for(int j = 0; j < round((edgeVector.size()/2)); j++) {
		std::vector<std::pair<double, double> > segment;
		for(int i=0; i < 2; i++) {
			segment.push_back(std::make_pair(
					nodes.at(edgeVector[counter])->getProperty("xcoord"), 
					nodes.at(edgeVector[counter])->getProperty("ycoord")));
			counter++;
		}
		all_segments.push_back(segment);
	}
    //gp << "set style line 3 lt 2 lc rgb 'yellow' lw 3";
    //gp << "set style line 1  linecolor rgb 'red'  linewidth 3.000 dashtype 2 pointtype 2 pointsize default pointinterval 0 pointnumber 0";
    //gp << "set linestyle  1 lt 1 lc 7 # black-solid";
	gp << "plot '-' with lines linecolor rgb 'blue' linewidth 1.5000 title'Pontok'\n";
	gp.send2d(all_segments);
}
void SegmentsAndShutdown::PlotSegmentGraph()
{
    std::vector<std::vector<std::pair<double, double> > > sub_sub_segments;
    counter = 0;  
    for (int i = 0; i < round(segmentEdgeVector.size()/2); ++i)
    {
        //counter = segmentEdgeVector[i];
        vector< pair<double, double> > segment;
        for (int j = 0; j < 2; ++j)
        {
            //cout << i <<". Szegmens X: " << segmentsCoordinatesX[segmentEdgeVector[counter]] << endl;
            //cout << i <<". Szegmens Y: " << segmentsCoordinatesY[segmentEdgeVector[counter]] << endl;
            segment.push_back(std::make_pair(
                    segmentsCoordinatesX[segmentEdgeVector[counter]],   
                    segmentsCoordinatesY[segmentEdgeVector[counter]]));
            counter++;
        }
        sub_sub_segments.push_back(segment);
    }
    gp << "plot '-' with lines linecolor rgb 'blue' linewidth 1.5000 title'Graph' set terminal pdf set output 'DataSOk.pdf'\n";
    gp.send2d(sub_sub_segments);
    cin.get();
}*/
double SegmentsAndShutdown::GlobalSegmentClustering(vector<int> segmentEdgeVector)
{
    igraph_t graph;
    igraph_vector_t v2;
    igraph_real_t edges[segmentEdgeVector.size()], res;
    for(int i = 0; i < segmentEdgeVector.size(); i++)
    {
    //    cout << i << " edge: " << segmentEdgeVector[i] << endl;
        edges[i] = segmentEdgeVector[i];
    //    cout << "edges: " << edges[i] << " edgeVector: " << segmentEdgeVector[i] << " v: "<<VECTOR(v)[i] << endl;
    }
    igraph_vector_view(&v2, edges, sizeof(edges)/sizeof(double));
    igraph_create(&graph, &v2, 0, IGRAPH_UNDIRECTED);
    igraph_transitivity_undirected(&graph, &res, IGRAPH_TRANSITIVITY_ZERO);
    cout << "Clustering: " << res << endl;
    return res;
}
vector<double> SegmentsAndShutdown::LocalSegmentClustering(vector<int> segmentEdgeVector)
{
    igraph_t graph; 
    igraph_vector_t v2, res;
    igraph_vs_t vs;
    igraph_real_t edges[segmentEdgeVector.size()];
    igraph_vector_init(&res, 0);
    igraph_integer_t size;
    vector<double> LocalClustering;
    for(int i = 0; i < segmentEdgeVector.size(); i++)
    {
    //    cout << i << " edge: " << segmentEdgeVector[i] << endl;
        edges[i] = segmentEdgeVector[i];
    //    cout << "edges: " << edges[i] << " edgeVector: " << segmentEdgeVector[i] << " v: "<< VECTOR(v)[i] << endl;
    }
    igraph_vector_view(&v2, edges, sizeof(edges)/sizeof(double));
    igraph_create(&graph, &v2, 0, IGRAPH_UNDIRECTED);
    igraph_transitivity_local_undirected(&graph, &res, igraph_vss_all(), IGRAPH_TRANSITIVITY_ZERO);
    cout << "EDDIG MEGVAN" << endl;
    igraph_vs_all(&vs);
    igraph_vs_size(&graph, &vs, &size);
    for (int i = 0; i < size; ++i)
    {
        cout << "local clustering[" << i << "] helyen: "<< VECTOR(res)[i] << endl;
    }
    return LocalClustering;
}

double SegmentsAndShutdown::AvPathLength(vector<int> segmentEdgeVector)
{
    double AvPath = 0.;
    igraph_t graph; 
    igraph_vector_t v2;
    igraph_vs_t vs;
    igraph_real_t edges[segmentEdgeVector.size()], res;
    //igraph_vector_init(&res, 0);
    igraph_integer_t size;
    for(int i = 0; i < segmentEdgeVector.size(); i++)
    {
    //    cout << i << " edge: " << segmentEdgeVector[i] << endl;
        edges[i] = segmentEdgeVector[i];
    //    cout << "edges: " << edges[i] << " edgeVector: " << segmentEdgeVector[i] << " v: "<< VECTOR(v)[i] << endl;
    }
    igraph_vector_view(&v2, edges, sizeof(edges)/sizeof(double));
    igraph_create(&graph, &v2, 0, IGRAPH_UNDIRECTED);
    igraph_average_path_length(&graph, &res, 0, 0);
    cout << "EDDIG MEGVAN" << endl;
    AvPath = res;
    return AvPath;
}

int SegmentsAndShutdown::GraphDiameter(vector<int> segmentEdgeVector)
{
    int Diameter;   
    igraph_t graph; 
    igraph_vector_t v2, res;
    igraph_vs_t vs;
    igraph_real_t edges[segmentEdgeVector.size()];
    igraph_vector_init(&res, 0);
    igraph_integer_t result;
    igraph_integer_t size;
    double AvPathLengthNum;
    for(int i = 0; i < segmentEdgeVector.size(); i++)
    {
    //    cout << i << " edge: " << segmentEdgeVector[i] << endl;
        edges[i] = segmentEdgeVector[i];
    //    cout << "edges: " << edges[i] << " edgeVector: " << segmentEdgeVector[i] << " v: "<< VECTOR(v)[i] << endl;
    }
    igraph_vector_view(&v2, edges, sizeof(edges)/sizeof(double));
    igraph_create(&graph, &v2, 0, IGRAPH_UNDIRECTED);
    igraph_diameter(&graph, &result, 0, 0, 0, IGRAPH_UNDIRECTED, 0);
    cout << "EDDIG MEGVAN" << endl;
    igraph_vs_all(&vs);
    igraph_vs_size(&graph, &vs, &size);
    Diameter = result;
    return Diameter;
}

vector<double> SegmentsAndShutdown::DegreeDistribution(vector<int> segmentEdgeVector, bool selfloop)
{
    igraph_t graph; 
    igraph_vector_t v2, res;
    igraph_vs_t vs;
    igraph_bool_t selector = selfloop;
    igraph_real_t edges[segmentEdgeVector.size()];
    igraph_vector_init(&res, 0);
    igraph_integer_t size;
    vector<double> DegreeDistribution;
    for(int i = 0; i < segmentEdgeVector.size(); i++)
    {
    //    cout << i << " edge: " << segmentEdgeVector[i] << endl;
        edges[i] = segmentEdgeVector[i];
    //    cout << "edges: " << edges[i] << " edgeVector: " << segmentEdgeVector[i] << " v: "<< VECTOR(v)[i] << endl;
    }
    igraph_vector_view(&v2, edges, sizeof(edges)/sizeof(double));
    igraph_create(&graph, &v2, 0, IGRAPH_UNDIRECTED);
    igraph_degree(&graph, &res, igraph_vss_all(), IGRAPH_ALL, selector);
    cout << "EDDIG MEGVAN" << endl;
    igraph_vs_all(&vs);
    igraph_vs_size(&graph, &vs, &size);
    for (int i = 0; i < size; ++i)
    {
  //      cout << "Degree distribution[" << i << "] helyen: "<< VECTOR(res)[i] << endl;
        DegreeDistribution.push_back(VECTOR(res)[i]);
    }
    return DegreeDistribution;
}

double SegmentsAndShutdown::Spectral_Gap_Calculator(vector<int> segmentEdgeVector, vector< vector<int> > Szegmensek)
{   
    double Spectral_Gap;
    vector<double> Eigen_Stack;
    complex<double> Local_Stack;
    MatrixXd ConnectionMatrix = MatrixXd::Zero(Szegmensek.size(), Szegmensek.size());
    for (int i = 0; i < segmentEdgeVector.size()-1; ++i)
    {
//      cout << " sajtos: " << i << " EV{i}: " << segmentEdgeVector[i] << " EV{i+1} " << segmentEdgeVector[i+1] << endl;
        ConnectionMatrix(segmentEdgeVector[i],segmentEdgeVector[i+1]) = 1;
        ConnectionMatrix(segmentEdgeVector[i+1],segmentEdgeVector[i]) = 1;
    }
    cout << ConnectionMatrix << endl << endl;
    EigenSolver<MatrixXd> es(ConnectionMatrix);
    cout << "425 megvan" << endl;
    for (int i = 0; i < es.eigenvalues().size(); ++i)
    {
        Local_Stack = es.eigenvalues()[i];
        Eigen_Stack.push_back(abs(Local_Stack.real()));
        cout << Eigen_Stack[i] <<endl;
        cout << Local_Stack.real() << endl;
        cout << Local_Stack << endl;
    }
    sort(Eigen_Stack.begin(), Eigen_Stack.end(), [](double & a, double & b){ return a > b; } );
    Spectral_Gap = Eigen_Stack[0] - Eigen_Stack[1];
    cout << endl << "Spectral_Gap: " << Spectral_Gap << endl;
  //  cin.get();
    return Spectral_Gap;
}

vector<double> SegmentsAndShutdown::Atriculation_point(vector<int> segmentEdgeVector)
{   
    vector<double> ArtPoints;
    igraph_t graph; 
    igraph_vector_t v2, res;
    igraph_vs_t vs;
    igraph_real_t edges[segmentEdgeVector.size()];
    igraph_vector_init(&res, 0);
    igraph_integer_t size;
    for(int i = 0; i < segmentEdgeVector.size(); i++)
    {
    //    cout << i << " edge: " << segmentEdgeVector[i] << endl;
        edges[i] = segmentEdgeVector[i];
    //    cout << "edges: " << edges[i] << " edgeVector: " << segmentEdgeVector[i] << " v: "<< VECTOR(v)[i] << endl;
    }
    igraph_vector_view(&v2, edges, sizeof(edges)/sizeof(double));
    igraph_create(&graph, &v2, 0, IGRAPH_UNDIRECTED);
    cout << "EDDIG MEGVAN" << endl;
    igraph_vs_all(&vs);
    igraph_articulation_points(&graph, &res); 
    for (int i = 0; i < size; ++i)
    {
    //    cout << "Art points[" << i << "] helyen: "<< VECTOR(res)[i] << endl;
        ArtPoints.push_back(VECTOR(res)[i]);
    }
    return ArtPoints;
}

void SegmentsAndShutdown::allGraphProperties(vector<int> segmentEdgeVector, vector< vector<int> > Szegmensek)
{
    ofstream stream1;
    stream1.open("/mnt/d/Linux/Staci_On_Windows/staci-master/staci-master/Staci@HT/VIZ-SOPTVR-F-47-input_modIES_V2.dat");
    double Sg = Spectral_Gap_Calculator(segmentEdgeVector, Szegmensek);
    vector<double> Loc_Clust = LocalSegmentClustering(segmentEdgeVector);
    double GlobalClustering = GlobalSegmentClustering(segmentEdgeVector);
    double AvPath = AvPathLength(segmentEdgeVector);
    int Diameter = GraphDiameter(segmentEdgeVector);
    vector<double> ArtPoints = Atriculation_point(segmentEdgeVector);
    cout << "AvPathLength" << " , " << AvPath << " , " << "GlobalClustering" << " , " << GlobalClustering << " , " << "Spectral_Gap" << " , " << Sg << " , " << "Diameter" << " , " << Diameter << " , " << "\n";
    cin.get();
    stream1 << "i" << " , " << "Loc_Clust[i]" << " , " << "ArtPoints[i]" << "\n";
    for (int i = 0; i < ArtPoints.size(); ++i)
    {
        if(i < Loc_Clust.size())
        {
            stream1 << " , " << i << " , " << Loc_Clust[i] << " , " << ArtPoints[i] << "\n";  
        }
        else
        {
            stream1 << " , " << i << " , " << "   " << " , " << " , " << ArtPoints[i] << "\n";
        }
    }
    stream1.close();
}

vector<double> SegmentsAndShutdown::Get_Segment_Demand(vector< vector<int> > Szegmensek)
{
    double Fulldemand = 0, sum = 0;
    vector<double> Segment_demand;
    vector<int> Used_element_list;
    int Enumerate;
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        for (int j = 0; j < Szegmensek[i].size(); ++j)
        {
            Enumerate = 0;
            for (int k = 0; k < Used_element_list.size(); ++k)
            {
                if(Used_element_list[k] != Szegmensek[i][j])
                {
                    Enumerate += 1;
                }
            }
            if(Enumerate == Used_element_list.size())
            {
                Fulldemand += nodes.at(Szegmensek[i][j])->getProperty("demand");
                Used_element_list.push_back(Szegmensek[i][j]);
            }
        }
        Segment_demand.push_back(Fulldemand);
        cout << i << ". szegmens, ossz fogy: " << Fulldemand << endl;
        Fulldemand = 0;
    }
    for (int i = 0; i < Segment_demand.size(); ++i)
    {
        sum += Segment_demand[i];
    }
    cout << "Osszesen: " << sum << endl;
    return Segment_demand;
}

vector<double> SegmentsAndShutdown::Get_Segment_Length(vector< vector<int> > Szegmensek)
{
    vector<double> SumLength;
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        SumLength.push_back(0);
        for (int j = 0; j < Szegmensek[i].size(); ++j)
        {
            for (int k = 0; k < edges.size()-1; ++k)
            {
                if(edges[k]->getEdgeStringProperty("type") == "Pipe")
                {
                    if ((Szegmensek[i][j] == edges[k]->getEdgeIntProperty("startNodeIndex") && Szegmensek[i][j+1] == edges[k]->getEdgeIntProperty("endNodeIndex")) || (Szegmensek[i][j+1] == edges[k]->getEdgeIntProperty("startNodeIndex") && Szegmensek[i][j] == edges[k]->getEdgeIntProperty("endNodeIndex")))
                    {
                        SumLength[i] += edges.at(k)->getDoubleProperty("length");
                    }
                }
            }
        }
        cout << i << ".szegmens hossza: " << SumLength[i] << endl;
    }
    return SumLength;
}

vector<double> SegmentsAndShutdown::Get_Segment_Volume(vector< vector<int> > Szegmensek)
{
    vector<double> SumVolume;
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        SumVolume.push_back(0);
        for (int j = 0; j < Szegmensek[i].size(); ++j)
        {
            for (int k = 0; k < edges.size()-1; ++k)
            {
                if(edges[k]->getEdgeStringProperty("type") == "Pipe")
                {
                    if ((Szegmensek[i][j] == edges[k]->getEdgeIntProperty("startNodeIndex") && Szegmensek[i][j+1] == edges[k]->getEdgeIntProperty("endNodeIndex")) || (Szegmensek[i][j+1] == edges[k]->getEdgeIntProperty("startNodeIndex") && Szegmensek[i][j] == edges[k]->getEdgeIntProperty("endNodeIndex")))
                    {
                        SumVolume[i] += edges.at(k)->getDoubleProperty("fluidVolume");
                    }
                }
            }
        }
        cout << i << ".szegmens terfogata: " << SumVolume[i] << endl;
    }
    return SumVolume;
}


vector< vector<int> > SegmentsAndShutdown::segmenter2(vector<int> SledgeVector)
{
    //tipp: A szukseges ciklusszam number of tolozarak+1
    vector<int> Segment;
    vector<int> edgeVector = SledgeVector;
    vector< vector<int> > Segments;
    int counter = 0, numerator = 0, collector = 0, First_round = 1;
    bool Goon = true, repeat = true, Last_round = false;
    while(edgeVector.size() != 0)//22//1433villasori//445//923ferto445
    {
        repeat = true;
        cout << " Szegmens meret " << Segment.size() << " Vektor meret: " << edgeVector.size() << endl;
        if(Segment.size() == 0)
        {
            Segment.push_back(edgeVector[0]);
            Segment.push_back(edgeVector[1]);
            edgeVector.erase(edgeVector.begin());
            edgeVector.erase(edgeVector.begin());
            cout << " Szegmens meret " << Segment.size() << " Vektor meret: " << edgeVector.size() << endl;
        }
        while(repeat == true)
        { 
            repeat = false;
            cout << " Szegmens meret " << Segment.size() << " Vektor meret: " << edgeVector.size() << endl;
            for (int j = 0; j < edgeVector.size(); j++)
            {
                for (int k = 0; k < Segment.size(); k++)
                {
                    if(repeat==true)
                    {
                        break;
                    }        
        //            cout << " belep 54 " << endl;
                    if(Segment[k] == edgeVector[j] && j % 2 != 0)
                    {
                        cout << " belep 60 " << "j size: " << j << endl;
                        cout << " el1: " << edgeVector[j-1] << " el2: ";
                        cout << edgeVector[j];
                        Segment.push_back(edgeVector[j-1]);
                        Segment.push_back(edgeVector[j]);
                        edgeVector.erase(edgeVector.begin() + j-1);
                        edgeVector.erase(edgeVector.begin() + j-1);
                        repeat = true;
                    }
                    if(Segment[k] == edgeVector[j] && j % 2 == 0)
                    {
                        cout << " belep 70 " << Segment.size() << " " << edgeVector.size() << " " << j << endl;
                        cout << " el2: " << edgeVector[j] << " el1: ";
                        cout << edgeVector[j+1];
                        Segment.push_back(edgeVector[j]);
                        Segment.push_back(edgeVector[j+1]);
                        edgeVector.erase(edgeVector.begin() + j);
                        edgeVector.erase(edgeVector.begin() + j);
                        repeat = true;
                    }
                }        
            }
        }
        if(edgeVector.size() > 2 && repeat != true)
        {
            Segments.push_back(Segment);
            Segment.clear();
        }
        if(edgeVector.size() == 2 && repeat != true)
        {
            Segments.push_back(Segment);
            Segment.clear();
            cout << " belep 100 " << endl;
            cout << " bemegy " << endl;
            Segment.push_back(edgeVector[0]);
            Segment.push_back(edgeVector[1]);
            Segments.push_back(Segment);
            edgeVector.clear();
            cout << " kijon " << endl;
        }
    }
    cout << " exiting ";
    return Segments;
}

vector< vector<int> > SegmentsAndShutdown::segmenter(vector<int> SledgeVector)
{
    //tipp: A szukseges ciklusszam number of tolozarak+1
    vector<int> Segment;
    vector<int> edgeVector = SledgeVector;
    vector< vector<int> > Segments;
    int counter = 0, numerator = 0, collector = 0, First_round = 1;
    bool Goon = true, ItWorks = true, Last_round = false;
    while(Goon == true)//22//1433villasori//445//923ferto445
    {
        ItWorks = false;
        cout << " Szegmens meret " << Segment.size() << " Vektor meret: " << edgeVector.size() << endl;
        if(Segment.size() == 0 && edgeVector.size() != 0)
        {
            Segment.push_back(edgeVector[0]);
            Segment.push_back(edgeVector[1]);
            edgeVector.erase(edgeVector.begin());
            edgeVector.erase(edgeVector.begin());
            cout << " Szegmens meret " << Segment.size() << " Vektor meret: " << edgeVector.size() << endl;
            First_round = 0;
        }
        else
        {
            cout << " Szegmens meret " << Segment.size() << " Vektor meret: " << edgeVector.size() << endl;
            for (int k = 0; k < Segment.size(); k++)
            {
        //        cout << " belep 51 " << endl;
                for (int j = 0; j < edgeVector.size(); j++)
                {
        //            cout << " belep 54 " << endl;
                    if(edgeVector[j] == Segment[k] && edgeVector.size() > 2)
                    {
        //                cout << " belep 57 " << endl;
                        if(j % 2 != 0)
                        {
                            cout << " belep 60 " << "j size: " << j << endl;
                            cout << " el1: " << edgeVector[j-1] << " el2: ";
                            cout << edgeVector[j];
                            Segment.push_back(edgeVector[j-1]);
                            Segment.push_back(edgeVector[j]);
                            edgeVector.erase(edgeVector.begin() + j-1);
                            edgeVector.erase(edgeVector.begin() + j-1);
                        }
                        else
                        {
                            cout << " belep 70 " << Segment.size() << " " << edgeVector.size() << " " << j << endl;
                            cout << " el2: " << edgeVector[j] << " el1: ";
                            cout << edgeVector[j+1];
                            Segment.push_back(edgeVector[j]);
                            Segment.push_back(edgeVector[j+1]);
                            edgeVector.erase(edgeVector.begin() + j);
                            edgeVector.erase(edgeVector.begin() + j);
                        }
                    }
                    else
                    {
                        collector += 1;
                    }
                }        
            }
            if(collector = edgeVector.size()-1)
            {
                cout << " belep 85 " << endl;
                cout << " Teljesul " << endl;
                Segments.push_back(Segment);
                Segment.clear();
            }
            if(edgeVector.size() == 2)
            {
                cout << " belep 100 " << endl;
                cout << " bemegy " << endl;
                Segment.push_back(edgeVector[0]);
                Segment.push_back(edgeVector[1]);
                Segments.push_back(Segment);
                edgeVector.clear();
                cout << " kijon " << endl;
            }
        }
        cout << "eddig oke";
        if (edgeVector.size() == 0)
        {
            cout << " EV. size() = " << edgeVector.size() << endl;
            Goon = false;
        }
    }
    cout << " exiting ";
    return Segments;
}

void SegmentsAndShutdown::planner(string Kizarando_ID, vector< vector<int> > segmentNest)
{
    int cspe = 0, cspv = 0;
    for (int i = 0; i < edges.size(); ++i)
    {
        //cout << " eddig megvan " << endl;
        if (edges[i]->getEdgeStringProperty("type") == "Pipe" && Kizarando_ID == edges[i]->getEdgeStringProperty("name"))
        {
             cspe = edges[i]->getEdgeIntProperty("startNodeIndex");
             cspv = edges[i]->getEdgeIntProperty("endNodeIndex");
        }
    }
    vector<string> Tolzar_list;
    vector<int> TolzarCSPID;
    vector<double> TolzarMeret;
    for (int k = 0; k < segmentNest.size(); ++k)
    {
        for (int l = 0; l < segmentNest[k].size(); ++l)
        {
            if ((segmentNest[k][l] == cspe && segmentNest[k][l+1] == cspv) || (segmentNest[k][l] == cspe && segmentNest[k][l-1] == cspv))
            {
                for (int i = 0; i < edges.size(); ++i)
                {
                    if (edges[i]->getEdgeStringProperty("type") == "Valve")
                    {
                        for (int j = 0; j < segmentNest[k].size(); ++j)
                        {
                           if (segmentNest[k][j] == edges[i]->getEdgeIntProperty("startNodeIndex") || segmentNest[k][j] == edges[i]->getEdgeIntProperty("endNodeIndex"))
                           {
                               Tolzar_list.push_back(edges[i]->getEdgeStringProperty("name"));
                               TolzarCSPID.push_back(edges[i]->getEdgeIntProperty("startNodeIndex"));
                           }
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < Tolzar_list.size(); ++i)
    {
        for (int j = 0; j < edges.size(); ++j)
        {
           if ((edges[j]->getEdgeStringProperty("type") == "Pipe" && TolzarCSPID[i] == edges[j]->getEdgeIntProperty("startNodeIndex")) || (edges[j]->getEdgeStringProperty("type") == "Pipe" && TolzarCSPID[j] == edges[j]->getEdgeIntProperty("endNodeIndex")))
           {
                cout << "megvan" << endl;
                TolzarMeret.push_back(edges[j]->getEdgeDoubleProperty("diameter"));
                break;
           }
        }
    }
    int sorszam = 0;
    vector< vector<string> > TeljesListaID;
    vector<double> TeljesListaAtm;
    for (int i = 0; i < Tolzar_list.size(); ++i)
    {
        sorszam = 0;
        for (int j = 0; j < Tolzar_list.size(); ++j)
        {
            if (TolzarMeret[i] < TolzarMeret[j])
            {
                sorszam += 1;
            }
        }
        while(TeljesListaID.size() <= sorszam)
        {
            TeljesListaID.push_back(vector<string>());
        }
        cout << "eddig oke" << endl;
        TeljesListaID[sorszam].push_back(Tolzar_list[i]);
        TeljesListaAtm.push_back(TolzarMeret[i]);
    }
    cout << endl;
    cout << "*** A kiválasztott csőszakasz ***" << endl;
    cout << "*          " << Kizarando_ID << "         *" << endl;
    cout << "*********************************" << endl;
    cout << "* A kizárandó tolózárak listája *" << endl;
    for (int i = 0; i < TeljesListaID.size(); ++i)
    {   
        for (int j = 0; j < TeljesListaID[i].size(); ++j)
        {
            cout << "*      " << TeljesListaID[i][j] << "   D: " << TeljesListaAtm[i] <<"     *" << endl;//<< "          *" << endl;//
        }
    }
    cout << "*********************************" << endl;
}

void SegmentsAndShutdown::Hunt_For_Critical_Segment(vector<int> segmentEdgeVector, vector<double> Get_Segment_Demand, vector< vector<int> > Szegmensek, vector<int> wells, vector<double> DegreeDistribution, vector<int> segmentsCoordinatesX, vector<int> segmentsCoordinatesY)
{
    ofstream stream1;
    vector<int> ForrasosSzegmensek;
    vector< vector<int> > edgeVector_Collector;
    vector< vector<int> > ConnectionMatrix( Szegmensek.size(), vector<int> ( Szegmensek.size(), 0 ) ), ConnectionMatrix_Save;
    vector<int> segmentEdgeVector_copy = segmentEdgeVector;
    vector<int> stacker, stacker_local, Guided_EV, Guided_EV_local, Guided_EV_save, Checked;
    vector<double> CritDemand( Szegmensek.size(), 0 ),Utkozbeni_Demand( Szegmensek.size(), 0 );
    double Crit_Demand_Local = 0.0;
    vector< vector<int> > Fuggosegi_Matrix;
    int CountIt = 0;
    bool Mehet = true;
    bool MoreCycleNeeded = true;
    vector<double> Segment_Length = Get_Segment_Length(Szegmensek);
    vector<double> Segment_Volume = Get_Segment_Volume(Szegmensek);
    vector<double> Kieso_Csohossz( Szegmensek.size(), 0 ), Eleresi_Csohossz( Szegmensek.size(), 0 );
    vector< vector<double> > Tarolo;
//    cout << "EDDIG MEGVAN " << wells.size() << endl;
    //Rank_Of_Segments.push_back(vector<int>());
    CountIt = 0;
    for (int i = 0; i < wells.size(); ++i)
    {
//        cout << "EDDIG MEGVAN" << endl;
        for (int j = 0; j < Szegmensek.size(); ++j)
        {
            for (int k = 0; k < Szegmensek[j].size(); ++k)
            {
                if (Szegmensek[j][k] == wells[i])
                {
                    ForrasosSzegmensek.push_back(j);
                    cout << " A forras helye: " << j << endl;
                }
            }
        }
    }
//    cout << "Szegmens forrassal: " << ForrasosSzegmensek.size() << endl;
    for (int i = 0; i < segmentEdgeVector.size()-1; ++i)
    {
        if (i != 0)
        {
            i++;
        }
//        cout << " sajtos: " << i << " EV{i}: " << segmentEdgeVector[i] << " EV{i+1} " << segmentEdgeVector[i+1] << endl;
        ConnectionMatrix[segmentEdgeVector[i]][segmentEdgeVector[i+1]] = 1;
        ConnectionMatrix[segmentEdgeVector[i+1]][segmentEdgeVector[i]] = 1;
    }
    for (int i = 0; i < ForrasosSzegmensek.size(); ++i)
    {
        stacker.clear();
        stacker.push_back(ForrasosSzegmensek[i]);
        MoreCycleNeeded = true;
        while(MoreCycleNeeded == true)
        {
            MoreCycleNeeded = false;
            for (int j = 0; j < stacker.size(); ++j)
            {
                for (int k = 0; k < ConnectionMatrix[stacker[j]].size(); ++k)
                {
                    if(ConnectionMatrix[stacker[j]][k] == 1)
                    {
                        ConnectionMatrix[k][stacker[j]] = 0;
                        stacker_local.push_back(k);
                        MoreCycleNeeded = true;
                        //cout << "valami tortent" << endl;
                        Guided_EV.push_back(stacker[j]);
                        Guided_EV.push_back(k);
                    }
                }
            }
            stacker.clear();
            for (int l = 0; l < stacker_local.size(); ++l)
            {
                stacker.push_back(stacker_local[l]);
            }
            stacker_local.clear();
        }
    }
/*    for (int ll = 0; ll < Guided_EV.size(); ++ll)
    {
        cout << "GEV el:" << Guided_EV[ll] << endl;
    }*/
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        Guided_EV_local.clear();
    //    cout << "belep 394" << endl;
        ConnectionMatrix_Save = ConnectionMatrix;
        Guided_EV_save = Guided_EV;
    //    cout << "belep 397" << endl;
        for (int j = 0; j < Szegmensek.size(); ++j)
        {
           ConnectionMatrix_Save[i][j] = 0;
           ConnectionMatrix_Save[j][i] = 0;
        }
    //    cout << "belep 403" << endl;
        for (int i = 0; i < ForrasosSzegmensek.size(); ++i)
        {
            stacker.clear();
            stacker_local.clear();
        //    cout << "belep 406" << endl;
            stacker.push_back(ForrasosSzegmensek[i]);
        //    cout << "belep 408" << endl;
            MoreCycleNeeded = true;
            while(MoreCycleNeeded == true)
            {
                MoreCycleNeeded = false;
                for (int j = 0; j < stacker.size(); ++j)
                {
                    for (int k = 0; k < ConnectionMatrix[stacker[j]].size(); ++k)
                    {
                        if(ConnectionMatrix_Save[stacker[j]][k] == 1)
                        {
                    //        cout << " TORTENT VALAMI $$%$" << endl;
                            stacker_local.push_back(k);
                            MoreCycleNeeded = true;
                    //        cout << "valami tortent2" << endl;
                            Guided_EV_local.push_back(stacker[j]);
                            Guided_EV_local.push_back(k);
                        //    cout << "GEV piszkalva: " << stacker[j] << " k: " << k << endl;
                        }
                    }
                }
                stacker.clear();
                for (int l = 0; l < stacker_local.size(); ++l)
                {
                    stacker.push_back(stacker_local[l]);
                }
                stacker_local.clear();
            }
        }
        Checked.clear();
        for (int o = 0; o < Guided_EV_local.size(); ++o)
        {
            Mehet = true;
            for (int p = 0; p < Checked.size(); ++p)
            {
                if (Checked[p] == Guided_EV_local[o])
                {
                    Mehet = false;
                }
            }
            if (Mehet == true)
            {
                Utkozbeni_Demand[i] += Get_Segment_Demand[Guided_EV_local[o]];
                Eleresi_Csohossz[i] += Segment_Length[Guided_EV_local[o]];
                Checked.push_back(Guided_EV_local[o]);
//                cout << "GEV: " << Guided_EV_save[o] << " fogy: " << Get_Segment_Demand[Guided_EV_save[o]] << endl;
            }
        }
/*        for (int ll = 0; ll < Guided_EV_local.size(); ++ll)
        {
        cout << "GEV el loc:" << Guided_EV_local[ll] << endl;
        }*/
//        cout << "       " << endl << endl;
        //vector<int> sajt = Possible_Ways(Guided_EV, ForrasosSzegmensek, 7);
        for (int m = 0; m < Guided_EV_local.size(); ++m)
        {
            for (int n = 0; n < Guided_EV_save.size(); ++n)
            {
                if ((Guided_EV_local[m] == Guided_EV_save[n] && Guided_EV_local[m+1] == Guided_EV_save[n+1]) && n < Guided_EV_save.size()-2)
                {
                    Guided_EV_save.erase(Guided_EV_save.begin() + n);
                    Guided_EV_save.erase(Guided_EV_save.begin() + n + 1);
                }
                if ((Guided_EV_local[m] == Guided_EV_save[n] && Guided_EV_local[m+1] == Guided_EV_save[n+1]) && n >= Guided_EV_save.size()-2)
                {
                    Guided_EV_save.pop_back();
                    Guided_EV_save.pop_back();
                }
            }
        }
//        cout << "meg magvan" << endl;
        Checked.clear();
        for (int o = 0; o < Guided_EV_save.size(); ++o)
        {
            Mehet = true;
            for (int p = 0; p < Checked.size(); ++p)
            {
                if (Checked[p] == Guided_EV_save[o])
                {
                    Mehet = false;
                }
            }
            if (Mehet == true)
            {
                CritDemand[i] += Get_Segment_Demand[Guided_EV_save[o]];
                Kieso_Csohossz[i] += Segment_Length[Guided_EV_save[o]];
                Checked.push_back(Guided_EV_save[o]);
//                cout << "GEV: " << Guided_EV_save[o] << " fogy: " << Get_Segment_Demand[Guided_EV_save[o]] << endl;
            }
        }
    }
    int sumSegmentDemand = 0;
    for (int i = 0; i < Get_Segment_Demand.size(); ++i)
    {
        sumSegmentDemand += Get_Segment_Demand[i];
    }
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        Tarolo.push_back(vector<double> ());
        Tarolo[i].push_back(CritDemand[i]);
        Tarolo[i].push_back(Segment_Length[i]);
        Tarolo[i].push_back(i);
    }
    sort(Tarolo.begin(), Tarolo.end(), [](const std::vector< double >& a, const std::vector< double >& b){ return a[0] < b[0]; } );
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        Tarolo[i].push_back(0);
        for (int j = 0; j <= i; ++j)
        {
            Tarolo[i][3] += Tarolo[j][1];
        }
    }
    stream1.open("/mnt/d/Linux/Staci_On_Windows/staci-master/staci-master/Staci@HT/VIZ-SOPTVR-F-47-input_modProperties_V2.dat");
    stream1 << "i" << " , " << "CritDemand[i]" << " , " << "Eleresi_Csohossz[i]" << " , " << "Utkozbeni_Demand[i]" << " , " << "Get_Segment_Demand[i]" << " , " << "DegreeDistribution[i]" << " , " << "segmentsCoordinatesX[i]" << " , " << "segmentsCoordinatesY[i]" << " , " << "Szegmens terfogat" << " , " << "Segment_Length[i]" << " , " << "Kieso_Csohossz[i]" << " , " << "Tarolo[i][0]" << " , " << "Tarolo[i][1]" << " , " << "Tarolo[i][2]" << " , " << "Tarolo[i][3]" << " , " << "\n";
    for (int i = 0; i < CritDemand.size(); ++i)
    {
        cout << "Az "<< i << ". szegmens kiesése esetén kieső fogyasztás nagysága: " << CritDemand[i] << ", " << Tarolo[i][0] << "  " << Tarolo[i][1] << "   " << Tarolo[i][2] << " [m3/h]" << endl;
        stream1 << i << " , " << CritDemand[i] << " , " << Eleresi_Csohossz[i] << " , " << Utkozbeni_Demand[i] << " , " << Get_Segment_Demand[i] << " , " << DegreeDistribution[i] << " , " << segmentsCoordinatesX[i] << " , " << segmentsCoordinatesY[i] << " , " << Segment_Volume[i] << " , " << Segment_Length[i] << " , " << Kieso_Csohossz[i] << " , " << Tarolo[i][0] << " , " << Tarolo[i][1] << " , " << Tarolo[i][2] << " , " << Tarolo[i][3] << "\n"; 
    }
    stream1.close();
}
vector< vector<double> > SegmentsAndShutdown::Segments_Flow_Matrix(vector< vector<int> > Szegmensek)
{
    double Rho = 1000.;
    vector< vector<double> > Segments_flow_matrix(Szegmensek.size(), vector<double> ( Szegmensek.size(), 0 ));
    for (int i = 0; i < Szegmensek.size(); i++)
    {
        for(int j = 0; j < Szegmensek[i].size(); j++)
        {
            for(int k = 0; k < edges.size(); k++)
            {
                if((edges[k]->getEdgeStringProperty("type") == "Jelleggorbes fojtas" && edges[k]->getEdgeIntProperty("startNodeIndex") == Szegmensek[i][j]))
                {
                    for (int l = 0; l < Szegmensek.size(); l++)
                    {
                        for(int m = 0; m < Szegmensek[l].size(); m++)
                        {
                            if(Szegmensek[l][m] == Szegmensek[i][j] && i != l)
                            {
                                for(int o = 0; o < edges.size(); o++)
                                {
                                    if(edges[o]->getEdgeStringProperty("type") == "Pipe" && edges[o]->getEdgeIntProperty("startNodeIndex") == Szegmensek[i][j])
                                    {
                                        Segments_flow_matrix[i][l] = -1*(edges[o]->getEdgeDoubleProperty("mass_flow_rate")/Rho)*3600;
                                        Segments_flow_matrix[l][i] = (edges[o]->getEdgeDoubleProperty("mass_flow_rate")/Rho)*3600;
                                    }
                                    else if(edges[k]->getEdgeStringProperty("type") == "Pipe" && edges[k]->getEdgeIntProperty("endNodeIndex") == Szegmensek[i][j])
                                    {
                                        Segments_flow_matrix[i][l] = -1*((edges[o]->getEdgeDoubleProperty("mass_flow_rate")/Rho)*3600);
                                    }
                                }
                            }
                        }
                    }
                }
                if((edges[k]->getEdgeStringProperty("type") == "Jelleggorbes fojtas" && edges[k]->getEdgeIntProperty("endNodeIndex") == Szegmensek[i][j]))
                {
                    for (int l = 0; l < Szegmensek.size(); l++)
                    {
                        for(int m = 0; m < Szegmensek[l].size(); m++)
                        {
                            if(Szegmensek[l][m] == Szegmensek[i][j] && i != l)
                            {
                                for(int o = 0; o < edges.size(); o++)
                                {
                                    if(edges[o]->getEdgeStringProperty("type") == "Pipe" && edges[o]->getEdgeIntProperty("startNodeIndex") == Szegmensek[i][j])
                                    {
                                        Segments_flow_matrix[i][l] = (edges[o]->getEdgeDoubleProperty("mass_flow_rate")/Rho)*3600;
                                    }
                                    else if(edges[k]->getEdgeStringProperty("type") == "Pipe" && edges[k]->getEdgeIntProperty("endNodeIndex") == Szegmensek[i][j])
                                    {
                                        Segments_flow_matrix[i][l] = -1*((edges[o]->getEdgeDoubleProperty("mass_flow_rate")/Rho)*3600);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Segments_flow_matrix;
}
void SegmentsAndShutdown::positionOfSources(vector<int> segmentsWithTanks)
{
    for (int i = 0; i < segmentsWithTanks.size(); ++i)
    {
        cout << "Vegakna ID: " << segmentsWithTanks[i] << "Mass flow rate: " << edges[segmentsWithTanks[i]]->getEdgeDoubleProperty("mass_flow_rate");
        cin.get();
    }
}
void SegmentsAndShutdown::Megzavart_szegmensek(vector<int> segmentEdgeVector, vector<double> Get_Segment_Demand, vector< vector<int> > Szegmensek)
{
    ofstream stream1;
    stream1.open("/mnt/d/Linux/Staci_On_Windows/staci-master/staci-master/Staci@HT/Segment_Inflow_Distribution_VIZ-SOPTVR-F-47-input_mod_V2.dat");
    vector< vector<double> > Flow_Matrix = Segments_Flow_Matrix(Szegmensek), Flow_Matrix_Local, Inflow_Rate(Szegmensek.size(), vector<double> ( Szegmensek.size(), 0 ));
    vector<double> Segment_Total_Inflow(Flow_Matrix.size()), Segment_Total_Inflow_Local(Flow_Matrix.size()), Total_Inflow_Local(Szegmensek.size()), Total_Inflow(Szegmensek.size());
    vector<double> Segment_Total_Outflow(Flow_Matrix.size()), Segment_Total_Outflow_Local(Flow_Matrix.size()), Total_Outflow_Local(Szegmensek.size()), Total_Outflow(Szegmensek.size());
    vector<double> Segment_Total_Flow_Scale(Flow_Matrix.size()), Segment_Total_Flow_Scale_Local(Flow_Matrix.size()), Plot_Stack(Flow_Matrix.size());
    vector<double> Utak_hossza = Get_Segment_Length(Szegmensek);
    vector< vector<double> > Tarolo;
    vector<int> Drop_Out_Segments, Thirsty_Segments;
    double Conti_Corrector, Flow_Correction, Norm, All_Segment_Demand = 0.;
    complex<double> Gap = 0;
    bool Calc_Again, Checked, First_round;
    //cout << endl << "eddig megvan 744" << endl;
    for (int i = 0; i < Flow_Matrix.size(); ++i)
    {
        for (int j = 0; j < Flow_Matrix[i].size(); ++j)
        {
    //        cout << Flow_Matrix[i][j] << " ; ";
            if(Flow_Matrix[i][j] > 0)
            {
                Segment_Total_Inflow[i] += Flow_Matrix[i][j];
                //cout << endl << " Total_inflow: " << Segment_Total_Inflow[i] << endl;
            }
            if(Flow_Matrix[i][j] < 0)
            {
                Segment_Total_Outflow[i] += -1*(Flow_Matrix[i][j]);
                //cout << endl << " Total_outflow: " << Segment_Total_Outflow[i] << endl;
            }
            Segment_Total_Flow_Scale[i] += Flow_Matrix[i][j];
        }
 //       cout << endl;
    }
    cout << endl << "eddig megvan 759" << endl;
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
       // stream1 << i << " , ";
        Flow_Matrix_Local = Flow_Matrix;
        for (int j = 0; j < Flow_Matrix.size(); ++j)
        {
            Flow_Matrix_Local[i][j] = 0;
            //Flow_Matrix_Local[j][i] = 0;
        }
        for (int k = 0; k < Flow_Matrix_Local.size(); ++k)
        {
            First_round = true;
            Checked = false;
            while(Checked == false)
            {
                Calc_Again = false;
                Segment_Total_Inflow_Local[k] = 0;
                Segment_Total_Outflow_Local[k] = 0;
                Segment_Total_Flow_Scale_Local[k] = 0;
                for (int l = 0; l < Flow_Matrix_Local[k].size(); ++l)
                {
                    if(Flow_Matrix_Local[k][l] > 0)
                    {
                        Segment_Total_Inflow_Local[k] += Flow_Matrix_Local[k][l];
                        //cout << endl << " Total_inflow_Local: " << Segment_Total_Inflow_Local[k] << endl;
                    }
                    if(Flow_Matrix[k][l] < 0)
                    {
                        Segment_Total_Outflow_Local[k] += -1*(Flow_Matrix_Local[k][l]);
                    }
                    Segment_Total_Flow_Scale_Local[k] += Flow_Matrix_Local[k][l];
                    if(First_round == true)
                    {
                        for (int n = 0; n < Flow_Matrix_Local[k].size(); ++n)
                        {
                            if(Segment_Total_Inflow[n] > 0 || Segment_Total_Inflow[n] < 0)
                            {
                                Flow_Correction = Segment_Total_Inflow_Local[n]/Segment_Total_Inflow[n];
                            }
                            else
                            {
                                Flow_Correction = 0;//Segment_Total_Inflow[n]/Segment_Total_Inflow_Local[n];
                            }
                            //cout << "FC: " << Flow_Correction << endl;
                            if(isnan(Flow_Correction) == true)
                            {
                                cout << "HIBAAAA........" << Segment_Total_Inflow[n] << " Local: " << Segment_Total_Inflow_Local[n] << endl;
                            }
                            if(Flow_Matrix_Local[k][n] < 0)
                            {
                                Flow_Matrix_Local[k][n] = Flow_Correction*Flow_Matrix_Local[k][n];
                            }
                        }
                        First_round = false;
                        Calc_Again = true;
                    }
                    if (Segment_Total_Inflow_Local[k] < (Segment_Total_Outflow_Local[k] + Get_Segment_Demand[k]) && Calc_Again == false)
                    {
                        Conti_Corrector = (Segment_Total_Inflow_Local[k])/(Segment_Total_Outflow_Local[k] + Get_Segment_Demand[k]);
                        for (int m = 0; m < Flow_Matrix_Local[k].size(); ++m)
                        {
                            if(Flow_Matrix_Local[k][m] < 0)
                            {
                                Flow_Matrix_Local[k][m] = Conti_Corrector*Flow_Matrix_Local[k][m];
                                Segment_Total_Outflow_Local[k] += -1*(Flow_Matrix_Local[k][m]);
                            }
                        }
                        Checked = false;
                    }
                    else
                    {
                        Checked = true;
                    }
                }
            }
        }
        for (int i = 0; i < Szegmensek.size(); ++i)
        {
            All_Segment_Demand += Get_Segment_Demand[i];
        }

     //   cout << endl << "eddig megvan 821" << endl;
        for (int m = 0; m < Segment_Total_Inflow_Local.size(); ++m)
        {
            Total_Inflow[i] += Segment_Total_Inflow[m];
            Total_Outflow[i] += Segment_Total_Outflow[m];
            Total_Inflow_Local[i] += Segment_Total_Inflow_Local[m];
            Total_Outflow_Local[i] += Segment_Total_Outflow_Local[m];
            if (Segment_Total_Inflow_Local[m] == 0)
            {
               Drop_Out_Segments.push_back(m);
            }
            if (Segment_Total_Inflow_Local[m] < Get_Segment_Demand[m])
            {
                Thirsty_Segments.push_back(m);
            }
            if (Segment_Total_Inflow[m] > 0)
            {
                Inflow_Rate[i][m] = (Segment_Total_Outflow_Local[m]/Segment_Total_Outflow[m]);
            //    cout << Inflow_Rate[i][m] << endl << Segment_Total_Inflow[m] << endl << Segment_Total_Inflow_Local[m] << endl;
                //stream1  << Inflow_Rate[i][m] << " , ";
            }
            else
            {
               //stream1 << 0 << " , ";
            }
        }
        if(Total_Inflow[i] != 0)
        {
            //cout << endl << "eddig megvan 837" << endl;
            Plot_Stack[i] = ((Total_Outflow_Local[i] + Get_Segment_Demand[i])/(Total_Outflow[i] + All_Segment_Demand));
           // stream1 << (1 - (Total_Inflow_Local[i]/Total_Inflow[i])) << "\n";
        }
        else
        {
            //cout << endl << "eddig megvan 837" << endl;
        //    stream1 << 0 << "\n";  
        }
    }
    //cout << endl << "eddig megvan 840" << endl;
    Norm = *max_element(Plot_Stack.begin(), Plot_Stack.end());
    for (int i = 0; i < Plot_Stack.size(); ++i)
    {
    //    stream1 << i << " , ";
        Plot_Stack[i] = (Plot_Stack[i] /Norm)*100;
    /*    for (int j = 0; j < Szegmensek.size(); ++j)
        {
            stream1  << Inflow_Rate[i][j] << " , ";
        }*/
        //stream1 << Plot_Stack[i] << " , " << Utak_hossza[i] << "\n";
    }
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        Tarolo.push_back(vector<double> ());
        Tarolo[i].push_back(Plot_Stack[i]);
        Tarolo[i].push_back(Utak_hossza[i]);
        Tarolo[i].push_back(i);
    }
    sort(Tarolo.begin(), Tarolo.end(), [](const std::vector< double >& a, const std::vector< double >& b){ return a[0] < b[0]; } );
    for (int i = 0; i < Szegmensek.size(); ++i)
    {
        Tarolo[i].push_back(0);
        for (int j = 0; j <= i; ++j)
        {
            Tarolo[i][3] += Tarolo[j][1];
        }
    }
    Gap = Spectral_Gap_Calculator(segmentEdgeVector, Szegmensek);
    cin.get();
    stream1 << "Plot_Stack[i]" << " , " << "Utak_hossza[i]" << " , " << " , " << "Tarolo[i][2]" << " , " << "Tarolo[i][1]" << " , " << "Tarolo[i][0]" << " , " << "Tarolo[i][3]" << " , " << "Spectral_Gap" << "\n";
    for (int i = 0; i < Plot_Stack.size(); ++i)
    {
        stream1 << i << " , " << " , ";
        //Plot_Stack[i] = (Plot_Stack[i]/Norm)*100;
       /* for (int j = 0; j < Szegmensek.size(); ++j)
        {
            stream1  << Inflow_Rate[i][j] << " , ";
        }*/
        stream1 << Plot_Stack[i] << " , " << Utak_hossza[i] << " , " << " , " << Tarolo[i][2] << " , " << Tarolo[i][1] << " , " << Tarolo[i][0] << " , " << Tarolo[i][3] << " , " << Gap << "\n";
    }
    stream1.close();
}