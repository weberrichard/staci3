#include <iostream>
#include "../../Sensitivity.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <complex> 
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include <igraph.h>
#include <stdlib.h>
#include <utility>
#include <boost/config.hpp>
#include <boost/tuple/tuple.hpp>
#include </mnt/d/Linux/LinuxPrograms/Gnuplot/gnuplot-iostream/gnuplot-iostream.h>
#include "/mnt/d/Linux/LinuxPrograms/Lemon/lemon-1.3.1/lemon/list_graph.h"
#include <lemon/dijkstra.h>
#include </usr/include/eigen3/Eigen/Eigen>
#include </usr/include/eigen3/Eigen/Dense>
#include <stdio.h>
#include <chrono> 
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>

#include "ACO.h"

//-------------------------------------------------------------------------------------//
// for testing on real: /home/namrak/HTamas/halozatok/Sopron/
// VIZ-SOPTVR-J-55-input_mod
//-------------------------------------------------------------------------------------//

using namespace std;
using namespace Eigen;
using namespace chrono;
using namespace lemon;
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//typedef std::pair<int, int> Edge;

Sensitivity *wds;
Gnuplot gp;

//------------------------------------------------------------------------Staci init----------------------------------------------------------------//
string case_folder = "/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Networks/"; // _halozatok/ name of //containing folder of staci file//Halozatok/Sopron/
string case_name = "balf";//name of staci file WITHOUT .inp_extended_Longest//balf
string critical_nodes = "nodelist_sop";
//-----------------------------------------------------------------------Staci init End-------------------------------------------------------------//


vector< vector<double> > SystemCharCurveOrig()//
{
    vector<double> PressurePoints, PressurePointsConnectingNodes;
    vector<double> ListOfDemands;
    int NumberOfDiscPoints = 10, counter;// the characteristic curve's disc points now, the demand will go from 10% to 200%
    double StepSize = 0.1, SumOfDemands, PressureScale = 300000, SummarizedPressureLoss, Territory, TerritoryModified;
    double EffectOfModification;
    bool first = true;
    int AKapcsolodCso;
    string PressurePointCsp;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    //wds->buildSystem();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        ListOfDemands.push_back(wds->nodes.at(i)->getProperty("demand"));
        SumOfDemands += wds->nodes.at(i)->getProperty("demand");
    }

    vector< vector<double> > CharacteristicCurve (NumberOfDiscPoints, vector<double> (2));
    for (int i = 0; i < CharacteristicCurve[0].size(); ++i)
    {
        CharacteristicCurve[0][i] = 0;
    }
    delete wds;
    for (int i = 1; i < NumberOfDiscPoints; ++i)
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
        //wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0));
    //    wds->Set_debug_level(0);
        //wds->buildSystem();
        SumOfDemands = 0;
        SummarizedPressureLoss = 0;
        wds->initialization();
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            wds->nodes[j]->setProperty("demand",(i*StepSize)*ListOfDemands[j]);
            SumOfDemands += wds->nodes.at(j)->getProperty("demand");
        }
        //wds->buildSystem();
        //wds->ini();
        wds->calculateSensitivity("demand");
        for (int j = 0; j < wds->edges.size(); ++j)
        {
            if (wds->edges[j]->getEdgeStringProperty("type") == "Pipe")
            {
                //cout << " belepett " << endl;
                //cout << " belepett " << j << " veszteseg: " << wds->edges.at(j)->getDoubleProperty("headloss") << endl;
                SummarizedPressureLoss += wds->edges.at(j)->getDoubleProperty("headLoss");
            }
        }
        CharacteristicCurve[i][0] = SumOfDemands/3600;//(wds->edges.at(AKapcsolodCso)->getProperty("mass_flow_rate")/1000)*3600;//(i*StepSize)*SumOfDemands;
        CharacteristicCurve[i][1] = SummarizedPressureLoss*1000*9.81;
        delete wds;
    }
    return CharacteristicCurve;
}

ACO::ACO (int nAnts, int nCities, double alpha, double beta, double q, double ro, double taumax, int initNode)
{
    NUMBEROFANTS    = nAnts;
    NUMBEROFNODES   = nCities;
    ALPHA           = alpha;
    BETA            = beta;
    Q               = q;
    RO              = ro;
    TAUMAX          = taumax;
    INITIALNODE     = initNode;

    randoms = new Randoms (21); 
}

ACO::~ACO ()
{
    for(int i=0; i<NUMBEROFNODES; i++)
    {
        delete [] GRAPH[i];
        delete [] CITIES[i];
        delete [] PHEROMONES[i];
        delete [] DELTAPHEROMONES[i];
        if(i < NUMBEROFNODES - 1)
        {
            delete [] PROBS[i];
        }
    }
    delete [] GRAPH;
    delete [] CITIES;
    delete [] PHEROMONES;
    delete [] DELTAPHEROMONES;
    delete [] PROBS;
}

void ACO::Search_Graph_Build_Up()
{
    bool check = true;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            check = true;
            for (int k = 0; k < wds->edges.size(); ++k)
            {
                if ((wds->nodes.at(i)->getName() == wds->edges.at(k)->startNodeName && wds->nodes.at(j)->getName() == wds->edges.at(k)->endNodeName) || (wds->nodes.at(i)->getName() == wds->edges.at(k)->endNodeName && wds->nodes.at(j)->getName() == wds->edges.at(k)->startNodeName))
                {
                    check = false;
                }
            }
            if( check != false)
            {
                connectNODES(i,j);
            }
        }
    }

}

void ACO::init()
{
    GRAPH           = new int*[NUMBEROFNODES];
    CITIES          = new double*[NUMBEROFNODES];
    PHEROMONES      = new double*[NUMBEROFNODES];
    DELTAPHEROMONES = new double*[NUMBEROFNODES];
    PROBS           = new double*[NUMBEROFNODES-1];
    for(int i=0; i<NUMBEROFNODES; i++)
    {
        GRAPH[i]            = new int[NUMBEROFNODES];
        CITIES[i]           = new double[2];
        PHEROMONES[i]       = new double[NUMBEROFNODES];
        DELTAPHEROMONES[i]  = new double[NUMBEROFNODES];
        PROBS[i]            = new double[2];
        for (int j=0; j<2; j++)
        {
            CITIES[i][j] = -1.0;
            PROBS[i][j]  = -1.0;
        }
        for (int j=0; j<NUMBEROFNODES; j++)
        {
            GRAPH[i][j]             = 0;
            PHEROMONES[i][j]        = 0.0;
            DELTAPHEROMONES[i][j]   = 0.0;
        }
    }   

    ROUTES = new int*[NUMBEROFANTS];
    for (int i=0; i<NUMBEROFANTS; i++)
    {
        ROUTES[i] = new int[NUMBEROFNODES];
        for (int j=0; j<NUMBEROFNODES; j++)
        {
            ROUTES[i][j] = -1;
        }
    }
    
    BESTLENGTH = (double) INT_MAX;
    BESTROUTE  = new int[NUMBEROFNODES];
    for (int i=0; i<NUMBEROFNODES; i++)
    {
        BESTROUTE[i] = -1;  
    }
}

void ACO::connectNODES (int node_in, int node_out) 
{
    GRAPH[node_in][node_out] = 1;
    PHEROMONES[node_in][node_out] = randoms -> Uniforme() * TAUMAX;
    GRAPH[node_out][node_in] = 1;
    PHEROMONES[node_out][node_in] = PHEROMONES[node_in][node_out];
}

double ACO::Fitness_Function_AVG_Sens (int node_in, int node_out)
{
    double PeakSensitivity = 0., l, Diameter = 0.1, AverageSensitivity = 0.;
    MatrixXd SensMatrix, RowSumMatrix;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    vector<double> LocalSensitivity(wds->nodes.size());
    l = sqrt(pow((wds->nodes.at(node_in)->getProperty("xPosition") - wds->nodes.at(node_out)->getProperty("xPosition")), 2) + pow((wds->nodes.at(node_in)->getProperty("yPosition") - wds->nodes.at(node_out)->getProperty("yPosition")), 2));
    wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes.at(node_in)->getName(), wds->nodes.at(node_out)->getName(), 1000, l, Diameter, 0.02, 0.0,0,2));      
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        AverageSensitivity += RowSumMatrix(i);
        if (abs(RowSumMatrix(i)) > abs(PeakSensitivity))
        {
            PeakSensitivity = RowSumMatrix(i);
        }
    }
    return AverageSensitivity;
}

double ACO::Fitness_Function_Peak_Sens (int node_in, int node_out)
{
    double LocalSensitivity = 0., PeakSensitivity = 0., l, Diameter = 0.1;
    MatrixXd SensMatrix, RowSumMatrix;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    l = sqrt(pow((wds->nodes.at(node_in)->getProperty("xPosition") - wds->nodes.at(node_out)->getProperty("xPosition")), 2) + pow((wds->nodes.at(node_in)->getProperty("yPosition") - wds->nodes.at(node_out)->getProperty("yPosition")), 2));
    wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes.at(node_in)->getName(), wds->nodes.at(node_out)->getName(), 1000, l, Diameter, 0.02, 0.0,0,2));      
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        if (abs(RowSumMatrix(i)) > abs(PeakSensitivity))
        {
            PeakSensitivity = RowSumMatrix(i);
        }
    }
    return PeakSensitivity;
}

bool ACO::exists (int node_in, int node_out) {
    return (GRAPH[node_in][node_out] == 1);
}

bool ACO::vizited (int antk, int c) {
    for (int l=0; l<NUMBEROFNODES; l++) {
        if (ROUTES[antk][l] == -1) {
            break;
        }
        if (ROUTES[antk][l] == c) {
            return true;
        }
    }
    return false;
}

double ACO::PHI (int node_in, int node_out, int antk)
{
    double ETAij = (double) pow (1 / Fitness_Function_AVG_Sens (node_in, node_out), BETA);
    double TAUij = (double) pow (PHEROMONES[node_in][node_out],   ALPHA);
    double sum = 0.0;
    for (int c=0; c<NUMBEROFNODES; c++) {
        if (exists(node_in, c)) {
            if (!vizited(antk, c)) {
                double ETA = (double) pow (1 / Fitness_Function_AVG_Sens (node_in, c), BETA);
                double TAU = (double) pow (PHEROMONES[node_in][c],   ALPHA);
                sum += ETA * TAU;
            }   
        }   
    }
    return (ETAij * TAUij) / sum;
}

double ACO::fitness_Sum (int antk) {
    double sum = 0.0;
    for (int j=0; j<NUMBEROFNODES-1; j++) {
        sum += Fitness_Function_AVG_Sens (ROUTES[antk][j], ROUTES[antk][j+1]);   
    }
    return sum;
}

int ACO::NODE ()
{
    double xi = randoms -> Uniforme();
    int i = 0;
    double sum = PROBS[i][0];
    while (sum < xi)
    {
        i++;
        sum += PROBS[i][0];
    }
    return (int) PROBS[i][1];
}

void ACO::route (int antk)
{
    ROUTES[antk][0] = INITIALNODE;
    for (int i=0; i<NUMBEROFNODES-1; i++) {        
        int node_in = ROUTES[antk][i];
        int count = 0;
        for (int c=0; c<NUMBEROFNODES; c++) {
            if (node_in == c) {
                continue;   
            }
            if (exists (node_in, c)) {
                if (!vizited (antk, c)) {
                    PROBS[count][0] = PHI (node_in, c, antk);
                    PROBS[count][1] = (double) c;
                    count++;
                }

            }
        }
        
        // deadlock
        if (0 == count) {
            return;
        }
        
        ROUTES[antk][i+1] = NODE();
    }
}

int ACO::valid (int antk, int iteration) {
    for(int i=0; i<NUMBEROFNODES-1; i++) {
        int node_in = ROUTES[antk][i];
        int node_out = ROUTES[antk][i+1];
        if (node_in < 0 || node_out < 0) {
            return -1;  
        }
        if (!exists(node_in, node_out)) {
            return -2;  
        }
        for (int j=0; j<i-1; j++) {
            if (ROUTES[antk][i] == ROUTES[antk][j]) {
                return -3;
            }   
        }
    }
    
    if (!exists (INITIALNODE, ROUTES[antk][NUMBEROFNODES-1])) {
        return -4;
    }
    
    return 0;
}

void ACO::updatePHEROMONES () {
    for (int k=0; k<NUMBEROFANTS; k++) {
        double rlength = fitness_Sum(k);
        for (int r=0; r<NUMBEROFNODES-1; r++) {
            int node_in = ROUTES[k][r];
            int node_out = ROUTES[k][r+1];
            DELTAPHEROMONES[node_in][node_out] += Q / rlength;
            DELTAPHEROMONES[node_out][node_in] += Q / rlength;
        }
    }
    for (int i=0; i<NUMBEROFNODES; i++) {
        for (int j=0; j<NUMBEROFNODES; j++) {
            PHEROMONES[i][j] = (1 - RO) * PHEROMONES[i][j] + DELTAPHEROMONES[i][j];
            DELTAPHEROMONES[i][j] = 0.0;
        }   
    }
}

void ACO::printRESULTS () {
    BESTLENGTH += Fitness_Function_AVG_Sens (BESTROUTE[NUMBEROFNODES-1], INITIALNODE);
    cout << " BEST ROUTE:" << endl;
    for (int i=0; i<NUMBEROFNODES; i++) {
        cout << BESTROUTE[i] << " ";
    }
    cout << endl << "length: " << BESTLENGTH << endl;
    
    cout << endl << " IDEAL ROUTE:" << endl;
    cout << "0 7 6 2 4 5 1 3" << endl;
    cout << "length: 127.509" << endl;
}

void ACO::printGRAPH () {
    cout << " GRAPH: " << endl;
    cout << "  | ";
    for( int i=0; i<NUMBEROFNODES; i++) {
        cout << i << " ";
    }
    cout << endl << "- | ";
    for (int i=0; i<NUMBEROFNODES; i++) {
        cout << "- ";
    }
    cout << endl;
    int count = 0;
    for (int i=0; i<NUMBEROFNODES; i++) {
        cout << i << " | ";
        for (int j=0; j<NUMBEROFNODES; j++) {
            if(i == j) {
                cout << "x ";   
            }
            else {
                cout << GRAPH[i][j] << " "; 
            }
            if (GRAPH[i][j] == 1) {
                count++;    
            }
        }
        cout << endl;
    }
    cout << endl;
    cout << "Number of connections: " << count << endl << endl;
}

void ACO::printPHEROMONES () {  
    cout << " PHEROMONES: " << endl;
    cout << "  | ";
    for (int i=0; i<NUMBEROFNODES; i++) {
        printf("%5d   ", i);
    }
    cout << endl << "- | ";
    for (int i=0; i<NUMBEROFNODES; i++) {
        cout << "--------";
    }
    cout << endl;
    for (int i=0; i<NUMBEROFNODES; i++) {
        cout << i << " | ";
        for (int j=0; j<NUMBEROFNODES; j++) {
            if (i == j) {
                printf ("%5s   ", "x");
                continue;
            }
            if (exists(i, j)) {
                printf ("%7.3f ", PHEROMONES[i][j]);
            }
            else {
                if(PHEROMONES[i][j] == 0.0) {
                    printf ("%5.0f   ", PHEROMONES[i][j]);
                }
                else {
                    printf ("%7.3f ", PHEROMONES[i][j]);
                }
            }
        }
        cout << endl;
    }
    cout << endl;
}

void ACO::optimize (int ITERATIONS) {
    for (int iterations=1; iterations<=ITERATIONS; iterations++) {
        cout << flush;
        cout << "ITERATION " << iterations << " HAS STARTED!" << endl << endl;

        for (int k=0; k<NUMBEROFANTS; k++) {
            cout << " : ant " << k << " has been released!" << endl;
            while (0 != valid(k, iterations)) {
                cout << "  :: releasing ant " << k << " again!" << endl;
                for (int i=0; i<NUMBEROFNODES; i++) {
                    ROUTES[k][i] = -1;  
                }
                route(k);
            }
            
            for (int i=0; i<NUMBEROFNODES; i++) {
                cout << ROUTES[k][i] << " ";    
            }
            cout << endl;
            
            cout << "  :: route done" << endl;
            double rlength = fitness_Sum(k);

            if (rlength < BESTLENGTH) {
                BESTLENGTH = rlength;
                for (int i=0; i<NUMBEROFNODES; i++) {
                    BESTROUTE[i] = ROUTES[k][i];
                }
            }
            cout << " : ant " << k << " has ended!" << endl;                
        }       

        cout << endl << "updating PHEROMONES . . .";
        updatePHEROMONES ();
        cout << " done!" << endl << endl;
        printPHEROMONES ();
        
        for (int i=0; i<NUMBEROFANTS; i++) {
            for (int j=0; j<NUMBEROFNODES; j++) {
                ROUTES[i][j] = -1;
            }
        }

        cout << endl << "ITERATION " << iterations << " HAS ENDED!" << endl << endl;
    }
}

vector<double> ShortestPathBetweenSources(vector< vector<int> > ListOfEdges, int NumberOfNodes, vector<double> PipeLengths, vector<int> sources)
{
    //auto start = high_resolution_clock::now();
    ListGraph Graph;
    typedef ListGraph::Node ND;
    vector< ListGraph::Node > NodeList;
    vector< ListGraph::Edge > EdgeList(ListOfEdges.size());
    ListGraph::EdgeMap<double> LengthMap(Graph);
    vector<double> ListOfSourceDistances;
    for (int i = 0; i < NumberOfNodes; ++i)
    {
        NodeList.push_back(Graph.addNode());
    }
    for (int i = 0; i < ListOfEdges.size(); ++i)
    {
        EdgeList[i] = Graph.addEdge(NodeList[ListOfEdges[i][0]], NodeList[ListOfEdges[i][1]]);
        LengthMap.set(EdgeList[i],PipeLengths[i]);
    }
    //cin.get();
    Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(Graph, LengthMap);
    for (int i = 0; i < sources.size(); ++i)
    {
        dijkstra.init();
        dijkstra.addSource(NodeList[sources[i]]);
        dijkstra.start();
        for (int j = 1; j < sources.size(); ++j)
        {
            if(i != j)
            {
                //cout << "távolság: " << dijkstra.dist(NodeList[sources[j]]) << endl;
                ListOfSourceDistances.push_back(dijkstra.dist(NodeList[j]));
            }
        }
    }
        //dijkstra.run(NodeList[73], NodeList[23]);
    /*auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "IDe Time taken by function: "
    << duration.count() << " microseconds" << endl;
    cout << "eljut ide" << endl;*/
    return ListOfSourceDistances;
}

void FullEvaluationCharCurve(bool PlotPlease, double Diameter)
{
    wds = new Sensitivity(case_folder + case_name + ".inp");
    wds->calculateSensitivity("demand");
    int NumberOfNodes = wds->nodes.size();
    vector< vector<int> > EdgeList;
    vector< vector<double> > CharacteristicCurveOriginal, CharacteristicCurveMoD;
    vector<double> LengthList;
    vector<double> OriginalSourceDifferences;
    vector<double> ModifiedSourceDifferences;
    vector<int> Sources;
    double SqDiff;
    bool mehet;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        EdgeList.push_back(vector<int> ());
        EdgeList[i].push_back(wds->edges.at(i)->getEdgeIntProperty("startNodeIndex"));
        EdgeList[i].push_back(wds->edges.at(i)->getEdgeIntProperty("endNodeIndex"));
        LengthList.push_back(500.);//wds->edges.at(i)->getDoubleProperty("length"));
    }
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if (wds->edges.at(i)->getEdgeIntProperty("numberNode") == 1)
        {
            Sources.push_back(i);
            cout << "catched" << endl;
        }
    }

    OriginalSourceDifferences = ShortestPathBetweenSources(EdgeList, NumberOfNodes, LengthList, Sources);
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
    << duration.count() << " microseconds" << endl;
    MatrixXd SensMatrix2, SensMatrix, RowSumMatrix,RowSumMatrix2;
    vector<double> LocalSensitivity(wds->nodes.size()), ModifiedLocalSensitivity(wds->nodes.size()), NodalPressuresOriginal(wds->nodes.size());
    double OriginalPressureDifference, PressureDifference, PeakSensitivityDifferenceProc, PeakSensitivity = 0., ModifiedPeakSensitivity, l, AverageSensitivity = 0., ModifiedAverageSensitivity, LocalSensitivityDifference, ModifiedLocalSensitivityDifference, AverageSensitivityDifference, AverageSensitivityDifferenceProc;
    ofstream stream1;
    ofstream stream2;
    stream1.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Optimisation/Anytown_FullEval_CharCurve.dat");
    stream2.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Optimisation/Anytown_FullEval_ProblematicNodes_kizarasos.dat");
   // stream1 << "i" << " , " << "j" << " , " << "CC[0]" << " , " << "CC[1]"<< " , " << "CC[2]" << " , " << "CC[3]" << " , " << "CC[4]" << " , " << "CC[5]" << " , " << "CC[6]" << " , " << "CC[7]" << " , " << "CC[8]" << " , " << "CC[9]" << " , " << "CC[10]" << "\n";
   // stream2 << "i" << " , " << "j" << " , " << "CC[0]" << " , " << "CC[1]"<< " , " << "CC[2]" << " , " << "CC[3]" << " , " << "CC[4]" << " , " << "CC[5]" << " , " << "CC[6]" << " , " << "CC[7]" << " , " << "CC[8]" << " , " << "CC[9]" << " , " << "CC[10]" << "\n";
    CharacteristicCurveOriginal = SystemCharCurveOrig();
    for (int i = 0; i < CharacteristicCurveOriginal.size(); ++i)
    {
        for (int j = 0; j < CharacteristicCurveOriginal[i].size(); ++j)
        {
            stream1 << CharacteristicCurveOriginal[i][j] << " , ";
        }
        stream1 << "\n";
    }
   /* for (int i = 0; i < wds->nodes.size(); ++i)
    {
        cout << "a ciklusba belep, itt tart: i[" << i << "] " << "j[" << 0 << "]" << endl;
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            mehet = false;
            SqDiff = 0.;
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            if(i != j )//&& mehet == true
            {
                CharacteristicCurveOriginal = SystemCharCurveOrig();
                if(PlotPlease == true)
                {
                    stream1 << 
                    stream2 << i  << " , " << j << " , " << CharacteristicCurveModified[][] << " , " << ;
                }
                //cin.get();
                //wds->edges.pop_back();
            }
        }
    }*/
    stream1.close();
    stream2.close();
}


void FullEvaluation(bool PlotPlease, double Diameter)
{
    wds = new Sensitivity(case_folder + case_name + ".inp");
    wds->calculateSensitivity("demand");
    int NumberOfNodes = wds->nodes.size();
    vector< vector<int> > EdgeList;
    vector<double> LengthList;
    vector<double> OriginalSourceDifferences;
    vector<double> ModifiedSourceDifferences;
    vector<int> Sources;
    double SqDiff;
    bool mehet;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        EdgeList.push_back(vector<int> ());
        EdgeList[i].push_back(wds->edges.at(i)->getEdgeIntProperty("startNodeIndex"));
        EdgeList[i].push_back(wds->edges.at(i)->getEdgeIntProperty("endNodeIndex"));
        LengthList.push_back(500.);//wds->edges.at(i)->getDoubleProperty("length"));
    }
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if (wds->edges.at(i)->getEdgeIntProperty("numberNode") == 1)
        {
            Sources.push_back(i);
            cout << "catched" << endl;
        }
    }
    /*for (int i = 0; i < wds->poolIndex.size(); ++i)
    {
        Sources.push_back(wds->poolIndex.at(i));
        cout << " medence: " << wds->poolIndex.at(i) << endl;
    }
    for (int i = 0; i < wds->pumpIndex.size(); ++i)
    {
        Sources.push_back(wds->pumpIndex.at(i));
        cout << " medence: " << wds->poolIndex.at(i) << endl;
    }
    for (int i = 0; i < wds->presIndex.size(); ++i)
    {
        Sources.push_back(wds->presIndex.at(i));
        cout << " medence: " << wds->poolIndex.at(i) << endl;
    }*/
    OriginalSourceDifferences = ShortestPathBetweenSources(EdgeList, NumberOfNodes, LengthList, Sources);
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
    << duration.count() << " microseconds" << endl;
    MatrixXd SensMatrix2, SensMatrix, RowSumMatrix,RowSumMatrix2;
    vector<double> LocalSensitivity(wds->nodes.size()), ModifiedLocalSensitivity(wds->nodes.size()), NodalPressuresOriginal(wds->nodes.size());
    double OriginalPressureDifference, PressureDifference, PeakSensitivityDifferenceProc, PeakSensitivity = 0., ModifiedPeakSensitivity, l, AverageSensitivity = 0., ModifiedAverageSensitivity, LocalSensitivityDifference, ModifiedLocalSensitivityDifference, AverageSensitivityDifference, AverageSensitivityDifferenceProc;
    ofstream stream1;
    ofstream stream2;
    stream1.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/FER-CSAPDV-1-input_mod_pressure_FullEval_kizarasos.dat");
    stream2.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/FER-CSAPDV-1-input_mod_pressure_FullEval_ProblematicNodes_kizarasos.dat");
    stream1 << "i" << " , " << "j" << " , " << "Node1" << " , " << "Node2"<< " , " << "AverageSensitivityDifferenceProc" << " , " << "AverageSensitivityDifference" << " , " << "PeakSensitivityDifferenceProc" << " , " << "LocalSensitivityDifference" << " , " << "ModifiedLocalSensitivityDifference" << " , " << "Pipelength" << " , " << "PressureDifference" << " , " << "OriginalPressureDifference" << " , " << "wds->nodes.at(i)->getProperty(pressure)" << " , " << "wds->nodes.at(i)->getProperty(pressure)" << " , " << "NodalPressuresOriginal[i]" << " , " << "NodalPressuresOriginal[j]" << " , " << "GH[i]" << " , " << "GH[j]" << " , " << "Flowrate" << " , " << "locsens[i]" << " , " << "locsens[j]" << "\n";
    stream2 << "i" << " , " << "j" << " , " << "Node1" << " , " << "Node2"<< " , " << "AverageSensitivityDifferenceProc" << " , " << "AverageSensitivityDifference" << " , " << "PeakSensitivityDifferenceProc" << " , " << "LocalSensitivityDifference" << " , " << "ModifiedLocalSensitivityDifference" << " , " << "Pipelength" << " , " << "PressureDifference" << " , " << "OriginalPressureDifference" << " , " << "wds->nodes.at(i)->getProperty(pressure)" << " , " << "wds->nodes.at(i)->getProperty(pressure)" << " , " << "NodalPressuresOriginal[i]" << " , " << "NodalPressuresOriginal[j]" << " , " << "GH[i]" << " , " << "GH[j]" << " , " << "Flowrate" << " , " << "locsens[i]" << " , " << "locsens[j]" << "\n";
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    //---------Original state sensitivity calculation------------//
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        AverageSensitivity += RowSumMatrix(i);
        NodalPressuresOriginal[i] = wds->nodes.at(i)->getProperty("pressure");
        if (abs(RowSumMatrix(i)) > abs(PeakSensitivity))
        {
            PeakSensitivity = RowSumMatrix(i);
        }
        //cout << "Av. Sens: " << AverageSensitivity << endl;
    }
    AverageSensitivity = AverageSensitivity / NumberOfNodes;

    //---------------Topology Modification full evaluation-----------------//
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        cout << "a ciklusba belep, itt tart: i[" << i << "] " << "j[" << 0 << "]" << endl;
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            mehet = false;
            SqDiff = 0.;
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            ModifiedAverageSensitivity = 0;
            ModifiedPeakSensitivity = 0;
            /*EdgeList.push_back(vector<int> ());
            EdgeList[EdgeList.size()-1].push_back(i);
            EdgeList[EdgeList.size()-1].push_back(j);
            LengthList.push_back(500.);//l);
            ModifiedSourceDifferences = ShortestPathBetweenSources(EdgeList, NumberOfNodes, LengthList, Sources);
            EdgeList.pop_back();
            LengthList.pop_back();
            for (int z = 0; z < OriginalSourceDifferences.size(); ++z)
            {
                SqDiff += (ModifiedSourceDifferences[z] - OriginalSourceDifferences[z])*(ModifiedSourceDifferences[z] - OriginalSourceDifferences[z]);
            }
            if (SqDiff == 0.)
            {
                mehet = true;
            }
            else
            {
                cout << "dobta" << endl;
            }*/
            if(i != j )//&& mehet == true
            {
                delete wds;
                wds = new Sensitivity(case_folder + case_name + ".inp");
                wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes.at(i)->getName(), wds->nodes.at(j)->getName(), 1000, l, Diameter, 0.02, 0.0,0,2));
                //---------Modified state sensitivity calculation------------//         
                wds->calculateSensitivity("demand");
                SensMatrix2 = wds->getPSensMatrix();
                RowSumMatrix2 = SensMatrix2.rowwise().sum();
                for (int i = 0; i < wds->nodes.size(); ++i)
                {
                    ModifiedLocalSensitivity[i] = RowSumMatrix2(i);
                    ModifiedAverageSensitivity += RowSumMatrix2(i);
                    if (abs(RowSumMatrix2(i)) > abs(ModifiedPeakSensitivity))
                    {
                        ModifiedPeakSensitivity = RowSumMatrix2(i);
                    }
                }
                OriginalPressureDifference = NodalPressuresOriginal[i] - NodalPressuresOriginal[j];
                PressureDifference = wds->nodes.at(i)->getProperty("pressure") - wds->nodes.at(j)->getProperty("pressure");
                ModifiedAverageSensitivity = ModifiedAverageSensitivity / NumberOfNodes;
                //cout << "Mod. Av. Sens: " << ModifiedAverageSensitivity << endl;
                LocalSensitivityDifference = abs(LocalSensitivity[i] - LocalSensitivity[j]);
                ModifiedLocalSensitivityDifference = abs(ModifiedLocalSensitivity[i] - ModifiedLocalSensitivity[j]);
                AverageSensitivityDifference = abs(AverageSensitivity) - abs(ModifiedAverageSensitivity);
                AverageSensitivityDifferenceProc = AverageSensitivityDifference/abs(AverageSensitivity);
                PeakSensitivityDifferenceProc = (abs(PeakSensitivity) - abs(ModifiedPeakSensitivity))/abs(PeakSensitivity);
                if(PlotPlease == true)
                {
                    stream1 << i  << " , " << j << " , " << wds->nodes.at(i)->getName() << " , " << wds->nodes.at(j)->getName() << " , " << AverageSensitivityDifferenceProc ,
                    stream1 << " , " << AverageSensitivityDifference << " , " << PeakSensitivityDifferenceProc << " , " <<  LocalSensitivityDifference << " , " ;
                    stream1 << ModifiedLocalSensitivityDifference << " , " << l << " , " << PressureDifference << " , " << OriginalPressureDifference << " , " ;
                    stream1 << wds->nodes.at(i)->getProperty("pressure") << " , " << wds->nodes.at(j)->getProperty("pressure") << " , " << NodalPressuresOriginal[i] << " , " ;
                    stream1 << NodalPressuresOriginal[j] << " , " << wds->nodes.at(i)->getProperty("height") << " , " << wds->nodes.at(j)->getProperty("height") << " , "; 
                    stream1 << wds->edges[wds->edges.size()-1]->getDoubleProperty("volumeFlowRate") << " , " << LocalSensitivity[i] << " , " << LocalSensitivity[j] << "\n";
                    cout << i  << " , " << j << " , " << wds->nodes.at(i)->getName() << " , " << wds->nodes.at(j)->getName() << " , " << AverageSensitivityDifferenceProc << " , " << AverageSensitivityDifference << " , " << PeakSensitivityDifferenceProc << " , " <<  LocalSensitivityDifference << " , " << ModifiedLocalSensitivityDifference << " , " << l << " , " << PressureDifference << " , " << OriginalPressureDifference << " , " << wds->nodes.at(i)->getProperty("pressure") << " , " << wds->nodes.at(j)->getProperty("pressure") << " , " << NodalPressuresOriginal[i] << " , " << NodalPressuresOriginal[j] << " , " << wds->nodes.at(i)->getProperty("height") << " , " << wds->nodes.at(j)->getProperty("height") << " , " << wds->edges[wds->edges.size()-1]->getDoubleProperty("volumeFlowRate") << "\n";
                    if (AverageSensitivityDifference < 0)
                    {
                        stream2 << i  << " , " << j << " , " << wds->nodes.at(i)->getName() << " , " << wds->nodes.at(j)->getName() << " , " << AverageSensitivityDifferenceProc ,
                        stream2 << " , " << AverageSensitivityDifference << " , " << PeakSensitivityDifferenceProc << " , " <<  LocalSensitivityDifference << " , " ;
                        stream2 << ModifiedLocalSensitivityDifference << " , " << l << " , " << PressureDifference << " , " << OriginalPressureDifference << " , " ;
                        stream2 << wds->nodes.at(i)->getProperty("pressure") << " , " << wds->nodes.at(j)->getProperty("pressure") << " , " << NodalPressuresOriginal[i] << " , " ;
                        stream2 << NodalPressuresOriginal[j] << " , " << wds->nodes.at(i)->getProperty("height") << " , " << wds->nodes.at(j)->getProperty("height") << " , "; 
                        stream2 << wds->edges[wds->edges.size()-1]->getDoubleProperty("volumeFlowRate") << " , " << LocalSensitivity[i] << " , " << LocalSensitivity[j] << "\n";
                    }
                }
                //cin.get();
                //wds->edges.pop_back();
            }
        }
    }
    stream1.close();
    stream2.close();
}

double CalculateAverageCapacity(string FileName, string case_name, bool PlotPlease)
{
    double dzeta = 0.5, AverageCapacity = 0.0;
    vector<double> firewater;
    //-----------------Read the selected nodes from txt---------------//
    ifstream ifile;
    ofstream wfile;
    ifile.open(FileName + ".txt");
    vector<string> hydrant_nodes;
    string line, temp;
    while(getline(ifile,line))
    {
        hydrant_nodes.push_back(line);
    }
    ifile.close();
    //----------------------------------------------------------------//
    if (PlotPlease == true)
    {
        wfile.open("firewater_"+case_name+".txt");
    }
    for (int i=0; i<hydrant_nodes.size(); i++)
    {
        wds = new Sensitivity(case_name + ".inp");
        int hydrant_idx = -1, idx=0;
        while(!(hydrant_idx+1) && idx<wds->nodes.size())
        {
          if(hydrant_nodes[i] == wds->nodes[idx]->getName())
            hydrant_idx = idx;
          idx++;
        }
        wds->nodes.push_back(new Node("Hout", 0., 0., wds->nodes[hydrant_idx]->getProperty("height"), 0., 0., 1000.));
        vector<double> e,zeta;
        e.push_back(0);  e.push_back(100);
        zeta.push_back(0);  zeta.push_back(dzeta);
        wds->edges.push_back(new ValveISO("Hvalve","Hout", hydrant_nodes[i], 1000., 1., 0.));
        wds->edges.push_back(new PressurePoint("Hpres", 1., hydrant_nodes[i], 1000., 100000., 500.,0.));

        //wds->ini();
        //wds->buildSystem();
        wds->solveSystem();

        firewater.push_back(wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("massFlowRate")*60.);
        if (PlotPlease == true)
        {
            wfile << firewater[i] << '\n';
        }

        delete wds;
    }
    wfile.close();
    for (int i = 0; i < firewater.size(); ++i)
    {
        AverageCapacity += firewater[i];
    }
    AverageCapacity /= firewater.size();
    return AverageCapacity; 
}

double CalculateAverageCapacity(double Pipelength, double Diameter, string StartNode, string EndNode, string FileName, string case_name, bool PlotPlease)
{
    double dzeta = 0.5, AverageCapacity;
    vector<double> firewater;
    //-----------------Read the selected nodes from txt---------------//
    ifstream ifile;
    ofstream wfile;
    ifile.open(FileName + ".txt");
    vector<string> hydrant_nodes;
    string line, temp;
    while(getline(ifile,line))
    {
        hydrant_nodes.push_back(line);
    }
    ifile.close();
    //----------------------------------------------------------------//
    if (PlotPlease == true)
    {
        wfile.open("firewater_"+case_name+".txt");
    }
    for (int i=0; i<hydrant_nodes.size(); i++)
    {
        wds = new Sensitivity(case_name + ".inp");
        wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0,0,0));
        int hydrant_idx = -1, idx=0;
        while(!(hydrant_idx+1) && idx<wds->nodes.size())
        {
          if(hydrant_nodes[i] == wds->nodes[idx]->getName())
            hydrant_idx = idx;
          idx++;
        }
        wds->nodes.push_back(new Node("Hout", 0., 0., wds->nodes[hydrant_idx]->getProperty("height"), 0., 0., 1000.));
        vector<double> e,zeta;
        e.push_back(0);  e.push_back(100);
        zeta.push_back(0);  zeta.push_back(dzeta);
        wds->addNewEdgeElement(new ValveISO("Hvalve","Hout", hydrant_nodes[i], 1000., 1., 0.));
        wds->addNewEdgeElement(new PressurePoint("Hpres", 1., hydrant_nodes[i], 1000., 100000., 500.,0.));

        //wds->ini();
        //wds->buildSystem();
        wds->solveSystem();

        firewater.push_back(wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("massFlowRate")*60.);
        if (PlotPlease == true)
        {
            wfile << firewater[i] << '\n';
        }

        delete wds;
    }
    wfile.close();
    for (int i = 0; i < firewater.size(); ++i)
    {
        AverageCapacity += firewater[i];
    }
    AverageCapacity /= firewater.size();
    return AverageCapacity; 
}

double CalculatePrice(double Pipelength, double Diameter, double NumberOfDitches)
{
    int CaseInt = Diameter*1000;
    double SummarizedPrice, DiameterSpecificPrice;
    //------Az ásás méterenkénti költsége--------//
    double PriceOfDigging = 0.; //HuF/m
    //-------------------------------------------//
    //--------Az ásási helyekhez felvonulás költsége-------------//
    double PriceOfDitch = 0.;
    //-----------------------------------------------------------//
    //----Árlista a Wavin kft. P6-os polipropilén ivóvízügyi nyomócsövekre vonatkozó árai alapján (2019.08.13.)-------------//
    switch (CaseInt)
    {
        case 25: DiameterSpecificPrice = 707;//Huf/fm//
        case 32: DiameterSpecificPrice = 858;//Huf/fm//
        case 40: DiameterSpecificPrice = 1281;//Huf/fm//
        case 50: DiameterSpecificPrice = 1990;//Huf/fm//
        case 63: DiameterSpecificPrice = 3081;//Huf/fm//
        case 75: DiameterSpecificPrice = 4459;//Huf/fm//
        case 90: DiameterSpecificPrice = 6311;//Huf/fm//
        case 110: DiameterSpecificPrice = 9504;//Huf/fm//
        case 125: DiameterSpecificPrice = 12460;//Huf/fm//
        case 140: DiameterSpecificPrice = 15664;//Huf/fm//
        case 160: DiameterSpecificPrice = 20355;//Huf/fm//
        case 200: DiameterSpecificPrice = 31822;//Huf/fm//
        case 250: DiameterSpecificPrice = 50773;//Huf/fm//
        case 315: DiameterSpecificPrice = 81811;//Huf/fm//
        case 355: DiameterSpecificPrice = 103411;//Huf/fm//
        case 400: DiameterSpecificPrice = 131425;//Huf/fm//
        case 450: DiameterSpecificPrice = 169609;//Huf/fm//
        case 500: DiameterSpecificPrice = 209310;//Huf/fm//
    }
    //---------------------------------------------------------------------------------------------------------------------//
    SummarizedPrice = (NumberOfDitches*PriceOfDitch) + (PriceOfDigging*Pipelength) + (DiameterSpecificPrice*Pipelength);
    return SummarizedPrice;
}

double IntegrateCurve(vector < vector<double> > CharacteristicCurve)
{
    double Territory = 0.;
    for (int i = 0; i < CharacteristicCurve.size()-1; i++)
    {
        Territory += ((CharacteristicCurve[i+1][0]-CharacteristicCurve[i][0])*CharacteristicCurve[i][1]) + ((CharacteristicCurve[i+1][1]-CharacteristicCurve[i][1]) * ((CharacteristicCurve[i+1][0]-CharacteristicCurve[i][0])))*0.5;
    }
    return Territory;
}



//-------Működése a rendszer jelleggörbe alatti terület kiszámításán alapszik-----------//
vector< vector<double> > SystemCharCurveMoD(double Pipelength, double Diameter, string StartNode, string EndNode, bool PlotPlease)//
{
    vector<double> PressurePoints, PressurePointsConnectingNodes;
    vector<double> ListOfDemands;
    int NumberOfDiscPoints = 100, counter;// the characteristic curve's disc points now, the demand will go from 10% to 200%
    double StepSize = 0.1, SumOfDemands, PressureScale = 300000, SummarizedPressureLoss, Territory, TerritoryModified;
    double EffectOfModification;
    bool first = true;
    int AKapcsolodCso;
    string PressurePointCsp;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    //wds->buildSystem();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        ListOfDemands.push_back(wds->nodes.at(i)->getProperty("demand"));
        SumOfDemands += wds->nodes.at(i)->getProperty("demand");
    }

    vector< vector<double> > CharacteristicCurve (NumberOfDiscPoints, vector<double> (2));
    for (int i = 0; i < CharacteristicCurve[0].size(); ++i)
    {
        CharacteristicCurve[0][i] = 0;
    }
    delete wds;
    for (int i = 1; i < NumberOfDiscPoints; ++i)
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
        //wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0));
    //    wds->Set_debug_level(0);
        //wds->buildSystem();
        SumOfDemands = 0;
        SummarizedPressureLoss = 0;
        wds->initialization();
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            wds->nodes[j]->setProperty("demand",(i*StepSize)*ListOfDemands[j]);
            SumOfDemands += wds->nodes.at(j)->getProperty("demand");
        }
        //wds->buildSystem();
        //wds->ini();
        wds->calculateSensitivity("demand");
        for (int j = 0; j < wds->edges.size(); ++j)
        {
            if (wds->edges[j]->getEdgeStringProperty("type") == "Pipe")
            {
                //cout << " belepett " << endl;
                //cout << " belepett " << j << " veszteseg: " << wds->edges.at(j)->getDoubleProperty("headloss") << endl;
                SummarizedPressureLoss += wds->edges.at(j)->getDoubleProperty("headloss");
            }
        }
        CharacteristicCurve[i][0] = SumOfDemands/3600;//(wds->edges.at(AKapcsolodCso)->getProperty("mass_flow_rate")/1000)*3600;//(i*StepSize)*SumOfDemands;
        CharacteristicCurve[i][1] = SummarizedPressureLoss*1000*9.81;
        delete wds;
    }
    //cout << endl << "eddig megvan2" << endl;
    /*for (int i = 0; i < CharacteristicCurve.size(); ++i)
    {
        for (int j = 0; j < CharacteristicCurve[i].size(); ++j)
        {
            cout << CharacteristicCurve[i][j] << " , ";
        }
        cout << endl;
    }*/
    //---------Plot the characteristic curve with gnuplot--------------//
    /*std::vector<std::pair<double, double> > curve;
    for(int j = 0; j < CharacteristicCurve.size(); ++j) {
        curve.push_back(std::make_pair(
                CharacteristicCurve[j][0], 
                CharacteristicCurve[j][1]));
    }
    gp << "plot '-' with points title 'pts_A'\n";
    gp.send1d(curve);*/
    //-----------------------------------------------------------------------------//*/
    Territory = IntegrateCurve(CharacteristicCurve);
    vector< vector<double> > CharacteristicCurveModified (NumberOfDiscPoints, vector<double> (2));
    /*for (int i = 0; i < CharacteristicCurveModified[0].size(); ++i)
    {
        CharacteristicCurveModified[0][i] = 0;
    }*/
    for (int i = 1; i < NumberOfDiscPoints; ++i)
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
    //    wds->Set_debug_level(0);
        wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0,0,0));
        wds->initialization();
        SumOfDemands = 0;
        SummarizedPressureLoss = 0;
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            wds->nodes[j]->setProperty("demand",(i*StepSize)*ListOfDemands[j]);
            SumOfDemands += wds->nodes.at(j)->getProperty("demand");
        }
        //wds->buildSystem();
        //wds->ini();
        wds->calculateSensitivity("demand");
        for (int j = 0; j < wds->edges.size(); ++j)
        {
            if (wds->edges[j]->getEdgeStringProperty("type") == "Pipe")
            {
                //cout << " belepett " << j << " veszteseg: " << wds->edges.at(j)->getDoubleProperty("headloss") << endl;
                SummarizedPressureLoss += wds->edges.at(j)->getDoubleProperty("headloss");
            }
        }
        CharacteristicCurveModified[i][0] = SumOfDemands/3600;//(wds->edges.at(AKapcsolodCso)->getProperty("mass_flow_rate")/1000)*3600;//(i*StepSize)*SumOfDemands;
        CharacteristicCurveModified[i][1] = SummarizedPressureLoss*1000*9.81;
        delete wds;
    }

    return CharacteristicCurveModified;
}

//-------Működése a rendszer jelleggörbe alatti terület kiszámításán alapszik-----------//
double CalculatSystemResistanceDifference(double Pipelength, double Diameter, string StartNode, string EndNode, bool PlotPlease)//
{
    vector<double> PressurePoints, PressurePointsConnectingNodes;
    vector<double> ListOfDemands;
    int NumberOfDiscPoints = 100, counter;// the characteristic curve's disc points now, the demand will go from 10% to 200%
    double StepSize = 0.1, SumOfDemands, PressureScale = 300000, SummarizedPressureLoss, Territory, TerritoryModified;
    double EffectOfModification;
    bool first = true;
    int AKapcsolodCso;
    string PressurePointCsp;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    //wds->buildSystem();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        ListOfDemands.push_back(wds->nodes.at(i)->getProperty("demand"));
        SumOfDemands += wds->nodes.at(i)->getProperty("demand");
    }

    vector< vector<double> > CharacteristicCurve (NumberOfDiscPoints, vector<double> (2));
    for (int i = 0; i < CharacteristicCurve[0].size(); ++i)
    {
        CharacteristicCurve[0][i] = 0;
    }
    delete wds;
    for (int i = 1; i < NumberOfDiscPoints; ++i)
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
        //wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0));
    //    wds->Set_debug_level(0);
        //wds->buildSystem();
        SumOfDemands = 0;
        SummarizedPressureLoss = 0;
        wds->initialization();
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            wds->nodes[j]->setProperty("demand",(i*StepSize)*ListOfDemands[j]);
            SumOfDemands += wds->nodes.at(j)->getProperty("demand");
        }
        //wds->buildSystem();
        //wds->ini();
        wds->calculateSensitivity("demand");
        for (int j = 0; j < wds->edges.size(); ++j)
        {
            if (wds->edges[j]->getEdgeStringProperty("type") == "Pipe")
            {
                //cout << " belepett " << endl;
                //cout << " belepett " << j << " veszteseg: " << wds->edges.at(j)->getDoubleProperty("headloss") << endl;
                SummarizedPressureLoss += wds->edges.at(j)->getDoubleProperty("headloss");
            }
        }
        CharacteristicCurve[i][0] = SumOfDemands/3600;//(wds->edges.at(AKapcsolodCso)->getProperty("mass_flow_rate")/1000)*3600;//(i*StepSize)*SumOfDemands;
        CharacteristicCurve[i][1] = SummarizedPressureLoss*1000*9.81;
        delete wds;
    }
    //cout << endl << "eddig megvan2" << endl;
    /*for (int i = 0; i < CharacteristicCurve.size(); ++i)
    {
        for (int j = 0; j < CharacteristicCurve[i].size(); ++j)
        {
            cout << CharacteristicCurve[i][j] << " , ";
        }
        cout << endl;
    }*/
    //---------Plot the characteristic curve with gnuplot--------------//
    /*std::vector<std::pair<double, double> > curve;
    for(int j = 0; j < CharacteristicCurve.size(); ++j) {
        curve.push_back(std::make_pair(
                CharacteristicCurve[j][0], 
                CharacteristicCurve[j][1]));
    }
    gp << "plot '-' with points title 'pts_A'\n";
    gp.send1d(curve);*/
    //-----------------------------------------------------------------------------//*/
    Territory = IntegrateCurve(CharacteristicCurve);
    vector< vector<double> > CharacteristicCurveModified (NumberOfDiscPoints, vector<double> (2));
    /*for (int i = 0; i < CharacteristicCurveModified[0].size(); ++i)
    {
        CharacteristicCurveModified[0][i] = 0;
    }*/
    for (int i = 1; i < NumberOfDiscPoints; ++i)
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
    //    wds->Set_debug_level(0);
        wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0,0,0));
        wds->initialization();
        SumOfDemands = 0;
        SummarizedPressureLoss = 0;
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            wds->nodes[j]->setProperty("demand",(i*StepSize)*ListOfDemands[j]);
            SumOfDemands += wds->nodes.at(j)->getProperty("demand");
        }
        //wds->buildSystem();
        //wds->ini();
        wds->calculateSensitivity("demand");
        for (int j = 0; j < wds->edges.size(); ++j)
        {
            if (wds->edges[j]->getEdgeStringProperty("type") == "Pipe")
            {
                //cout << " belepett " << j << " veszteseg: " << wds->edges.at(j)->getDoubleProperty("headloss") << endl;
                SummarizedPressureLoss += wds->edges.at(j)->getDoubleProperty("headloss");
            }
        }
        CharacteristicCurveModified[i][0] = SumOfDemands/3600;//(wds->edges.at(AKapcsolodCso)->getProperty("mass_flow_rate")/1000)*3600;//(i*StepSize)*SumOfDemands;
        CharacteristicCurveModified[i][1] = SummarizedPressureLoss*1000*9.81;
        delete wds;
    }
    //cout << endl << "eddig megvan22 PL: " << wds->edges.at(wds->edges.size()-1)->getDoubleProperty("headloss") << endl;
    /*for (int i = 0; i < CharacteristicCurveModified.size(); ++i)
    {
        //for (int j = 0; j < CharacteristicCurveModified[i].size(); ++j)
        {
            cout << CharacteristicCurveModified[i][j] << " , ";
        }
        cout << endl;
    }*/
    /*//---------Plot the characteristic curve with gnuplot--------------//
    std::vector<std::pair<double, double> > curve1;
    for(int j = 0; j < CharacteristicCurveModified.size(); ++j) {
        curve1.push_back(std::make_pair(
                CharacteristicCurveModified[j][0], 
                CharacteristicCurveModified[j][1]));
    }
    gp << "plot '-' with points title 'pts_A'\n";
    gp.send1d(curve1);
    //-----------------------------------------------------------------------------//*/
    TerritoryModified = IntegrateCurve(CharacteristicCurveModified);
    EffectOfModification = (1-(TerritoryModified/Territory))*100;
    cout << "A valtozas szazalekos erteke: " << EffectOfModification <<endl;
    //cout << Territory << " , " << TerritoryModified << endl;
    return EffectOfModification;
}

double CalculatePriceReduction(double Pipelength, double Diameter, string StartNode, string EndNode)//
{
    double EffectOfModification, Pressure, Flowrate, Performance, ModifiedPressure, ModifiedFlowrate, ModifiedPerformance;
    bool first = true;
    int counter;
    vector<int> PPPos;
    //cout << "ideér" << endl;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    wds->initialization();
    wds->calculateSensitivity("demand");
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        //cout << " type: " << wds->edges[i]->getEdgeStringProperty("type") << endl;
        if (wds->edges[i]->getEdgeStringProperty("type") == "PressurePoint" || wds->edges[i]->getEdgeStringProperty("type") == "Pool")
        {
            Pressure = wds->edges[i]->getDoubleProperty("pressure");
            Flowrate = wds->edges[i]->getDoubleProperty("volumeFlowRate");
            Performance += Pressure*Flowrate;
        }
    }
    //cout << "ideér" << endl;
    
    Performance = Pressure*Flowrate;
    //cout << "ideér" << endl;
    delete wds;

    //cout << "ideér" << endl;
    wds = new Sensitivity(case_folder + case_name + ".inp");
//    wds->Set_debug_level(0);
    wds->addNewEdgeElement(new Pipe("extra_cso", StartNode, EndNode, 1000, Pipelength, Diameter, -0.02, 0.0,0,0));
    wds->initialization();
    wds->calculateSensitivity("demand");
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        //cout << " type: " << wds->edges[i]->getEdgeStringProperty("type") << endl;
        if (wds->edges[i]->getEdgeStringProperty("type") == "PressurePoint" || wds->edges[i]->getEdgeStringProperty("type") == "Pool")
        {
            ModifiedPressure = wds->edges[i]->getDoubleProperty("pressure");
            ModifiedFlowrate = wds->edges[i]->getDoubleProperty("volumeFlowRate");
            ModifiedPerformance += ModifiedPressure*ModifiedFlowrate;
        }
    }

    delete wds;
    EffectOfModification = (1-(ModifiedPerformance/Performance))*100;
    //cout << " ATTENTION: " << EffectOfModification << endl;
    return EffectOfModification;
}

vector< vector<double> > LocalSensitivityListCreator(vector< vector<string> > AlreadyInsertedPipeNodes, vector<double> AlreadyInsertedPipeLengths, vector<double> AlreadyInsertedPipeDiameters, double MaximalPickableLength)
{
    wds = new Sensitivity(case_folder + case_name + ".inp");
    vector< vector<double> > LocalSensitivityList;
    MatrixXd SensMatrix, RowSumMatrix;
    vector<double> LocalSensitivity(wds->nodes.size());
    double l, LocalSensitivityDifference;
//    wds->Set_debug_level(0);
    if (AlreadyInsertedPipeNodes.size() != 0 && AlreadyInsertedPipeLengths.size() != 0 && AlreadyInsertedPipeDiameters.size() != 0)
    {
        for (int i = 0; i < AlreadyInsertedPipeNodes.size(); ++i)
        {
            wds->edges.push_back(new Pipe("extra_cso", AlreadyInsertedPipeNodes[i][0], AlreadyInsertedPipeNodes[i][1], 1000, AlreadyInsertedPipeLengths[i], AlreadyInsertedPipeDiameters[i], -0.02, 0.0,0,0));
        }
    }
    //wds->buildSystem();
    ////wds->ini();
    //wds->solveSystem();
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
    }
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            if ( l <= MaximalPickableLength)
            {
                LocalSensitivityDifference = abs(LocalSensitivity[i] - LocalSensitivity[j]);
                LocalSensitivityList.push_back(vector<double> ());
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(i);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(j);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(LocalSensitivityDifference);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(l);
            }
        }
    }
    sort(LocalSensitivityList.begin(), LocalSensitivityList.end(), [](const std::vector< double >& a, const std::vector< double >& b){ return a[2] > b[2]; } );
    delete wds;
    return LocalSensitivityList;
}

vector< vector<double> > LocalSensitivityListCreator(double MaximalPickableLength)
{
    wds = new Sensitivity(case_folder + case_name + ".inp");
    vector< vector<double> > LocalSensitivityList;
    vector<double> LocalSensitivity(wds->nodes.size());
    double l, LocalSensitivityDifference;
    MatrixXd SensMatrix, RowSumMatrix;
    ofstream stream2;
    stream2.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/balf_Map.txt");
    //stream1 << "Node1" << " , " << "Node2"<< " , " << "AverageSensitivityDifferenceProc" << " , " << "AverageSensitivityDifference" << " , " <<  "LocalSensitivityDifference" << " , " << "ModifiedLocalSensitivityDifference" << " , " << "Pipelength" << "\n";
//    wds->Set_debug_level(0);
    //wds->buildSystem();
    //wds->ini();
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        stream2 << wds->nodes.at(i)->getProperty("xPosition") << "," << wds->nodes.at(i)->getProperty("yPosition") << "," << LocalSensitivity[i] << endl;
    }
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            if ( l <= MaximalPickableLength)
            {
                LocalSensitivityDifference = abs(LocalSensitivity[i] - LocalSensitivity[j]);
                LocalSensitivityList.push_back(vector<double> ());
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(i);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(j);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(LocalSensitivityDifference);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(l);
            }
        }
    }
    sort(LocalSensitivityList.begin(), LocalSensitivityList.end(), [](const std::vector< double >& a, const std::vector< double >& b){ return a[2] > b[2]; } );
    cout << " loc sens dif: " << LocalSensitivityList[0][0] << " , " << LocalSensitivityList[0][1] << " , " << LocalSensitivityList[0][2] << " , " << LocalSensitivityList[0][3] << endl;
    /*for (int i = 0; i < 100; ++i)
    {
        cout << LocalSensitivityList[i][0] << " , " << LocalSensitivityList[i][1] << " , "  << LocalSensitivityList[i][2] << " , "  << LocalSensitivityList[i][3] << " , " << wds->nodes.at(LocalSensitivityList[i][0])->getProperty("xPosition") << " , " << wds->nodes.at(LocalSensitivityList[i][0])->getProperty("yPosition") << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getProperty("xPosition") << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getProperty("yPosition") << endl;
    }*/
    //cin.get();
    delete wds;
    stream2.close();
    return LocalSensitivityList;
}

vector< vector<double> > LocalSensitivityListCreator(double MaximalPickableLength, vector< vector<string> > InsertedNode, vector<double> InsertedPipeLength, vector<double> PipeDiameter)
{
    wds = new Sensitivity(case_folder + case_name + ".inp");
    vector< vector<double> > LocalSensitivityList;
    vector<double> LocalSensitivity(wds->nodes.size());
    double l, LocalSensitivityDifference;
    MatrixXd SensMatrix, RowSumMatrix;
    stringstream ss;
    //ofstream stream2;
    //stream2.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/_Map.txt");
    if (InsertedNode.size() != 0 && InsertedPipeLength.size() != 0 && PipeDiameter.size() != 0)
    {
        for (int i = 0; i < InsertedPipeLength.size(); ++i)
        {
            ss << i;
            wds->addNewEdgeElement(new Pipe("ExtraPipe_"+ss.str(), InsertedNode[i][0], InsertedNode[i][1], 1000, InsertedPipeLength[i], PipeDiameter[i], -0.02, 0.0,0,0));
            wds->calculateSensitivity("demand");
            //cout << "hozza adta: " << endl;
        }
    }
    //stream1 << "Node1" << " , " << "Node2"<< " , " << "AverageSensitivityDifferenceProc" << " , " << "AverageSensitivityDifference" << " , " <<  "LocalSensitivityDifference" << " , " << "ModifiedLocalSensitivityDifference" << " , " << "Pipelength" << "\n";
//    wds->Set_debug_level(0);
    //wds->buildSystem();
    //wds->ini();
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        //stream2 << wds->nodes.at(i)->getProperty("xPosition") << "," << wds->nodes.at(i)->getProperty("yPosition") << "," << LocalSensitivity[i] << endl;
    }
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            if ( l <= MaximalPickableLength)
            {
                LocalSensitivityDifference = abs(LocalSensitivity[i] - LocalSensitivity[j]);
                LocalSensitivityList.push_back(vector<double> ());
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(i);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(j);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(LocalSensitivityDifference);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(l);
            }
        }
    }
    sort(LocalSensitivityList.begin(), LocalSensitivityList.end(), [](const std::vector< double >& a, const std::vector< double >& b){ return a[2] > b[2]; } );
    cout << " loc sens dif: " << LocalSensitivityList[0][0] << " , " << LocalSensitivityList[0][1] << " , " << LocalSensitivityList[0][2] << " , " << LocalSensitivityList[0][3] << endl;
    /*for (int i = 0; i < 100; ++i)
    {
        cout << LocalSensitivityList[i][0] << " , " << LocalSensitivityList[i][1] << " , "  << LocalSensitivityList[i][2] << " , "  << LocalSensitivityList[i][3] << " , " << wds->nodes.at(LocalSensitivityList[i][0])->getProperty("xPosition") << " , " << wds->nodes.at(LocalSensitivityList[i][0])->getProperty("yPosition") << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getProperty("xPosition") << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getProperty("yPosition") << endl;
    }*/
    //cin.get();
    delete wds;
    //stream2.close();
    return LocalSensitivityList;
}

vector<double> GenerateResults(double NodalNumber, vector< vector<string> > InsertedNode, vector<double> InsertedPipeLength, vector<double> PipeDiameter, vector< vector<double> > LocalSensitivityList, bool PlotAll)
{
    cout << "ide eljut1" << endl;
    //delete wds; 
    //wds = new Sensitivity(case_folder + case_name + ".inp");
    cout << "ide eljut22" << endl;
    //------------------------------------------------------------//
    int NumberOfNodes1 = NodalNumber;
    vector<double> LocalSensitivity(NodalNumber), ModifiedLocalSensitivity(NodalNumber), Results;
    cout << "ide eljut2" << endl;
    double AverageSensitivity = 0., ModifiedAverageSensitivity = 0., LocalSensitivityDifference, ModifiedLocalSensitivityDifference, AverageSensitivityDifference, AverageSensitivityDifferenceProc;
    double PeakSensitivityDifferenceProc, PeakSensitivityDifference, ModifiedPeakSensitivity, PeakSensitivity, PriceOfAll = 0., TravelTimeDifference, AverageTravelTime, ModifiedAverageTravelTime, AverageCapacity, ModifiedAverageCapacity, AverageCapacityDifference, AverageCapacityDifferenceProc, CharacteristicCurveDifference = 0.;
    //-----------------------------------------------------------//
    MatrixXd SensMatrix, RowSumMatrix, ModRowSumMatrix;
    MatrixXd SensMatrixMod, RowSumMatrixMod;
    cout << "ide eljut3" << endl;
   /* wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    //---------Original state sensitivity calculation------------//
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        AverageSensitivity += RowSumMatrix(i);
        if (abs(RowSumMatrix(i)) > abs(PeakSensitivity))
        {
            PeakSensitivity = RowSumMatrix(i);
        }
        //AverageTravelTime  += wds->nodes.at(i)->getProperty("tt");
    }*/
    //---------------------------------------------------------------//
    //AverageSensitivity = -176.008;
    //PeakSensitivity = -238.632;
    //---------------------------------------------------------------//
    //cout << " ez még megvan: !!! " << endl;
   // AverageSensitivity = AverageSensitivity / NumberOfNodes1;
    cout << "Av. Sens: " << AverageSensitivity << endl;
    cout << "Peak. Sens: " << PeakSensitivity << endl;
    //cout << "Av. Sens Met2. : " << wds->CalculateAverageSensitivity() << endl;
    //cin.get();
    //AverageTravelTime = AverageTravelTime / NumberOfNodes1;
    //AverageCapacity = CalculateAverageCapacity(critical_nodes, case_name, PlotAll);
    cout << " ez még megvan: !!! " << endl;
    //-----------------------------------------------------------//
    //delete wds;
    //---------Modified state sensitivity calculation------------//
    wds = new Sensitivity(case_folder + case_name + ".inp");
//    wds->Set_debug_level(0);
    cout << "ide eljut" << endl;
    if (InsertedNode.size() != 0 && InsertedPipeLength.size() != 0 && PipeDiameter.size() != 0)
    {
        for (int i = 0; i < InsertedPipeLength.size(); ++i)
        {
            cout << " bepakolt csovek: " << InsertedNode[i][0] << " és " << InsertedNode[i][1] << endl;
            wds->addNewEdgeElement(new Pipe("Extra", InsertedNode[i][0], InsertedNode[i][1], 1000, InsertedPipeLength[i], PipeDiameter[i], -0.02, 0.0,0,0));
            //cout << "hozza adta: " << endl;
        }
    }
    cout << "ide eljut" << endl;
    /*for (int i = 0; i < wds->edges.size(); ++i)
    {
        cout << "Pipe name: " << wds->edges.at(i)->getEdgeStringProperty("name") << " type" << wds->edges.at(i)->getEdgeStringProperty("type") << "Startnode: "<< " Node_From:" << wds->edges.at(i)->getEdgeStringProperty("startNodeName") << " Node_To:" << wds->edges.at(i)->getEdgeStringProperty("endNodeName") << endl;
    }*/
    //cout << "a ciusba belep, itt tart: i[" << i << "] " << "j[" << j << "]" << endl;
    //wds->buildSystem();
    //wds->ini();
    cout << " ez még megvan: !!! " << endl;
    wds->calculateSensitivity("demand");
    //wds->listSystem();
 //  cin.get();
 //   cout << "Nev: " << wds->edges[wds->edges.size()-50]->getEdgeStringProperty("name") << endl;
 //   cout << "Nev: " << wds->edges[wds->edges.size()-50]->getEdgeStringProperty("startNodeName") << endl;
 //   cout << "Van terfar: " << wds->edges[wds->edges.size()-50]->getDoubleProperty("volumeFlowRate") << endl;
 //   cout << "Nev: " << wds->edges[wds->edges.size()-1]->getEdgeStringProperty("name") << endl;
 //   cout << "Nev: " << wds->edges[wds->edges.size()-1]->getEdgeStringProperty("startNodeName") << endl;
 //   cout << "Van terfar: " << wds->edges[wds->edges.size()-1]->getDoubleProperty("volumeFlowRate") << endl;
    SensMatrixMod = wds->getPSensMatrix();
    ModRowSumMatrix = SensMatrixMod.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        ModifiedLocalSensitivity[i] = ModRowSumMatrix(i);
        ModifiedAverageSensitivity += ModRowSumMatrix(i);
        if (abs(ModRowSumMatrix(i)) > abs(ModifiedPeakSensitivity))
        {
            ModifiedPeakSensitivity = ModRowSumMatrix(i);
        }
        //ModifiedAverageTravelTime  += wds->nodes.at(i)->getProperty("tt");
    }
    //cout << " ez megvan 1" << endl;
    //ModifiedAverageTravelTime = ModifiedAverageTravelTime / NumberOfNodes1;
    ModifiedAverageSensitivity = ModifiedAverageSensitivity / NumberOfNodes1;
    cout << "Peak. Sens: " << PeakSensitivity << endl;
    cout << "Mod. Peak. Sens: " << ModifiedPeakSensitivity << endl;
    LocalSensitivityDifference = abs(LocalSensitivity[LocalSensitivityList[0][0]] - LocalSensitivity[LocalSensitivityList[0][1]]);
    ModifiedLocalSensitivityDifference = abs(ModifiedLocalSensitivity[LocalSensitivityList[0][0]] - ModifiedLocalSensitivity[LocalSensitivityList[0][1]]);
    AverageSensitivityDifference = (1 - (ModifiedAverageSensitivity/AverageSensitivity));
    cout << "Av. Sens: " << AverageSensitivity << endl;
    cout << "Mod. Av. Sens: " << ModifiedAverageSensitivity << endl;
    //cout << "Mod. Av. Sens Met2. : " << wds->CalculateAverageSensitivity() << endl;
    AverageSensitivityDifferenceProc = (1 - (ModifiedAverageSensitivity/AverageSensitivity))*100;
    cout << "Av. Sens. Diff. Proc: " << AverageSensitivityDifferenceProc << endl;
    PeakSensitivityDifference = (1 - (ModifiedPeakSensitivity/PeakSensitivity));
    PeakSensitivityDifferenceProc = (1 - (ModifiedPeakSensitivity/PeakSensitivity))*100;
    cout << "Peak. Sens. Diff. Proc: " << PeakSensitivityDifferenceProc << endl;
    //ModifiedAverageCapacity = CalculateAverageCapacity(InsertedPipeLength[0], PipeDiameter[0], wds->nodes.at(LocalSensitivityList[0][0])->getName(), wds->nodes.at(LocalSensitivityList[0][1])->getName(), critical_nodes, case_name, PlotAll);
    //AverageCapacityDifference = (1 - (AverageCapacity/ModifiedAverageCapacity));
    //AverageCapacityDifferenceProc = (1 - (AverageCapacity/ModifiedAverageCapacity))*100;
    //PriceOfAll = CalculatePriceReduction(InsertedPipeLength[0], PipeDiameter[0], wds->nodes.at(LocalSensitivityList[0][0])->getName(), wds->nodes.at(LocalSensitivityList[0][1])->getName());
    //cout << " ez megvan 3" << endl;
    //CharacteristicCurveDifference = CalculatSystemResistanceDifference(InsertedPipeLength[0], PipeDiameter[0], wds->nodes.at(LocalSensitivityList[0][0])->getName(), wds->nodes.at(LocalSensitivityList[0][1])->getName(), PlotAll);
    //cout << " ez megvan 4" << endl;
    //PriceOfAll = CalculatePriceReduction(InsertedPipeLength[0], PipeDiameter[0], wds->nodes.at(LocalSensitivityList[0][0])->getName(), wds->nodes.at(LocalSensitivityList[0][1])->getName());
    //cout << " ez megvan 5" << endl;
    //TravelTimeDifference = (1-(ModifiedAverageTravelTime/AverageTravelTime));
    //cout << " TT : " << TravelTimeDifference << endl;
    //-----------------------------------------------------------//
    //stream1 << wds->nodes.at(LocalSensitivityList[i][0])->getName() << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getName() << " , " << AverageSensitivityDifferenceProc << " , " << AverageSensitivityDifference << " , " <<  LocalSensitivityDifference << " , " << ModifiedLocalSensitivityDifference << " , " << AverageCapacityDifferenceProc << " , " << l << " , " << CharacteristicCurveDifference << " , " << TravelTimeDifference << "\n";
    //cout << wds->nodes.at(LocalSensitivityList[0][0])->getName() << " , " << wds->nodes.at(LocalSensitivityList[0][1])->getName() << " , " << AverageSensitivityDifferenceProc << " , " << AverageSensitivityDifference << " , " <<  LocalSensitivityDifference << " , " << ModifiedLocalSensitivityDifference << " , " << InsertedPipeLength[0] << " , " << CharacteristicCurveDifference << "\n";//<< TravelTimeDifference
    Results.push_back(AverageSensitivityDifferenceProc);
    Results.push_back(PeakSensitivityDifferenceProc);
    //Results.push_back(AverageCapacityDifferenceProc);
    //Results.push_back(CharacteristicCurveDifference);
    //Results.push_back(TravelTimeDifference);
    Results.push_back(PriceOfAll);
    //-----------------------------------------------------------//
    //wds->edges.pop_back();
    cout << "megvan 3" << endl;
    //wds->popLastEdgeElement();
    //delete wds;
    return Results;
}

double CalculateSummarizedPipeLength()
{
    double SummarizedPipeLength = 0.;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if (wds->edges[i]->getEdgeStringProperty("type") == "Pipe")
        {
            SummarizedPipeLength += wds->edges.at(i)->getDoubleProperty("length");
        }
    }
    return SummarizedPipeLength;
}

double CalculateSummarizedNodalDemand()
{
    double SummarizedNodalDemand = 0.;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        SummarizedNodalDemand += wds->nodes.at(i)->getProperty("demand");
    }
    return SummarizedNodalDemand;
}

void MapsForPlot()
{
    wds = new Sensitivity(case_folder + case_name + ".inp");
    vector< vector<double> > LocalSensitivityList;
    vector<double> LocalSensitivity(wds->nodes.size());
    double l, LocalSensitivityDifference;
    MatrixXd SensMatrix, RowSumMatrix;
    ofstream stream2;
    ofstream stream3;
    ofstream stream4;
    ofstream stream5;
    stream2.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/balf_Sensitivity_Map.txt");
    stream3.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/balf_Pressure_Map.txt");
    stream4.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/balf_GeodeticHeight_Map.txt");
    stream5.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/balf_Demand_Map.txt");
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        stream2 << wds->nodes.at(i)->getProperty("xPosition") << "," << wds->nodes.at(i)->getProperty("yPosition") << "," << LocalSensitivity[i] << endl;
        stream3 << wds->nodes.at(i)->getProperty("xPosition") << "," << wds->nodes.at(i)->getProperty("yPosition") << "," << wds->nodes.at(i)->getProperty("head") << endl;
        stream4 << wds->nodes.at(i)->getProperty("xPosition") << "," << wds->nodes.at(i)->getProperty("yPosition") << "," << wds->nodes.at(i)->getProperty("height") << endl;
        stream5 << wds->nodes.at(i)->getProperty("xPosition") << "," << wds->nodes.at(i)->getProperty("yPosition") << "," << wds->nodes.at(i)->getProperty("demand") << endl;
    }
    delete wds;
    stream2.close();
    stream3.close();
    stream4.close();
    stream5.close();
}

/*Graph GenerateGraph(vector< vector<string> > EdgeList, int NumberOfNodes)
{
    cout << " Graph started " << endl;
    int Begin, End;
    for (int i = 0; i < EdgeList.size(); ++i)
    {
        Begin = 1;//EdgeList[i][0];
        End = 2; //EdgeList[i][1];
        Edge edge_array[i] = Edge(Begin,End); 
    }

    const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);

    // declare a graph object
    Graph g(edge_array, edge_array + sizeof(edge_array) / sizeof(Edge), NumberOfNodes);
    cout << " Graph completed " << endl;
    return g;
}*/



int main(int argc, char** argv) {
    FullEvaluationCharCurve(true, 0.1);
    cout << "SAJT" <<endl;
    cin.get();
    wds = new Sensitivity(case_folder + case_name + ".inp");
    cout << "Optimisation started" << endl;
    ACO *ANTS = new ACO (5, wds->nodes.size(), 
                    0.5, 0.8, 800, 0.2, 2,
                    0);
    ANTS -> init();
    ANTS -> Search_Graph_Build_Up();
    ANTS -> printGRAPH ();
    ANTS -> printPHEROMONES ();
    ANTS -> optimize (10);
    ANTS -> printRESULTS ();
    cout << "optimalizáció vége" << endl;
    cin.get();
    MapsForPlot();
    cout << "OK" << endl;
    cin.get();
    cout << " Calculation started in the case of " << case_folder + case_name + ".inp" << " network..." << endl;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    /*for (int i = 0; i < wds->edges.size(); ++i)
    {
        cout << "ID: " << wds->edges.at(i)->getEdgeStringProperty("name") << " Flowrate: " << wds->edges.at(i)->getEdgeDoubleProperty("volumeFlowRate") << endl;
        cin.get();
    }*/
    //----------------------------------------------------//
    cout << " eddig elmegy" << endl;
    double MaximalPickableLength, MaximalPickableLengthLocal; //In meters
    double Diameter = 0.1, demandOrig;
    bool SolveByMultiplePipes = true;
    int NumberOfPossiblePipes = 5, id = 0;
    //----------------------------------------------------//
    FullEvaluation(true, Diameter);
    //----------------------------------------------------//
    for (int j = 0; j < wds->nodes.size(); ++j)
    {
        cout << "elem " << j << ", " << wds->nodes.at(j)->getName() << endl;
        if(wds->nodes.at(j)->getName() == "node_371")
        {
            id = j;
            cout << "megvan" << id << " " << wds->nodes.at(j)->getName() << "  demand: " << wds->nodes.at(id)->getProperty("demand") << endl;
        }
    }
    delete wds;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    //demandOrig = wds->nodes.at(id)->getProperty("demand");
    //cout << "demandOrig " << demandOrig << endl;
    //wds->nodes.at(id)->setProperty("demand",(1+i*0.5)*2);
    wds->calculateSensitivity("demand");
    //wds->fillNodalSensitivityForPlotNotNormalized();
    wds->saveResult("volumeFlowRate", "Pipe");//Psensitivity
    //cout << " Plot mehet " << wds->nodes.at(id)->getProperty("demand") << endl;
    cin.get();
    delete wds;
    bool PlotAll = true, again = true;  
    int NumberOfNodes = wds->nodes.size(), selector;
    vector<double> LocalSensitivity(wds->nodes.size()), ModifiedLocalSensitivity(wds->nodes.size());
    vector< vector<double> > NodeCoordinates (wds->nodes.size(), vector<double> (2)), Results, LocalSensitivityList;
    double SummarizedPipeLength, SummarizedNodalDemand;
    double l, AverageSensitivity, ModifiedAverageSensitivity, LocalSensitivityDifference, ModifiedLocalSensitivityDifference, AverageSensitivityDifference, AverageSensitivityDifferenceProc;
    double SummarizedLocalPipeLength, TravelTimeDifference, AverageTravelTime, ModifiedAverageTravelTime, AverageCapacity, ModifiedAverageCapacity, AverageCapacityDifference, AverageCapacityDifferenceProc, CharacteristicCurveDifference;
    vector<double> LocalResult;
    vector<double> InsertedPipeLength;
    vector<double> InsertedPipeDiameter;
    vector< vector<string> > NodeInput;
    vector<string> NodeLine;
    string Nev1, Nev2;
    //delete wds;
    //-----------------------------------------------------------//
    SummarizedPipeLength = CalculateSummarizedPipeLength();
    cout << "Sum pipe Length: " << SummarizedPipeLength << endl;
    SummarizedNodalDemand = CalculateSummarizedNodalDemand();
    cout << "Sum nodal Demand: " << SummarizedNodalDemand << endl;
    //cin.get();
    //-----------------------------------------------------------//
    //cout << " LocSensList! " << "Nev:" << wds->nodes.at(LocalSensitivityList[0][0])->getName() << " Nev2: " << wds->nodes.at(LocalSensitivityList[0][1])->getName() << endl;
    ofstream stream1;
    stream1.open("/mnt/d/Linux/Staci_Complete/Staci_Full_Version/staci3-master/Staci@HuT/Sensitivity/balf_2_Cso.dat");
    stream1 << "Node1" << " , " << "Node2"<< " , " << "AverageSensitivityDifferenceProc" << " , " <<  "PeakSensitivityDifferenceProc" << " , " << "Pipelength" << "," << "Price" << "\n";
    //----------------------------------------------------------//  //round(SummarizedPipeLength/10
    for (int i = 0; i < SummarizedPipeLength/10; ++i) 
    {
        //----------------------------------------------------//
        MaximalPickableLength = 10+i*10;//10 + i*10; //In meters
        Diameter = 0.1;
        //----------------------------------------------------//
        cout << "L: " << MaximalPickableLength << endl;
        //-----------------------------------------------------------//
        for (int j = 1; j <= NumberOfPossiblePipes; ++j)
        {
            cout <<"begincyclops" << endl;
            if (j==1)
            {
                cout << " Csak egy cso van bent .... " << endl;
                LocalSensitivityList = LocalSensitivityListCreator(MaximalPickableLength);
                wds = new Sensitivity(case_folder + case_name + ".inp");
                NodeLine.push_back(wds->nodes.at(LocalSensitivityList[0][0])->getName());
                NodeLine.push_back(wds->nodes.at(LocalSensitivityList[0][1])->getName());
                NodeInput.push_back(NodeLine);
                InsertedPipeLength.push_back(LocalSensitivityList[0][3]);
                InsertedPipeDiameter.push_back(Diameter);
                LocalResult = GenerateResults(NumberOfNodes, NodeInput, InsertedPipeLength, InsertedPipeDiameter, LocalSensitivityList, PlotAll);
                stream1 << wds->nodes.at(LocalSensitivityList[0][0])->getName() << "," << wds->nodes.at(LocalSensitivityList[0][1])->getName() << "," << LocalResult[0] << "," << LocalResult[1] << "," << LocalSensitivityList[0][3] << "," << LocalResult[2] << endl;
                LocalResult.clear();
                LocalSensitivityList.clear();
                NodeLine.clear();
                delete wds;   
            }
            else
            {
                cout << "megvan00" << endl;
                MaximalPickableLengthLocal = MaximalPickableLength/NumberOfPossiblePipes;
                NodeLine.clear();
                LocalSensitivityList = LocalSensitivityListCreator(MaximalPickableLengthLocal, NodeInput, InsertedPipeLength, InsertedPipeDiameter);
                selector = 0;
                again = true;
                while ( again == true )
                {
                    NodeLine.clear();
                    again = false;
                    for (int l = 0; l < NodeInput.size(); ++l)
                    {
                        NodeLine.push_back(wds->nodes.at(LocalSensitivityList[selector][0])->getName());
                        NodeLine.push_back(wds->nodes.at(LocalSensitivityList[selector][1])->getName());
                        if ((NodeLine[0] == NodeInput[l][0] && NodeLine[1] == NodeInput[l][1]) || (NodeLine[1] == NodeInput[l][0] && NodeLine[0] == NodeInput[l][1]))
                        {
                            selector += 1;
                            again = true;
                            NodeLine.clear();
                        }
                    }
                }
                cout << "megvan44 " << NodeLine.size() << endl;
                NodeInput.push_back(NodeLine);
                NodeLine.clear();
                InsertedPipeLength.push_back(LocalSensitivityList[selector][3]);
                InsertedPipeDiameter.push_back(Diameter);
                LocalSensitivityList.clear();
                cout << "megvan55" << endl;
                cout << "A csovek szama: " << j << "/" << NodeInput.size() << endl;
                cout << "Ez meg megvan" << endl;
            }
            cout << "Ez meg megvan" << endl;
        }
        //-----------------------------------------------------------//
        for (int j = 0; j < InsertedPipeLength.size(); ++j)
        {
            SummarizedLocalPipeLength += InsertedPipeLength[i];
        }
        LocalSensitivityList = LocalSensitivityListCreator(MaximalPickableLengthLocal, NodeInput, InsertedPipeLength, InsertedPipeDiameter);
        cout << "ide eljut" << endl;
        LocalResult = GenerateResults(NumberOfNodes, NodeInput, InsertedPipeLength, InsertedPipeDiameter, LocalSensitivityList, PlotAll);
        cout << "ide eljut" << endl;
        stream1 << wds->nodes.at(LocalSensitivityList[0][0])->getName() << "," << wds->nodes.at(LocalSensitivityList[0][1])->getName() << "," << LocalResult[0] << "," << LocalResult[1] << "," << LocalSensitivityList[0][3] << "," << LocalResult[2] << endl;
        NodeInput.clear();
        InsertedPipeLength.clear();
        InsertedPipeDiameter.clear();
        LocalResult.clear();
        LocalSensitivityList.clear();
        NodeLine.clear();
    }
    stream1.close();
    //-----------------------------------------------------------//
} 

/*    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if (wds->edges[i]->getEdgeStringProperty("type") == "Konstans nyomas" && first == true)
        {
            first = true;
            for (int j = 0; j < wds->nodes.size(); ++j)
            {
                if(wds->nodes.at(j)->getName() == wds->edges[i]->Get_Cspe_Nev())
                {
                    PressurePointsConnectingNodes.push_back(j);
                    first == false;
                }
            }
            PressurePoints.push_back(i);
            cout << endl << i;
        }
    }*/
/*
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            if ( l <= MaximalPickableLength)
            {
                LocalSensitivityDifference = abs(LocalSensitivity[i] - LocalSensitivity[j]);
                LocalSensitivityList.push_back(vector<double> ());
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(i);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(j);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(LocalSensitivityDifference);
                LocalSensitivityList[LocalSensitivityList.size()-1].push_back(l);
            }
        }
    }
    sort(LocalSensitivityList.begin(), LocalSensitivityList.end(), [](const std::vector< double >& a, const std::vector< double >& b){ return a[2] > b[2]; } );
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < LocalSensitivityList[i].size(); ++j)
        {
            cout << LocalSensitivityList[i][j] << " , ";
        }
        cout << endl;
        //cin.get();
    }
    //---------Original state sensitivity calculation------------//
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = wds->SM_row_sum_Pressures[i];
        AverageSensitivity += wds->SM_row_sum_Pressures[i];
        AverageTravelTime  += wds->nodes.at(i)->getProperty("tt");
    }
    AverageSensitivity = AverageSensitivity / NumberOfNodes;
    AverageTravelTime = AverageTravelTime / NumberOfNodes;
    AverageCapacity = CalculateAverageCapacity(critical_nodes, case_name, PlotAll);
    //delete wds;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    l = LocalSensitivityList[i][3];
//    wds->Set_debug_level(0);
    wds->edges.push_back(new Pipe("extra_cso", wds->nodes.at(LocalSensitivityList[i][0])->getName(), wds->nodes.at(LocalSensitivityList[i][1])->getName(), 1000, l, Diameter, -0.02, 0, 0, 0));
    //cout << "a ciklusba belep, itt tart: i[" << i << "] " << "j[" << j << "]" << endl;
    //wds->buildSystem();
    //wds->ini();
    wds->solveSystem();
    //---------Modified state sensitivity calculation------------//
    wds->Compute_Sensitivity_Matrix("demand", 0);
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        ModifiedLocalSensitivity[i] = wds->SM_row_sum_Pressures[i];
        ModifiedAverageSensitivity += wds->SM_row_sum_Pressures[i];
        ModifiedAverageTravelTime  += wds->nodes.at(i)->getProperty("tt");
    }
    ModifiedAverageTravelTime = ModifiedAverageTravelTime / NumberOfNodes;
    ModifiedAverageSensitivity = ModifiedAverageSensitivity / NumberOfNodes;
    LocalSensitivityDifference = abs(LocalSensitivity[LocalSensitivityList[i][0]] - LocalSensitivity[LocalSensitivityList[i][1]]);
    ModifiedLocalSensitivityDifference = abs(ModifiedLocalSensitivity[LocalSensitivityList[i][0]] - ModifiedLocalSensitivity[LocalSensitivityList[i][1]]);
    AverageSensitivityDifference = (1 - (ModifiedAverageSensitivity/AverageSensitivity));
    AverageSensitivityDifferenceProc = (1 - (ModifiedAverageSensitivity/AverageSensitivity))*100;
    //ModifiedAverageCapacity = CalculateAverageCapacity(l, Diameter, wds->nodes.at(LocalSensitivityList[i][0])->getName(), wds->nodes.at(LocalSensitivityList[i][1])->getName(), critical_nodes, case_name, PlotAll);
    AverageCapacityDifference = (1 - (AverageCapacity/ModifiedAverageCapacity));
    AverageCapacityDifferenceProc = (1 - (AverageCapacity/ModifiedAverageCapacity))*100;
    CharacteristicCurveDifference = CalculatSystemResistanceDifference(l, Diameter, wds->nodes.at(LocalSensitivityList[i][0])->getName(), wds->nodes.at(LocalSensitivityList[i][1])->getName(), PlotAll);
    TravelTimeDifference = (1-(ModifiedAverageTravelTime/AverageTravelTime))*100;
    cout << " TT : " << TravelTimeDifference << endl;
    //-----------------------------------------------------------//
    //stream1 << wds->nodes.at(LocalSensitivityList[i][0])->getName() << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getName() << " , " << AverageSensitivityDifferenceProc << " , " << AverageSensitivityDifference << " , " <<  LocalSensitivityDifference << " , " << ModifiedLocalSensitivityDifference << " , " << AverageCapacityDifferenceProc << " , " << l << " , " << CharacteristicCurveDifference << " , " << TravelTimeDifference << "\n";
    //cout << wds->nodes.at(LocalSensitivityList[i][0])->getName() << " , " << wds->nodes.at(LocalSensitivityList[i][1])->getName() << " , " << AverageSensitivityDifferenceProc << " , " << AverageSensitivityDifference << " , " <<  LocalSensitivityDifference << " , " << ModifiedLocalSensitivityDifference << " , " << AverageCapacityDifferenceProc  << " , " << l << " , " << CharacteristicCurveDifference << " , " << TravelTimeDifference << "\n";
    //-----------------------------------------------------------//
    wds->edges.pop_back();
    ////wds->buildSystem();
    */