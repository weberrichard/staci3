#include <fstream>
#include <sstream>
#include <string>
#include "../../Leakage.h"


using namespace std;
using namespace Eigen;
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//typedef pair<int, int> Edge;

Leakage *wds;


int main(int argc, char *argv[])
{
    string caseFolder = "Networks/";

    double nominalISODiameter = 0.2;
    double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;
    double setting = 10.;
    double minorLoss = 0.02;
    string caseName = argv[1];
    string inFileName = argv[2];

    ofstream wFile(caseName+"_edgelist.txt");
    ofstream wFile2(caseName+"_valvelist.txt");
    vector<double> settings, minorLosses;
    cout << "[*]Network: " << caseFolder << caseName << ".inp" << endl;
    wds = new Leakage(caseFolder +caseName + ".inp");
    if(inFileName == "EdgeList")
    {
        for (int i = 0; i < wds->edges.size(); ++i)
        {
            if (wds->edges.at(i)->getEdgeStringProperty("type") != "Pool" || wds->edges.at(i)->getEdgeStringProperty("type") == "Pressure")
            {
                wFile << wds->edges.at(i)->getEdgeIntProperty("startNodeIndex") << "," << wds->edges.at(i)->getEdgeIntProperty("endNodeIndex") << '\n';
            }
            if (wds->edges.at(i)->getEdgeStringProperty("type") == "ValveISO")
            {
                wFile2 << wds->edges.at(i)->getEdgeIntProperty("startNodeIndex") << "," << wds->edges.at(i)->getEdgeIntProperty("endNodeIndex") << '\n';
            }
        }
    }
}
