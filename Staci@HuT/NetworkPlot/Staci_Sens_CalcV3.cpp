#include "iostream"
#include <fstream>
#include <sstream>
#include "../../Sensitivity.h"
#include <cmath>
#include <algorithm>
#include <complex> 
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include <igraph.h>
#include <stdlib.h>
#include <utility>
#include </usr/include/eigen3/Eigen/Eigen>
#include </usr/include/eigen3/Eigen/Dense>
#include <stdio.h>
#include <chrono> 
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>
//-------------------------------------------------------------------------------------//
// for testing on real: /home/namrak/HTamas/halozatok/Sopron/
// VIZ-SOPTVR-J-55-input_mod
//-------------------------------------------------------------------------------------//

using namespace std;
using namespace Eigen;
using namespace chrono;
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//typedef pair<int, int> Edge;

Sensitivity *wds;

vector<vector<double>> parse2DCsvFile(string inputFileName)
{ 
    vector<vector<double> > data;
    //cout << "elso sor megvan..." << endl;
    ifstream inputFile;
    inputFile.open(inputFileName);
    //cout << "elso sor megvan..." << endl;
    int l = 0;
    while (inputFile)
    {
        //cout << "elso sor megvan2..." << endl;
        l++;
        string s;
        if (!getline(inputFile, s))
        {
            break;   
        }
        if (s[0] != ' ') 
        {
            istringstream ss(s);
            vector<double> record;
            while (ss) 
            {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try 
                {
                    //cout << "record..." << endl;
                    record.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
            data.push_back(record);
        }
        //cout << "elso sor megvan2..." << endl;
    }
    /*if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }*/
 
    return data;
}

int main(int argc, char *argv[])
{
    cout << "[*]Hydraulic solver started...." << endl;
    //cout << "ide belepett..." << endl;
    stringstream ss;
    string Network = argv[1];
    //------------------------------------------------------------------------Staci init----------------------------------------------------------------//
    string case_folder = "..\\Networks\\";
    string case_name = Network;
    string critical_nodes = Network + "_Nodelist";
    //-----------------------------------------------------------------------Staci init End-------------------------------------------------------------//
    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    cout << "[*]Generating plots...." << endl;
    MatrixXd SensMatrix, RowSumMatrix, DiagonalSensMatrix;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    //cout << " node1 name: " << wds->nodes[262]->getName() <<" node2 name: " << wds->nodes[2054]->getName() << endl;
    //cin.get();
    /*for (int k = 0; k < wds->edges.size(); ++k)
    {
        if(wds->edges[k]->getEdgeStringProperty("name") == "PRESS_Cap")
        {
            wds->edges[k]->setDoubleProperty("startHeight", wds->nodes[545]->getProperty("geodeticHeight"));
            wds->edges[k]->setDoubleProperty("head", 0.);
        }
    }
    for (int k = 0; k < wds->nodes.size(); ++k)
    {
        if(wds->nodes[k]->getName() == "PRESS_Cap")
        {
            wds->nodes[k]->setProperty("height",wds->nodes[545]->getProperty("geodeticHeight"));
        }
    }*/
   // wds->addNewEdgeElement(new Pipe("OutflowPipe", wds->nodes[545]->getName(), "PRESS_Cap", 1000., 0.000000001, 0.1, 0.02, 0.0, 0, 2));
    //wds->solveSystem();
    //wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes[620]->getName(), wds->nodes[1715]->getName(), 1000, 100.785, 0.1, 0.02, 0.0,0,2));
   // wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes[2132]->getName(), wds->nodes[6618]->getName(), 1000, 1111.97, 0.1, 0.02, 0.0,0,2));
    vector<double> LocalSensitivity(wds->nodes.size());
    wds->calculateSensitivity("demand");
    SensMatrix = wds->getPSensMatrix();
    DiagonalSensMatrix = SensMatrix.diagonal();
    RowSumMatrix = SensMatrix.rowwise().sum();
    double PeakSensitivity = 0.;
    //---------Original state sensitivity calculation------------//
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        LocalSensitivity[i] = RowSumMatrix(i);
        if (abs(RowSumMatrix(i)) > abs(PeakSensitivity))
        {
            PeakSensitivity = RowSumMatrix(i);
        }
        cout << DiagonalSensMatrix(i) << endl;
        //cout << "Av. Sens: " << AverageSensitivity << endl;
    }
    wds->fillNodalSensitivityForPlotNotNormalized();
    //cin.get();
    wds->saveResult("Psensitivity", "All");//Psensitivity
    //cout << " Plot mehet " << wds->nodes.at(id)->getProperty("demand") << endl;
    cout << "[*]Plot generation ended succesfully...." << endl;
}

