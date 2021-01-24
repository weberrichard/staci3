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
    cout << "[*]Generating search space started...." << endl;
    ofstream stream1, stream2;
    MatrixXd SensMatrix, RowSumMatrix;
    double DS = 0., l = 0., SumLength = 0., SumDemand = 0., limit = 0.;
    stream2.open("SummarizedParameters_"+case_name+".txt");
    wds = new Sensitivity(case_folder + case_name + ".inp");
    wds->initialization();
    wds->calculateSensitivity("demand");
    vector<double> LocalSensitivity(wds->nodes.size());
    SensMatrix = wds->getPSensMatrix();
    RowSumMatrix = SensMatrix.rowwise().sum();
    string StartNodeName, EndNodeName;
    stream1.open("SearchSpace_"+case_name+".txt");
    if(argc < 3)
    {
        limit = 5.;
    }
    else
    {
        limit = stod(argv[2]);
    }
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        SumDemand += wds->nodes.at(i)->getProperty("demand");
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            if(i != j && l <= limit)
            {
                DS = RowSumMatrix(i) - RowSumMatrix(j);
                StartNodeName = wds->nodes[i]->getName();
                EndNodeName = wds->nodes[j]->getName();
                stream1 << i << "," << j << "," << StartNodeName << "," << EndNodeName << "," << DS << "," << l << "\n";
            }
        }
    }
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        //cout << wds->edges[i]->getEdgeStringProperty("type") << endl;
        if (wds->edges[i]->getEdgeStringProperty("type") == "Pipe")
        {
            SumLength += wds->edges[i]->getDoubleProperty("length");
        }
    }
    stream2 << "Summarized Pipe Length: " << SumLength << "\n" << "Summarized Demand: " << SumDemand << "\n";
    stream1.close();
    stream2.close();
    cout << "[*]Generating search space completed...." << endl;
}

