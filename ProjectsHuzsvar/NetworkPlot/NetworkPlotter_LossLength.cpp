#include "iostream"
#include <fstream>
#include <sstream>
#include "../../Vulnerability.h"
#include "../../Statistic.h"
#include <cmath>
#include <algorithm>
#include <complex> 
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
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

Vulnerability *wds;

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
    stringstream ss;
    string Network = argv[1];
    string ID1 = argv[2];
    string ID2 = argv[3];
    //------------------------------------------------------------------------Staci init----------------------------------------------------------------//
    string case_folder = "../../Networks/Sopron/";
    string case_name = Network;
    //-----------------------------------------------------------------------Staci init End-------------------------------------------------------------//
    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    cout << "[*]Generating plots...." << endl;
    wds = new Vulnerability(case_folder + case_name + ".inp");
    wds->solveSystem();
    double EdgeLoss = 0., EdgeLossMax = 0.;
    cout << "node_1: " << wds->nodes.at(stoi(ID1))->name << " node_2: " << wds->nodes.at(stoi(ID2))->name << endl;
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        wds->nodes.at(i)->setProperty("userOutput",0.);
        if(i == stoi(ID1) || i == stoi(ID2)) // pipe, pipeCV
        {
            wds->nodes.at(i)->setProperty("userOutput",99.);
            cout << "check: " << i << endl;
        }
    }
    cin.get();
    int EdgeLossID = 0;
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if(wds->edges[i]->typeCode == 1 || wds->edges[i]->typeCode == 0) // pipe, pipeCV
        {
            EdgeLoss = wds->edges.at(i)->getDoubleProperty("volumeFlowRate")*wds->edges.at(i)->getDoubleProperty("volumeFlowRate")*wds->edges.at(i)->getDoubleProperty("pipeConst");
            if(EdgeLoss > EdgeLossMax)
                EdgeLossMax = EdgeLoss;
                EdgeLossID = i;
            
        }
    }
    cout << "Edge loss max: " << EdgeLossMax << " , " << wds->edges.at(EdgeLossID)->getEdgeIntProperty("startNodeIndex") << " , " << wds->edges.at(EdgeLossID)->getEdgeIntProperty("endNodeIndex") << " , " <<wds->edges.at(EdgeLossID)->getDoubleProperty("length")<< endl; 
    cout << EdgeLossMax << endl;
    double Input = 0.;
    vector<int> IndexList;
    vector<double> LossLength; 
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if(wds->edges[i]->typeCode == 1 || wds->edges[i]->typeCode == 0) // pipe, pipeCV
        {
            EdgeLoss = wds->edges.at(i)->getDoubleProperty("volumeFlowRate")*wds->edges.at(i)->getDoubleProperty("volumeFlowRate")*wds->edges.at(i)->getDoubleProperty("pipeConst");
            //cout << "Edge loss: " << EdgeLoss << endl; 
            LossLength.push_back(EdgeLoss/wds->edges.at(i)->getDoubleProperty("length"));
            Input = (EdgeLoss)/EdgeLossMax;///wds->edges.at(i)->getDoubleProperty("length")
            //Input = 2;
            if(Input < 0.015q)
            {
                wds->edges.at(i)->setEdgeDoubleProperty("userOutput",Input);
                IndexList.push_back(i);
            }
        }
    }
    /*quickSort(LossLength, IndexList, 0, wds->edges.size());
    //wds->solveSystem();
    for (int i = 0; i < IndexList.size(); ++i)
    {
        if(wds->edges[i]->typeCode == 1 || wds->edges[i]->typeCode == 0) // pipe, pipeCV
        {
            EdgeLoss = wds->edges.at(IndexList[i])->getDoubleProperty("volumeFlowRate")*wds->edges.at(IndexList[i])->getDoubleProperty("volumeFlowRate")*wds->edges.at(IndexList[i])->getDoubleProperty("pipeConst");
            
            wds->edges.at(IndexList[i])->setEdgeDoubleProperty("userOutput",i);
        }
        cout << "Edge loss: " << LossLength[i] << " , " << IndexList[i] << endl; 
    }*/
    wds->saveResult("userOutput", "All");
    cout << "[*]Plot generation ended succesfully...." << endl;
}
