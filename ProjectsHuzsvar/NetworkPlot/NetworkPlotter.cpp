#include "iostream"
#include <fstream>
#include <sstream>
#include "../../Vulnerability.h"
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
    double nominalISODiameter = 0.2;
    double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;
    //------------------------------------------------------------------------Staci init----------------------------------------------------------------//
    string case_folder = "../../Networks/Sopron/";
    string case_name = Network;
    //-----------------------------------------------------------------------Staci init End-------------------------------------------------------------//
    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    cout << "[*]Generating plots...." << endl;
    wds = new Vulnerability(case_folder + case_name + ".inp");
    vector<int> ISOValvesToDelete;
    for(unsigned int i=0; i<wds->valveISOIndex.size(); i++)
    {
     ISOValvesToDelete.push_back(wds->valveISOIndex[i]);
    }
    wds->deleteISOValves(ISOValvesToDelete);
    // todo string helyett index

    // loading iso valve positions
    string inFileName = argv[3];
    vector<string> fileData = readVectorString(inFileName);
    vector<int> addISOPipe;
    vector<bool> isStart;
    vector<string> addISOName;
    vector<double> refCros;
    for(int i=0; i<fileData.size(); i++)
    {
     stringstream ss(fileData[i]);
     vector<string> sv;
     cout << fileData[i] << endl;
     while(ss.good())
     {
        string substr;
        getline(ss,substr,',');
        sv.push_back(substr);
     }
     cout << sv[0] << "," << sv[1] << endl;
     addISOPipe.push_back(stoi(sv[0]));
     isStart.push_back(stoi(sv[1]));
     addISOName.push_back("ISO_" + to_string(i));
     refCros.push_back(refA);
    }

    // adding the new iso valves
    wds->addNewISOValves(addISOName, addISOPipe, isStart, 1000., refCros, 0.);
    wds->buildSegmentGraph();
    //wds->buildSegmentGraph();
    wds->solveSystem();
    double sum_PL = 0.;
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if(wds->edges[i]->typeCode == 1 || wds->edges[i]->typeCode == 0) // pipe, pipeCV
        {
            sum_PL += wds->edges.at(i)->getDoubleProperty("length");
            wds->edges.at(i)->setEdgeDoubleProperty("userOutput",wds->edges.at(i)->segment);
            cout << "segment" << wds->edges.at(i)->getEdgeDoubleProperty("userOutput") << endl; 
        }
    }
    int numberofvalves = 0;
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if(wds->edges[i]->typeCode == 9) // pipe, pipeCV
        {
            numberofvalves += 1;
            //cout << "segment" << wds->edges.at(i)->getEdgeDoubleProperty("userOutput") << endl; 
        }
    }
    double consumption = 0.;
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        consumption += wds->nodes.at(i)->getProperty("consumption");
    }
    cout << "number of nodes: " << wds->nodes.size() << endl;
    cout << "Sum. length: " << sum_PL << endl;
    cout << "Sum. consumption: " << consumption << endl;
    cout << "Number of valves: " << numberofvalves << endl;
    wds->saveResult("userOutput", "Pipe");
    cout << "[*]Plot generation ended succesfully...." << endl;
}
