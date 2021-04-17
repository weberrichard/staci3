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
    //------------------------------------------------------------------------Staci init----------------------------------------------------------------//
    string case_folder = "../../Networks/Sopron/";
    string case_name = Network;
    //-----------------------------------------------------------------------Staci init End-------------------------------------------------------------//
    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    cout << "[*]Generating plots...." << endl;
    wds = new Vulnerability(case_folder + case_name + ".inp");
    wds->buildSegmentGraph();
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if(wds->edges[i]->typeCode == 1 || wds->edges[i]->typeCode == 0) // pipe, pipeCV
        {
            wds->edges.at(i)->setEdgeDoubleProperty("userOutput",wds->edges.at(i)->segment);
            cout << "segment" << wds->edges.at(i)->getEdgeDoubleProperty("userOutput") << endl; 
        }
    }
    wds->saveResult("userOutput", "Pipe");
    cout << "[*]Plot generation ended succesfully...." << endl;
}
