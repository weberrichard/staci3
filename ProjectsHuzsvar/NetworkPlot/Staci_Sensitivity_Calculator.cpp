#include "iostream"
#include "../../Sensitivity.h"
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
    ofstream stream1;
    vector<double> LengthList;
    vector< vector<int> > InputEdges;
    int localCounter = 0, counter = 0, NumberOfNodes = 0;
    double AverageSensitivity = 0., PeakSensitivity = 0.;
    vector<vector<double>> data = parse2DCsvFile("ListOfNewEdges.txt");
    //cout << "ez még megvan.222...." << endl;
    for (auto l : data)
    {
        localCounter = 0;
        InputEdges.push_back(vector<int> ());
        for (auto x : l)
        {
            ss << x;
            if(localCounter < 2)
            {
                InputEdges[counter].push_back(stoi(ss.str()));
            }
            else
            {
                LengthList.push_back(stof(ss.str()));
            }
            localCounter += 1;
            ss.str("");
        }
        counter += 1;
    }    
    auto start = high_resolution_clock::now();
    double Diameter = 0.1;
    if(argc >= 4 && string(argv[2]) == "D")
    {
        Diameter = stod(argv[3]);
        cout << "[*]Diameter modified: " << argv[3] << endl;
    }
    if(argc >= 6 && string(argv[4]) == "BS")
    {
        AverageSensitivity = stod(argv[5]);
    }
    else
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
        wds->calculateSensitivity("demand");
        NumberOfNodes = wds->nodes.size();
        MatrixXd SensMatrix, RowSumMatrix;
        vector<double> LocalSensitivity(wds->nodes.size());
        
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
        }
        AverageSensitivity = AverageSensitivity / NumberOfNodes;
        cout << "[*]AverageSensitivity (Original network): " << AverageSensitivity << endl;
    }
    stream1.open("ListOfNewEdges.txt");
    //---------------Topology Modification full evaluation-----------------//
    for (int i = 0; i < InputEdges.size(); ++i)
    {
        delete wds;
        wds = new Sensitivity(case_folder + case_name + ".inp");
        NumberOfNodes = wds->nodes.size();
        MatrixXd SensMatrix2, RowSumMatrix2;
        vector<double> ModifiedLocalSensitivity(wds->nodes.size());
        double PeakSensitivityDifferenceProc, ModifiedPeakSensitivity, ModifiedAverageSensitivity, AverageSensitivityDifference, AverageSensitivityDifferenceProc;
        cout << "[*]Insert new pipeline: " << wds->nodes[InputEdges[i][0]]->getName() << " and " << wds->nodes[InputEdges[i][1]]->getName() << " with the length of : " << LengthList[i] << endl;
        wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes[InputEdges[i][0]]->getName(), wds->nodes[InputEdges[i][1]]->getName(), 1000, LengthList[i], Diameter, 0.02, 0.0,0,2));
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
        ModifiedAverageSensitivity = ModifiedAverageSensitivity / NumberOfNodes;
        cout << "[*]AverageSensitivity (Modified network): " << ModifiedAverageSensitivity << endl;
        AverageSensitivityDifference = abs(AverageSensitivity) - abs(ModifiedAverageSensitivity);
        AverageSensitivityDifferenceProc = (AverageSensitivityDifference/abs(AverageSensitivity))*100;
        PeakSensitivityDifferenceProc = (abs(PeakSensitivity) - abs(ModifiedPeakSensitivity))/abs(PeakSensitivity);
        cout << "[*]Simulation results: " << endl << "   Start node: " << InputEdges[i][0] << endl << "   End node: " 
        << InputEdges[i][1] << endl << "   Inserted pipelength [m]: " << LengthList[i] << endl << "   Avg. Sens. Difference [%]: " << AverageSensitivityDifferenceProc << "\n";
        stream1 << InputEdges[i][0] << "," << InputEdges[i][1] << "," << LengthList[i] << "," << AverageSensitivityDifferenceProc << "\n";
    }
    stream1.close();
    //cout << "ez megvan..." << endl;
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "[*]Runtime: " << duration.count()/1000/1000 << " [s]" << endl;
}

