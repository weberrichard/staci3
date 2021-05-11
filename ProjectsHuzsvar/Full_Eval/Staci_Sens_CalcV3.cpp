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

void CalculateCapacity(string FileName, string case_folder, string case_name, bool PlotPlease, bool from_file)
{
    Sensitivity *wds;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    int localindex = 0;
    wds->printLevel = 3;
    vector<double> LengthList;
    LengthList.push_back(2.);
    double dzeta = 0.5, AverageCapacity;
    vector<double> firewater;
    vector<int> ID;
    ifstream ifile;
    string local;
    ofstream wfile;
    vector<string> hydrant_nodes;
    if(from_file == true)
    {        
        //-----------------Read the selected nodes from txt---------------//
        ifile.open(FileName + ".txt");
        string line, temp, localstack;
        while(getline(ifile,line))
        {
            temp = line.substr(0,line.length());
            cout << line.substr(0,line.length()) << endl;
            for (int j = 0; j < wds->nodes.size(); ++j)
            {
               if (temp == wds->nodes[j]->getName())
               {
                   ID.push_back(j);
               }
            }
            hydrant_nodes.push_back(temp);
            cout << "hydrant_nodes " << hydrant_nodes.back() << endl;
        }
        ifile.close();
        //----------------------------------------------------------------//
    }
    else
    {
        for (int i = 0; i < wds->nodes.size(); i+=1)
        {
            hydrant_nodes.push_back(wds->nodes.at(i)->getName());
        }
    }
    if (PlotPlease == true)
    {
        wfile.open("firewater_"+case_name+".txt");
    }

    for (int i=0; i<hydrant_nodes.size()-1; i++)
    {
        //wds = new Sensitivity(case_folder + case_name + ".inp");
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            if (hydrant_nodes[i] == wds->nodes[j]->getName())
            {
                local = wds->nodes[j]->getName();
                if (local.substr(0, 4) == "pres" || local.substr(0, 4) == "pool" || local.substr(0, 4) == "PRES" || local.substr(0, 4) == "POOL")
                {
                    cout << "reservoir: " << local << endl;
                }
                else
                {
                    localindex = j;
                }
            }
        }     
        cout << "HN: " << localindex << " " << wds->nodes[localindex]->getName() << endl;
        cout << " megvan2 " << endl;
        vector<double> e,zeta;
        e.push_back(0);  e.push_back(100);
        cout << " megvan3 " << endl;
        zeta.push_back(0);
        zeta.push_back(dzeta);
        cout << " megvan3 " << endl;
        const double a = 1.;
        for (int k = 0; k < wds->edges.size(); ++k)
        {
            if(wds->edges[k]->getEdgeStringProperty("name") == "PRESS_Cap")
            {
                wds->edges[k]->setDoubleProperty("startHeight", wds->nodes[localindex]->getProperty("geodeticHeight"));
                wds->edges[k]->setDoubleProperty("head", 0.);
            }
        }
        for (int k = 0; k < wds->nodes.size(); ++k)
        {
            if(wds->nodes[k]->getName() == "PRESS_Cap")
            {
                wds->nodes[k]->setProperty("height",wds->nodes[localindex]->getProperty("geodeticHeight"));
            }
        }
        //wds->addNewEdgeElement(new Pipe("ExtraPipe", wds->nodes[2731]->getName(), wds->nodes[1003]->getName(), 1000., 505.584, 0.1, 0.02, 0.0, 0, 2));
        if(i == 0)
        {
            wds->addNewEdgeElement(new Pipe("OutflowPipe", wds->nodes[localindex]->getName(), "PRESS_Cap", 1000., 0.000000001, 0.1, 0.02, 0.0, 0, 2));
        }
        else
        {
            wds->modLastEdgeElement(wds->nodes[localindex]->getName(), "PRESS_Cap");
        }
        wds->calculateSensitivity("demand");
        cout << "FR: " << wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("volumeFlowRate")*3600 << endl;
        firewater.push_back(wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("volumeFlowRate")*3600);
        if (PlotPlease == true)
        {
            wfile << firewater[i] << '\n';
        }
    }
    for (int i = 0; i < firewater.size(); ++i)
    {
        AverageCapacity += firewater[i];
    }
    AverageCapacity /= firewater.size()     ;
    cout << " AC: " << AverageCapacity << endl;
    wfile.close();
}

void FullEvaluation(string case_folder, string case_name, bool PlotPlease, double Diameter)
{
    Sensitivity *wds;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    wds->calculateSensitivity("demand");
    int NumberOfNodes = wds->nodes.size();
    vector< vector<int> > EdgeList;
    vector<double> LengthList;
    vector<double> OriginalSourceDifferences;
    vector<double> ModifiedSourceDifferences;
    vector<int> Sources;
    double SqDiff;  
    bool mehet = true, PlotAll = false;
    auto start = high_resolution_clock::now();
    cout << " eddig eljutott" << endl;
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        EdgeList.push_back(vector<int> ());
        EdgeList[i].push_back(wds->edges.at(i)->getEdgeIntProperty("startNodeIndex"));
        EdgeList[i].push_back(wds->edges.at(i)->getEdgeIntProperty("endNodeIndex"));
        LengthList.push_back(wds->edges.at(i)->getEdgeDoubleProperty("length"));
    }
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if (wds->edges.at(i)->getEdgeIntProperty("numberNode") == 1)
        {
            Sources.push_back(i);
            cout << "catched" << endl;
        }
    }
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
    << duration.count() << " microseconds" << endl;
    MatrixXd SensMatrix2, SensMatrix, RowSumMatrix,RowSumMatrix2;
    vector<double> LocalSensitivity(wds->nodes.size()), ModifiedLocalSensitivity(wds->nodes.size()), NodalPressuresOriginal(wds->nodes.size());
    double OriginalPressureDifference, PressureDifference, PeakSensitivityDifferenceProc, PeakSensitivity = 0., ModifiedPeakSensitivity, l, AverageSensitivity = 0., ModifiedAverageSensitivity, LocalSensitivityDifference, ModifiedLocalSensitivityDifference, AverageSensitivityDifference, AverageSensitivityDifferenceProc;
    ofstream stream1;
    ofstream stream2;
    stream1.open("Full_Evaluation"+case_name+".txt");
    stream2.open("Full_Evaluation"+case_name+".dat");
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
    }
    AverageSensitivity = AverageSensitivity / NumberOfNodes;

    //---------------Topology Modification full evaluation-----------------//
    for (int i = 0; i < wds->nodes.size(); ++i)
    {
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            cout << "a ciklusba belep, itt tart: i[" << i << "] " << "j[" << j << "]" << endl;
            mehet = false;
            SqDiff = 0.;
            l = sqrt(pow((wds->nodes.at(i)->getProperty("xPosition") - wds->nodes.at(j)->getProperty("xPosition")), 2) + pow((wds->nodes.at(i)->getProperty("yPosition") - wds->nodes.at(j)->getProperty("yPosition")), 2));
            ModifiedAverageSensitivity = 0;
            ModifiedPeakSensitivity = 0;
            if(i != j)
            {
                if(i == 0 && j == 1)
                {
                    wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes.at(i)->getName(), wds->nodes.at(j)->getName(), 1000, l, Diameter, 0.02, 0.0,0,2));
                }
                else
                {
                    wds->modLastEdgeElement(wds->nodes.at(i)->getName(), wds->nodes.at(j)->getName(), l);
                }
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
                LocalSensitivityDifference = abs(LocalSensitivity[i] - LocalSensitivity[j]);
                ModifiedLocalSensitivityDifference = abs(ModifiedLocalSensitivity[i] - ModifiedLocalSensitivity[j]);
                AverageSensitivityDifference = abs(AverageSensitivity) - abs(ModifiedAverageSensitivity);
                AverageSensitivityDifferenceProc = AverageSensitivityDifference/abs(AverageSensitivity);
                PeakSensitivityDifferenceProc = (abs(PeakSensitivity) - abs(ModifiedPeakSensitivity))/abs(PeakSensitivity);
                cout << "checkpoint 1" << endl;
                if(PlotPlease == true)
                {
                    stream1 << i  << " , " << j << " , " << wds->nodes.at(i)->getName() << " , " << wds->nodes.at(j)->getName() << " , " << AverageSensitivityDifferenceProc;
                    stream1 << " , " << AverageSensitivityDifference << " , " << PeakSensitivityDifferenceProc << " , " <<  LocalSensitivityDifference << " , " ;
                    stream1 << ModifiedLocalSensitivityDifference << " , " << l << " , " << PressureDifference << " , " << OriginalPressureDifference << " , " ;
                    stream1 << wds->nodes.at(i)->getProperty("pressure") << " , " << wds->nodes.at(j)->getProperty("pressure") << " , " << NodalPressuresOriginal[i] << " , " ;
                    stream1 << NodalPressuresOriginal[j] << " , " << wds->nodes.at(i)->getProperty("height") << " , " << wds->nodes.at(j)->getProperty("height") << " , "; 
                    stream1 << wds->edges[wds->edges.size()-1]->getDoubleProperty("volumeFlowRate") << " , " << LocalSensitivity[i] << " , " << LocalSensitivity[j] << " , " << "\n";
                    if (AverageSensitivityDifference < 0)
                    {
                        stream2 << i  << " , " << j << " , " << wds->nodes.at(i)->getName() << " , " << wds->nodes.at(j)->getName() << " , " << AverageSensitivityDifferenceProc;
                        stream2 << " , " << AverageSensitivityDifference << " , " << PeakSensitivityDifferenceProc << " , " <<  LocalSensitivityDifference << " , " ;
                        stream2 << ModifiedLocalSensitivityDifference << " , " << l << " , " << PressureDifference << " , " << OriginalPressureDifference << " , " ;
                        stream2 << wds->nodes.at(i)->getProperty("pressure") << " , " << wds->nodes.at(j)->getProperty("pressure") << " , " << NodalPressuresOriginal[i] << " , " ;
                        stream2 << NodalPressuresOriginal[j] << " , " << wds->nodes.at(i)->getProperty("height") << " , " << wds->nodes.at(j)->getProperty("height") << " , "; 
                        stream2 << wds->edges[wds->edges.size()-1]->getDoubleProperty("volumeFlowRate") << " , " << LocalSensitivity[i] << " , " << LocalSensitivity[j] << " , " << "\n";
                    }
                }
                cout << "checkpoint 2" << endl;
            }
        }
    }
    stream1.close();
    stream2.close();
    cout << "ez megvan..." << endl;
}

void CalculateCapacity(string FileName, string case_folder, string case_name, bool PlotPlease, bool from_file, int StartNode, int EndNode, double Length)
{
    Sensitivity *wds;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    int localindex = 0;
    vector<double> LengthList;
    LengthList.push_back(2.);
    double dzeta = 0.5, AverageCapacity;
    vector<double> firewater;
    vector<int> ID;
    ifstream ifile;
    ofstream wfile;
    string local;
    vector<string> hydrant_nodes;
    if(from_file == true)
    {        
        //-----------------Read the selected nodes from txt---------------//
        ifile.open(FileName + ".txt");
        string line, temp, localstack;
        while(getline(ifile,line))
        {
            temp = line.substr(0,line.length());
            cout << line.substr(0,line.length()) << endl;
            for (int j = 0; j < wds->nodes.size(); ++j)
            {
               if (temp == wds->nodes[j]->getName())
               {
                   ID.push_back(j);
               }
            }
            hydrant_nodes.push_back(temp);
            cout << "hydrant_nodes " << hydrant_nodes.back() << endl;
        }
        ifile.close();
        //----------------------------------------------------------------//
    }
    else
    {
        for (int i = 0; i < wds->nodes.size()-1; i+=1)
        {
            hydrant_nodes.push_back(wds->nodes.at(i)->getName());
        }
    }
    if (PlotPlease == true)
    {
        wfile.open("firewater_"+case_name+".txt");
    }

    for (int i=0; i<hydrant_nodes.size(); i++)
    {
        //wds = new Sensitivity(case_folder + case_name + ".inp");
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
            if (hydrant_nodes[i] == wds->nodes[j]->getName())
            {
                local = wds->nodes[j]->getName();
                if (local.substr(0, 4) == "pres" || local.substr(0, 4) == "pool" || local.substr(0, 4) == "PRES" || local.substr(0, 4) == "POOL")
                {
                    cout << "reservoir" << local << endl;
                }
                else
                {
                    localindex = j;
                }
            }
        }      
        cout << "HN: " << localindex << " " << wds->nodes[localindex]->getName() << endl;
        cout << " megvan2 " << endl;
        vector<double> e,zeta;
        e.push_back(0);  e.push_back(100);
        cout << " megvan3 " << endl;
        zeta.push_back(0);
        zeta.push_back(dzeta);
        cout << " megvan3 " << endl;
        const double a = 1.;
        for (int k = 0; k < wds->edges.size(); ++k)
        {
            if(wds->edges[k]->getEdgeStringProperty("name") == "PRESS_Cap")
            {
                wds->edges[k]->setDoubleProperty("startHeight", wds->nodes[localindex]->getProperty("geodeticHeight"));
                wds->edges[k]->setDoubleProperty("head", 0.);
            }
        }
        for (int k = 0; k < wds->nodes.size(); ++k)
        {
            if(wds->nodes[k]->getName() == "PRESS_Cap")
            {
                wds->nodes[k]->setProperty("height",wds->nodes[localindex]->getProperty("geodeticHeight"));
            }
        }
        if(i == 0)
        {   
            wds->addNewEdgeElement(new Pipe("ExtraPipe", wds->nodes[StartNode]->getName(), wds->nodes[EndNode]->getName(), 1000., Length, 0.1, 0.02, 0.0, 0, 2));
            wds->addNewEdgeElement(new Pipe("OutflowPipe", wds->nodes[localindex]->getName(), "PRESS_Cap", 1000., 0.000000001, 0.1, 0.02, 0.0, 0, 2));
        }
        else
        {
            wds->modLastEdgeElement(wds->nodes[localindex]->getName(), "PRESS_Cap");
        }
        wds->calculateSensitivity("demand");
        cout << wds->edges[wds->edges.size()-1]->getEdgeStringProperty("name") << " FR: " << wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("volumeFlowRate")*3600 << endl;
        firewater.push_back(wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("volumeFlowRate")*3600);
        if (PlotPlease == true)
        {
            wfile << firewater[i] << '\n';
        }
    }
    for (int i = 0; i < firewater.size(); ++i)
    {
        AverageCapacity += firewater[i];
    }
    AverageCapacity /= firewater.size()     ;
    cout << " AC: " << AverageCapacity << endl;
    wfile.close();
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
    int StartNode, EndNode;
    double Length;
    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    cout << "[*]Generating full evaluation...." << endl;
    FullEvaluation(case_folder, case_name, true, 0.025);
    cout << "[*]Capacity distribution calculated sucessfully...." << endl;
}

