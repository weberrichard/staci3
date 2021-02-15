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

void main(int argc, char* argv[])
{
   // Name of containing folder of staci file
   string caseFolder = "../Networks/";

   vector<string> everyCase;
   //everyCase.push_back("nagycenk");
   everyCase.push_back("villasor");
   everyCase.push_back("ferto");
   everyCase.push_back("sanchegy");
   everyCase.push_back("buk");
   everyCase.push_back("lovo");
   everyCase.push_back("nagycenk");
   everyCase.push_back("vashegy");
   everyCase.push_back("varis");
   everyCase.push_back("becsidomb");
   everyCase.push_back("tomalom");
   everyCase.push_back("szakov");
   everyCase.push_back("kohegy");
   everyCase.push_back("harka");
   everyCase.push_back("pozsonyiut");
   everyCase.push_back("sopronkovesd");
   everyCase.push_back("dudlesz");
   everyCase.push_back("ivan");
   everyCase.push_back("agyagosszergeny");
   everyCase.push_back("kofejto");
   everyCase.push_back("simasag");
   everyCase.push_back("acsad");
   everyCase.push_back("csaford");
   everyCase.push_back("nagylozs");
   everyCase.push_back("balf");
   everyCase.push_back("csapod");
   everyCase.push_back("und");
   everyCase.push_back("rojtokmuzsaj");
   everyCase.push_back("brennberg");
   everyCase.push_back("pusztacsalad");
   everyCase.push_back("kutyahegy");
   everyCase.push_back("nyarliget");
   everyCase.push_back("meszlen");
   everyCase.push_back("fertoujlak");
   everyCase.push_back("gorbehalom");
   everyCase.push_back("tozeggyarmajor");
   everyCase.push_back("ebergoc");
   everyCase.push_back("csillahegy");
   everyCase.push_back("jerevan");
   everyCase.push_back("gloriette");
   everyCase.push_back("alomhegy");
   everyCase.push_back("ohermes");
   everyCase.push_back("ujhermes");

   int nCases = everyCase.size();
   cout << endl << "   CASES\n***********\n";
   for(int i=0; i<nCases; i++)
      cout << i+1 << "  " << everyCase[i] << endl;
   
   srand( (unsigned)time(NULL) );

   vector<vector<double> > everyLocalGamma(nCases);
   for(int i=0; i<nCases; i++)
   {
      printf("\n[*] %15s\n", everyCase[i].c_str());

      string caseName = everyCase[i];
      Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");

      //double pat = 1.0;
      //for(int j=0;j<wds->nodes.size(); j++)
      //{
      //   double d = wds->nodes[j]->getProperty("demand");
      //   double s = (rand()%1000)/1000.-0.5;
      //   wds->nodes[j]->setProperty("demand",pat*d*s);
      //}

      wds->calculateVulnerability();

      everyLocalGamma[i].resize(wds->getNumberSegment());
      for(int j=0; j<wds->getNumberSegment(); j++)
      {
         everyLocalGamma[i][j] = wds->localGamma[j];
      }

      for(int j=0; j<wds->nodes.size(); j++)
      {
         int s = wds->nodes[j]->getProperty("segment");
         double v = wds->localGamma[s];
         wds->nodes[j]->setProperty("userOutput",v);
      }
      for(int j=0; j<wds->edges.size(); j++)
      {
         if(wds->edges[j]->typeCode == 1)
         {
            int s = wds->edges[j]->getEdgeIntProperty("segment");
            double v = wds->localGamma[s];
            wds->edges[j]->setEdgeDoubleProperty("userOutput",v);
         }
      }
      wds->saveResult("userOutput","Node");
      wds->saveResult("userOutput","Pipe");
   }

   ofstream wFile;
   wFile.open("localGamma.txt");
   for(int i=0; i<nCases; i++)
   {
      wFile << everyLocalGamma[i][0];
      for(int j=0; j<everyLocalGamma[i].size(); j++)
      {
        wFile << ", " << everyLocalGamma[i][j];
      }
      wFile << ";" << endl;
   }
   wFile.close();

   cout << endl << endl;
   return 0;
}

void CalculateCapacity(string FileName, string case_folder, string case_name, bool PlotPlease, bool from_file)
{
    cout << "belép" << endl;
    wds = new Sensitivity(case_folder + case_name + ".inp");
    int localindex = 0;
    vector<double> LengthList;
    LengthList.push_back(2.);
    double dzeta = 0.5, AverageCapacity;
    vector<double> firewater;
    vector<int> ID;
    ifstream ifile;
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
        delete wds;
        //----------------------------------------------------------------//
    }
    else
    {
        for (int i = 0; i < wds->nodes.size(); i+=100)
        {
            hydrant_nodes.push_back(wds->nodes.at(i)->getName());
        }
    }
    if (PlotPlease == true)
    {
        wfile.open("firewater_extension_"+case_name+".txt");
    }

    for (int i=0; i<hydrant_nodes.size(); i++)
    {
        wds = new Sensitivity(case_folder + case_name + ".inp");
        for (int j = 0; j < wds->nodes.size(); ++j)
        {
           if (hydrant_nodes[i] == wds->nodes[j]->getName())
           {
               localindex = j;
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
                wds->edges[k]->setDoubleProperty("head",wds->nodes[localindex]->getProperty("height"));
            }
        }
        //wds->addNewEdgeElement(new Pipe("extra_cso", wds->nodes[1293]->getName(), wds->nodes[1110]->getName(), 1000, 749.51, 0.1, 0.02, 0.0,0,2));
        //wds->addNewEdgeElement(new Pipe("OutflowPipe", wds->nodes[localindex]->getName(), "PRESS_Cap", 1000., 0.000000001, 0.1, 0.02, 0.0, 0, 2));
        cout << " megvan5 " << endl;
        wds->calculateSensitivity("demand");
        cout << "FR: " << wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("volumeFlowRate")*3600 << endl;
        firewater.push_back(wds->edges[wds->edges.size()-1]->getEdgeDoubleProperty("volumeFlowRate")*3600);
        if (PlotPlease == true)
        {
            wfile << firewater[i] << '\n';
        }
        delete wds;
    }
    for (int i = 0; i < firewater.size(); ++i)
    {
        AverageCapacity += firewater[i];
    }
    AverageCapacity /= firewater.size();
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
    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    for (int i = 0; i < wds->edges.size(); ++i)
    {
        if(wds->edges.at(i)->getEdgeStringProperty("type") == "Valve")
        {   
            wds->edges.erase(i);
            wds->addNewPipe(wds->edges.at(i)->getEdgeStringProperty("name"),wds->edges.at(i)->getEdgeStringProperty("startNodeName"),wds->edges.at(i)->getEdgeStringProperty("endNodeName"), 1000.,1.,0.1,0.02,0.,false,1,0.);
        }
    }
    cout << "[*]Generating capacity distribution...." << endl;
    CalculateCapacity("Varis_Hydrants", case_folder, case_name, true, false);
    cout << "[*]Capacity distribution calculated sucessfully...." << endl;
}

