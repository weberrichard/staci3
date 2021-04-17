#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>

#include "../../DMA.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
   ofstream wFile("results.txt");
   double nominalISODiameter = 0.2;
   double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;

   // Name of containing folder of staci file
   //string case_folder = "../../Networks/Sopron/";
   //string case_folder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";

   string caseName = argv[1];
   string runType = argv[2];

   //cout << endl << "Case: " << caseName << endl;

   DMA *wds;
   wds = new DMA(caseName);

   // simple stuff for input data
   if(runType == "Min") // min index of pipes
   {
      double n = wds->poolIndex.size() + wds->presIndex.size();
      wFile << n << '\n';
   }
   else if(runType == "Max") // max index of pipes
   {
      double n = wds->poolIndex.size() + wds->presIndex.size() + wds->pipeIndex.size();
      wFile << n << '\n';
   }
   else if(runType == "ID")
   {
      string inFileName = argv[3];
      vector<string> fileData = readVectorString(inFileName);
      for(int i=0; i<fileData.size(); i++)
      {
         int index = stoi(fileData[i]);
         string ID = wds->edges[index]->name;
         wFile << ID << '\n';
      }
   }

   if(runType == "nDMA" || runType == "demand")
   {
      // adding flow meters
      string inFileName = argv[3];
      vector<string> fileData = readVectorString(inFileName);
      vector<int> addFMPipe;
      vector<bool> isStart;
      vector<string> addFMName;
      vector<double> refCros;
      for(int i=0; i<fileData.size(); i++)
      {
         stringstream ss(fileData[i]);
         vector<string> sv;

         while(ss.good())
         {
            string substr;
            getline(ss,substr,',');
            sv.push_back(substr);
         }
         addFMPipe.push_back(stoi(sv[0]));
         isStart.push_back(1);
         addFMName.push_back("FM_" + to_string(i));
         refCros.push_back(refA);
      }
      wds->addNewFlowMeter(addFMName, addFMPipe, isStart, 1000., refCros, 0.);

      wds->determineDMAZones();

      if(runType == "nDMA")
      {
         wFile << wds->numberDMAZones;
      }
      else if(runType == "demand")
      {
         // calculating zone demands
         vector<double> zoneDemand(wds->numberDMAZones,0.);
         for(int i=0; i<wds->numberNodes; i++)
         {
            int z = wds->nodes[i]->DMAZone;
            zoneDemand[z] += wds->nodes[i]->demand;
         }

         // write to file
         for(int i=0; i<wds->numberDMAZones; i++)
         {
            wFile << zoneDemand[i] << '\n';
         }
      }
   }

   //cout << endl << endl;
   return 0;
}