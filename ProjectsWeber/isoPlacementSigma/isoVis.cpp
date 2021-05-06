#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../Vulnerability.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
   // Name of containing folder of staci file
   string caseFolder = "../../Networks/Sopron/";
   ofstream wFile("results_vis.txt");

   double nominalISODiameter = 0.2;
   double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;

   if(argc < 2)
   {
      cout << " !ERROR! At least two input arguments is required" << endl;
      exit(-1);
   }

   string caseName = argv[1];

   Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");
   wds->buildSegmentGraph();
   wds->calculateVulnerability();
   double sigma = standardDeviation(wds->localGamma);
   wFile << wds->globalGamma << ", " << sigma << ", ";

   // deleting every existing 
   vector<int> ISOValvesToDelete;
   for(unsigned int i=0; i<wds->valveISOIndex.size(); i++)
   {
      ISOValvesToDelete.push_back(wds->valveISOIndex[i]);
   }
   wds->deleteISOValves(ISOValvesToDelete);

   // loading iso valve positions
   string inFileName = argv[2];
   vector<string> fileData = readVectorString(inFileName);
   vector<int> addISOPipe;
   vector<bool> isStart;
   vector<string> addISOName;
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

      addISOPipe.push_back(stoi(sv[0]));
      isStart.push_back(stoi(sv[1]));
      addISOName.push_back("ISO_" + to_string(i));
      refCros.push_back(refA);
   }

   // adding the new iso valves
   wds->addNewISOValves(addISOName, addISOPipe, isStart, 1000., refCros, 0.);
   wds->buildSegmentGraph();
   wds->calculateVulnerability();
   sigma = standardDeviation(wds->localGamma);
   wFile << wds->globalGamma << ", " << sigma << '\n';

   // saving the new optimized system
   string folder = "Sopron_sigma_opt";
   makeDirectory(folder);
   wds->saveSystem(folder + '/' + caseName + "_sigma_opt.inp");

   wFile.close();
   //cout << endl << endl;
   return 0;
}

