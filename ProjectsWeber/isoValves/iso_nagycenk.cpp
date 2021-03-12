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

   string caseName = "nagycenk_mod";

   srand( (unsigned)time(NULL) );
   makeDirectory("Network Data");

   // for writing to files
   ofstream wFile;

   makeDirectory(("Network Data/" + caseName).c_str());

   cout << "\n\n Original vulnerability calculation";
   cout << "\n------------------------------------";
   Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");

   vector<string> addISOName;
   vector<string> addISOPipe;
   vector<bool> isStart;
   vector<double> refCros;

   double nominalISODiameter = 0.2;
   double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;

   /*addISOName.push_back("iso_1"); 
   addISOPipe.push_back("PIPE_2669853");
   isStart.push_back(true);
   refCros.push_back(refA);

   addISOName.push_back("iso_2");
   addISOPipe.push_back("PIPE_540915");
   isStart.push_back(true);
   refCros.push_back(refA);

   addISOName.push_back("iso_3");
   addISOPipe.push_back("PIPE_ext_1");
   isStart.push_back(false);
   refCros.push_back(refA);

   addISOName.push_back("iso_4");
   addISOPipe.push_back("PIPE_ext_2");
   isStart.push_back(true);
   refCros.push_back(refA);*/

   wds->addNewISOValves(addISOName, addISOPipe, isStart, 1000., refCros, 0.);
   wds->buildSegmentGraph();
   wds->calculateVulnerability();

   // saving original vulnerability
   vector<double> vulner;
   for(unsigned int j=0; j<wds->getNumberSegment(); j++)
   {
      vulner.push_back(wds->localGamma[j]);
   }

   // saving relative demand loss
   for(unsigned int j=0; j<wds->edges.size(); j++)
   {  
      int seg = wds->edges[j]->segment;
      double vul = wds->localGamma[seg];
      wds->edges[j]->setEdgeDoubleProperty("userOutput",vul);
   }
   wds->saveResult("userOutput","Pipe");

   cout << endl << "Gamma: " << wds->globalGamma << endl;

   cout << endl << endl;
   return 0;
}

