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
   string caseFolder = "../../Networks/Network_technion/";
   //string caseFolder = "../../Networks/Sopron/";

   vector<string> everyCase;
   everyCase.push_back("Network_3");
   everyCase.push_back("Network_9_ISO");
   everyCase.push_back("Network_13_ISO");

   int nCases = everyCase.size();
   cout << endl << "   CASES\n***********\n";
   for(int i=0; i<nCases; i++)
      cout << i+1 << "  " << everyCase[i] << endl;
   
   srand( (unsigned)time(NULL) );

   // for writing to files
   ofstream wFile;

   vector<vector<double> > everyLocalGamma(nCases);
   for(int i=0; i<nCases; i++)
   {
      printf("\n[*] %15s\n", everyCase[i].c_str());

      string caseName = everyCase[i];
      Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");
      cout << "LOAD OK " << endl;

      wds->calculateVulnerabilitySeries();
      cout << "RUN OK " << endl;

      double sumDem=0.;
      vector<double> betaSegment(wds->numberSegment,0.);
      for(int j=0; j<wds->nodes.size(); j++)
      {
         sumDem += wds->nodes[j]->demand;
         int seg = wds->nodes[j]->segment;
         betaSegment[seg] += wds->nodes[j]->demand;
      }
      for(int j=0; j<wds->numberSegment; j++)
      {
         betaSegment[j] /= sumDem;
      }

      // writing to file for further analysis
      wFile.open("Network Data/" + caseName + "/beta_segment.txt");
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << betaSegment[j] << '\n';
      }
      wFile.close();

      wFile.open("Network Data/" + caseName + "/beta_full.txt");
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << wds->relativeDemandLoss[j] << '\n';
      }
      wFile.close();

      wFile.open("Network Data/" + caseName + "/alfa_rel.txt");
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << wds->relativePipeLength[j] << '\n';
      }
      wFile.close();
   }
   
   cout << endl << endl;
   return 0;
}