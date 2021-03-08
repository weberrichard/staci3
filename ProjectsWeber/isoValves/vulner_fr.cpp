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

double randomDouble(double a, double b)
{
   double random = ((double) rand()) / (double) RAND_MAX;
   double diff = b - a;
   double r = random * diff;
   return a + r;
}

int main(int argc, char* argv[])
{
   // Name of containing folder of staci file
   string caseFolder = "../../Networks/Sopron/";

   vector<string> everyCase;
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
   makeDirectory("Network Data");

   // for writing to files
   ofstream wFile;

   for(int i=0; i<nCases; i++)
   {
      printf("\n[*] %15s\n", everyCase[i].c_str());
      string caseName = everyCase[i];

      Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");
      wds->buildSegmentGraph();

      vector<double> failureRate(wds->numberSegment);
      double norm=0.;
      for(int j=0; j<wds->numberSegment; j++)
      {
         double r = randomDouble(0.,1.);
         failureRate[j] = wds->relativePipeLength[j]*r;
         norm += failureRate[j];
      }
      for(int j=0; j<wds->numberSegment; j++)
      {
         failureRate[j] /= norm;
      }

      wds->calculateVulnerability(failureRate);

     // saving original vulnerability
      vector<double> vulner;
      for(unsigned int j=0; j<wds->getNumberSegment(); j++)
      {
         vulner.push_back(wds->localGamma[j]);
      }

      // writing to file for further analysis
      wFile.open("Network Data/" + caseName + "/vulner_orig_rand.txt");
      for(int j=0; j<vulner.size(); j++)
      {
         wFile << vulner[j] << '\n';
      }
      wFile.close();
      
   }

   cout << endl << endl;
   return 0;
}

