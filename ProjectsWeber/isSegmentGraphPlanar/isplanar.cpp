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

   vector<vector<double> > everyLocalGamma(nCases);
   for(int i=0; i<nCases; i++)
   {
      printf("\n[*] %15s\n", everyCase[i].c_str());

      string caseName = everyCase[i];
      Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");

      wds->calculateVulnerability();   

      vector<int> sev = wds->getSegmentEdgeVector();
      vector<double> bi = wds->getRelativeDemandLoss();

      vector<double> di(wds->getNumberSegment(),0.);
      double sumD=0.;
      for(int j=0; j<wds->nodes.size(); j++)
      {
         int idx = wds->nodes[j]->getProperty("segment");
         double d = wds->nodes[j]->getProperty("demand");
         di[idx] += d;
         sumD += d;
      }
      for(int j=0; j<di.size(); j++)
      {
         di[j] /= sumD;
      }

      vector<int> pi;
      for(int j=0; j<wds->presIndex.size(); j++)
      {
      	pi.push_back(wds->edges[wds->presIndex[j]]->getEdgeIntProperty("segment"));
      }
      for(int j=0; j<wds->poolIndex.size(); j++)
      {
      	pi.push_back(wds->edges[wds->poolIndex[j]]->getEdgeIntProperty("segment"));
      }

      makeDirectory("Network Data");
      makeDirectory("Network Data/" + caseName);

      ofstream wFile;
      wFile.open("Network Data/" + caseName + "/segment_graph.txt");
      for(int j=0; j<sev.size() ; j=j+2)
      {
         wFile << sev[j] << ", " << sev[j+1] << '\n';
      } 
      wFile.close();

      wFile.open("Network Data/" + caseName + "/relative_demand_loss.txt");
      for(int j=0; j<bi.size() ; j++)
      {
         wFile << bi[j] << '\n';
      } 
      wFile.close();

      wFile.open("Network Data/" + caseName + "/relative_demand.txt");
      for(int j=0; j<di.size() ; j++)
      {
         wFile << di[j] << '\n';
      } 
      wFile.close();

      wFile.open("Network Data/" + caseName + "/pressure_points.txt");
      for(int j=0; j<pi.size() ; j++)
      {
      	wFile << pi[j] << '\n';
      } 
      wFile.close();
   }

   cout << endl << endl;
   return 0;
}