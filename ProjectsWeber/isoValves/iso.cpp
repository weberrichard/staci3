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

   int numberN1 = 1; // number of random N-1 iso placing
   double nominalISODiameter = 0.2;
   double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;

   vector<string> everyCase;
   everyCase.push_back("villasor");
   //everyCase.push_back("ferto");
   /*everyCase.push_back("sanchegy");
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
   everyCase.push_back("ujhermes");*/

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

      // for containing number of ISO valves/segments: original, N, N-1
      vector<int> numberOfISO;
      vector<int> numberOfSegment;

      string caseName = everyCase[i];
      makeDirectory(("Network Data/" + caseName).c_str());

      cout << "\n\n Original vulnerability calculation";
      cout << "\n------------------------------------";
      Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");
      wds->calculateVulnerability();
      numberOfISO.push_back(wds->valveISOIndex.size());
      numberOfSegment.push_back(wds->getNumberSegment());

      // saving original vulnerability
      vector<double> vulner;
      for(unsigned int j=0; j<wds->getNumberSegment(); j++)
      {
         vulner.push_back(wds->localGamma[j]);
      }

      // writing to file for further analysis
      wFile.open("Network Data/" + caseName + "/vulner_orig.txt");
      for(int j=0; j<vulner.size(); j++)
      {
         wFile << vulner[j] << '\n';
      }
      wFile.close();

      // rank of segments
      wFile.open("Network Data/" + caseName + "/rank_orig.txt");
      for(int j=0; j<wds->segmentRank; j++)
      {
         wFile << wds->segmentRank[j] << '\n';
      }
      wFile.close();

      // deleting every existing 
      vector<string> ISOValvesToDelete;
      for(unsigned int j=0; j<wds->valveISOIndex.size(); j++)
      {
         ISOValvesToDelete.push_back(wds->edges[wds->valveISOIndex[j]]->name);
      }
      wds->deleteISOValves(ISOValvesToDelete);

      // determine the rank of nodes
      vector<int> rank(wds->numberNodes,0);
      for(unsigned int j=0; j<wds->numberNodes; j++)
      {
         rank[j] += wds->nodes[j]->edgeIn.size();
         rank[j] += wds->nodes[j]->edgeOut.size();
      }

      // ------------------------
      //          N RULE
      // ------------------------
      // adding iso valves at intersections according to N rule
      cout << "\n\n N rule vulnerability calculation";
      cout << "\n------------------------------------";
      int counter=0;
      vector<string> addISOName;
      vector<string> addISOPipe;
      vector<bool> isStart;
      vector<double> refCros;
      for(unsigned int j=0; j<rank.size(); j++)
      {
         if(rank[j]>=3) // intersections
         {
            for(unsigned int k=0; k<wds->nodes[j]->edgeIn.size(); k++)
            {
               int tc = wds->edges[wds->nodes[j]->edgeIn[k]]->typeCode;
               if(tc != -1 && tc != -2)
               {
                  addISOName.push_back("ISO_" + to_string(counter));
                  string pn = wds->edges[wds->nodes[j]->edgeIn[k]]->name;
                  addISOPipe.push_back(pn);
                  isStart.push_back(false);
                  refCros.push_back(refA);
                  counter++;
               }
            }
            for(unsigned int k=0; k<wds->nodes[j]->edgeOut.size(); k++)
            {
               int tc = wds->edges[wds->nodes[j]->edgeOut[k]]->typeCode;
               if(tc != -1 && tc != -2)
               {
                  addISOName.push_back("ISO_" + to_string(counter));
                  string pn = wds->edges[wds->nodes[j]->edgeOut[k]]->name;
                  addISOPipe.push_back(pn);
                  isStart.push_back(true);
                  refCros.push_back(refA);
                  counter++;
               }
            }
         }
      }

      wds->addNewISOValves(addISOName, addISOPipe, isStart, 1000., refCros, 0.);
      wds->buildSegmentGraph();
      wds->calculateVulnerability();
      numberOfISO.push_back(wds->valveISOIndex.size());
      numberOfSegment.push_back(wds->getNumberSegment());

      // saving original vulnerability
      vulner.clear();
      for(unsigned int j=0; j<wds->getNumberSegment(); j++)
      {
         vulner.push_back(wds->localGamma[j]);
      }

      // writing to file for further analysis
      wFile.open("Network Data/" + caseName + "/vulner_N.txt");
      for(int j=0; j<vulner.size(); j++)
      {
         wFile << vulner[j] << '\n';
      }
      wFile.close();

      // rank of segments
      wFile.open("Network Data/" + caseName + "/rank_N.txt");
      for(int j=0; j<wds->segmentRank; j++)
      {
         wFile << wds->segmentRank[j] << '\n';
      }
      wFile.close();

      // ------------------------
      //         N-1 RULE
      // ------------------------

      vulner.clear();
      // placing iso valves numberN1 times
      for(unsigned int I=0; I<numberN1; I++)
      {
         cout << "\n\n N-1 rule vulnerability calculation (" << I << "/" << numberN1 << ")";
         cout << "\n------------------------------------";
         // deleting every existing 
         ISOValvesToDelete.clear();
         for(unsigned int j=0; j<wds->valveISOIndex.size(); j++)
         {
            ISOValvesToDelete.push_back(wds->edges[wds->valveISOIndex[j]]->name);
         }
         wds->deleteISOValves(ISOValvesToDelete);

         // clearing vars
         addISOName.clear();
         addISOPipe.clear();
         isStart.clear();
         refCros.clear();
         counter = 0;

         for(unsigned int j=0; j<rank.size(); j++)
         {
            if(rank[j]>=3) // intersections
            {
               // not placing one ISO out of N
               int r = rand() % rank[j];
               vector<bool> placeISO(rank[j],true);
               placeISO[r] = false;
               for(unsigned int k=0; k<wds->nodes[j]->edgeIn.size(); k++)
               {
                  if(placeISO[k])
                  {
                     int tc = wds->edges[wds->nodes[j]->edgeIn[k]]->typeCode;
                     if(tc != -1 && tc != -2)
                     {
                        addISOName.push_back("ISO_" + to_string(counter));
                        string pn = wds->edges[wds->nodes[j]->edgeIn[k]]->name;
                        addISOPipe.push_back(pn);
                        isStart.push_back(false);
                        refCros.push_back(refA);
                        counter++;
                     }
                  }
               }
               for(unsigned int k=0; k<wds->nodes[j]->edgeOut.size(); k++)
               {
                  if(placeISO[k+wds->nodes[j]->edgeIn.size()])
                  {
                     int tc = wds->edges[wds->nodes[j]->edgeOut[k]]->typeCode;
                     if(tc != -1 && tc != -2)
                     {
                        addISOName.push_back("ISO_" + to_string(counter));
                        string pn = wds->edges[wds->nodes[j]->edgeOut[k]]->name;
                        addISOPipe.push_back(pn);
                        isStart.push_back(true);
                        refCros.push_back(refA);
                        counter++;
                     }
                  }
               }
            }
         }
         wds->addNewISOValves(addISOName, addISOPipe, isStart, 1000., refCros, 0.);
         wds->buildSegmentGraph();
         wds->calculateVulnerability();
         numberOfISO.push_back(wds->valveISOIndex.size());
         numberOfSegment.push_back(wds->getNumberSegment());

         // saving vulnerability
         for(unsigned int j=0; j<wds->getNumberSegment(); j++)
         {
            vulner.push_back(wds->localGamma[j]);
         }
      }
      // writing to file for further analysis
      wFile.open("Network Data/" + caseName + "/vulner_Nm1.txt");
      for(int j=0; j<vulner.size(); j++)
      {
         wFile << vulner[j] << '\n';
      }
      wFile.close();

      wFile.open("Network Data/" + caseName + "/number_of_valves.txt");
      for(unsigned int j=0; j<numberOfISO.size(); j++)
      {
         wFile << numberOfISO[j] << ',' << numberOfSegment[j] << '\n';
      }
      wFile.close();

      // rank of segments
      wFile.open("Network Data/" + caseName + "/rank_Nm1.txt");
      for(int j=0; j<wds->segmentRank; j++)
      {
         wFile << wds->segmentRank[j] << '\n';
      }
      wFile.close();
   }

   cout << endl << endl;
   return 0;
}

