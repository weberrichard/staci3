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
   string caseFolder = "../../Networks/ky/";

   int numberN1 = 1; // number of random N-1 iso placing
   double nominalISODiameter = 0.2;
   double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;

   vector<string> everyCase;
   //everyCase.push_back("Anytown");
   everyCase.push_back("ky1");
   everyCase.push_back("ky2");
   everyCase.push_back("ky3");
   everyCase.push_back("ky4");
   everyCase.push_back("ky5");
   everyCase.push_back("ky6");
   everyCase.push_back("ky7");
   everyCase.push_back("ky8");
   everyCase.push_back("ky9");
   everyCase.push_back("ky10");
   everyCase.push_back("ky11");
   everyCase.push_back("ky12");
   everyCase.push_back("ky13");
   everyCase.push_back("ky14");

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
      for(int j=0; j<wds->nodes.size(); j++)
      {
         wds->nodes[j]->demand = wds->nodes[j]->demand*1.0;
      }

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
               if(tc == 1 || tc ==0)
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
               if(tc == 1 || tc ==0)
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

     // saving original vulnerability
      vector<double> vulner;
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

         // writing the relative demand losses to file
      writeVectorDouble("Network Data/" + caseName + "/relativeDemandLoss.txt",wds->relativeDemandLoss);
      wds->calculateVulnerabilityApprox();
      writeVectorDouble("Network Data/" + caseName + "/sumSegments.txt", wds->relativeSegmentLoss);
      writeVectorDouble("Network Data/" + caseName + "/sumPipes.txt", wds->relativeLengthLoss);
      writeVectorDouble("Network Data/" + caseName + "/sumCons.txt", wds->relativeDemandLossApprox);

      // ------------------------
      //         N-1 RULE
      // ------------------------

      /*vulner.clear();
      // placing iso valves numberN1 times
      for(unsigned int I=0; I<numberN1; I++)
      {
         cout << "\n\n N-1 rule vulnerability calculation (" << I << "/" << numberN1 << ")";
         cout << "\n------------------------------------";
         // deleting every existing 
         vector<int> ISOValvesToDelete;
         for(unsigned int j=0; j<wds->valveISOIndex.size(); j++)
         {
            ISOValvesToDelete.push_back(wds->valveISOIndex[j]);
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

         vulner.clear();
         for(unsigned int j=0; j<wds->getNumberSegment(); j++)
         {
            vulner.push_back(wds->localGamma[j]);
         }

         // writing to file for further analysis
         wFile.open("Network Data/" + caseName + "/vulner_Nm1.txt");
         for(int j=0; j<vulner.size(); j++)
         {
            wFile << vulner[j] << '\n';
         }
         wFile.close();
      }*/
   }

   cout << endl << endl;
   return 0;
}

