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
      for(int j=0; j<wds->segmentRank.size(); j++)
      {
         wFile << wds->segmentRank[j] << '\n';
      }
      wFile.close();

      // edge vector of segment graph
      wFile.open("Network Data/" + caseName + "/segment_edge_orig.txt");
      for(int j=0; j<wds->segmentEdgeVector.size(); j+=2)
      {
         wFile << wds->segmentEdgeVector[j] << ',' << wds->segmentEdgeVector[j+1] << '\n';
      }
      wFile.close();

      // edge vector of segment graph
      wFile.open("Network Data/" + caseName + "/input_segment_orig.txt");
      for(int j=0; j<wds->presIndex.size(); j++)
      {
         int idx = wds->presIndex[j];
         wFile << wds->edges[idx]->segment << '\n';
      }
      for(int j=0; j<wds->poolIndex.size(); j++)
      {
         int idx = wds->poolIndex[j];
         wFile << wds->edges[idx]->segment << '\n';
      }
      wFile.close();

      // demand loss
      wFile.open("Network Data/" + caseName + "/demand_loss_orig.txt");
      for(int j=0; j<wds->relativeDemandLoss.size(); j++)
      {
         wFile << wds->relativeDemandLoss[j] << '\n';
      }
      wFile.close();

      // segment demand
      wFile.open("Network Data/" + caseName + "/segment_demand_orig.txt");
      vector<double> segmentDemand(wds->numberSegment,0.);
      for(int j=0; j<wds->nodes.size(); j++)
      {
         int seg = wds->nodes[j]->segment;
         segmentDemand[seg] += wds->nodes[j]->demand;
      }
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << segmentDemand[j] << '\n';
      }
      wFile.close();

      // absolute pipeline length
      wFile.open("Network Data/" + caseName + "/absolute_segment_length_orig.txt");
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << wds->absolutePipeLength[j] << '\n';
      }
      wFile.close();

      // network vulnerability
      wFile.open("Network Data/" + caseName + "/network_vulner_orig.txt");
      wFile << wds->globalGamma << '\n';
      wFile.close();

      // topological fitness function
      double ns = wds->numberSegment;
      double f1 = 1./ns;
      // finding Lmax
      double Lmax = 0.;
      for(int j=0; j<ns; j++)
      {  
         double Li = wds->absolutePipeLength[j];
         if(Lmax<Li)
         {
            Lmax = Li;
         }
      }

      // calculating lamda for fitness function
      double f2=0.;
      for(int j=0; j<ns; j++)
      {
         double Li = wds->absolutePipeLength[j];
         f2 -= Li/Lmax;
      }
      f2 += ns;

      wFile.open("Network Data/" + caseName + "/Lambda_orig.txt");
      wFile << f1 << '\n';
      wFile << f2 << '\n';
      wFile.close();

      wds->isPressureDemand=false;
      wds->solveSystem();
      wds->isPressureDemand=true;
      // saving relative demand loss
      for(unsigned int j=0; j<wds->nodes.size(); j++)
      {  
         int seg = wds->nodes[j]->segment;
         double dem_loss = wds->relativeDemandLoss[seg];
         double vul = wds->localGamma[seg];
         wds->nodes[j]->setProperty("userOutput",vul);
      }
      wds->saveResult("head","Node");

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
      /*cout << "\n\n N rule vulnerability calculation";
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
      for(int j=0; j<wds->segmentRank.size(); j++)
      {
         wFile << wds->segmentRank[j] << '\n';
      }
      wFile.close();

      // edge vector of segment graph
      wFile.open("Network Data/" + caseName + "/segment_edge_N.txt");
      for(int j=0; j<wds->segmentEdgeVector.size(); j+=2)
      {
         wFile << wds->segmentEdgeVector[j] << ',' << wds->segmentEdgeVector[j+1] << '\n';
      }
      wFile.close();

      // edge vector of segment graph
      wFile.open("Network Data/" + caseName + "/input_segment_N.txt");
      for(int j=0; j<wds->presIndex.size(); j++)
      {
         int idx = wds->presIndex[j];
         wFile << wds->edges[idx]->segment << '\n';
      }
      for(int j=0; j<wds->poolIndex.size(); j++)
      {
         int idx = wds->poolIndex[j];
         wFile << wds->edges[idx]->segment << '\n';
      }
      wFile.close();

      // demand loss
      wFile.open("Network Data/" + caseName + "/demand_loss_N.txt");
      for(int j=0; j<wds->relativeDemandLoss.size(); j++)
      {
         wFile << wds->relativeDemandLoss[j] << '\n';
      }
      wFile.close();

      // network vulnerability
      wFile.open("Network Data/" + caseName + "/network_vulner_N.txt");
      wFile << wds->globalGamma << '\n';
      wFile.close();

      // segment demand
      wFile.open("Network Data/" + caseName + "/segment_demand_N.txt");
      segmentDemand.clear();
      segmentDemand.resize(wds->numberSegment);
      for(int j=0; j<wds->nodes.size(); j++)
      {
         int seg = wds->nodes[j]->segment;
         segmentDemand[seg] += wds->nodes[j]->demand;
      }
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << segmentDemand[j] << '\n';
      }
      wFile.close();

      // absolute pipeline length
      wFile.open("Network Data/" + caseName + "/absolute_segment_length_N.txt");
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << wds->absolutePipeLength[j] << '\n';
      }
      wFile.close();

      // topological fitness function
      ns = wds->numberSegment;
      f1 = 1./ns;
      // finding Lmax
      Lmax = 0.;
      for(int j=0; j<ns; j++)
      {  
         double Li = wds->absolutePipeLength[j];
         if(Lmax<Li)
         {
            Lmax = Li;
         }
      }

      // calculating lamda for fitness function
      f2=0.;
      for(int j=0; j<ns; j++)
      {
         double Li = wds->absolutePipeLength[j];
         f2 -= Li/Lmax;
      }
      f2 += ns;

      wFile.open("Network Data/" + caseName + "/Lambda_N.txt");
      wFile << f1 << '\n';
      wFile << f2 << '\n';
      wFile.close();

      // saving relative demand loss
      /*for(unsigned int j=0; j<wds->nodes.size(); j++)
      {  
         int seg = wds->nodes[j]->segment;
         double dem_loss = wds->relativeDemandLoss[seg];
         double vul = wds->localGamma[seg];
         wds->nodes[j]->setProperty("userOutput",vul);
      }
      wds->saveResult("userOutput","Node");*/

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
      for(int j=0; j<wds->segmentRank.size(); j++)
      {
         wFile << wds->segmentRank[j] << '\n';
      }
      wFile.close();

      // edge vector of segment graph
      wFile.open("Network Data/" + caseName + "/segment_edge_Nm1.txt");
      for(int j=0; j<wds->segmentEdgeVector.size(); j+=2)
      {
         wFile << wds->segmentEdgeVector[j] << ',' << wds->segmentEdgeVector[j+1] << '\n';
      }
      wFile.close();

      // edge vector of segment graph
      wFile.open("Network Data/" + caseName + "/input_segment_Nm1.txt");
      for(int j=0; j<wds->presIndex.size(); j++)
      {
         int idx = wds->presIndex[j];
         wFile << wds->edges[idx]->segment << '\n';
      }
      for(int j=0; j<wds->poolIndex.size(); j++)
      {
         int idx = wds->poolIndex[j];
         wFile << wds->edges[idx]->segment << '\n';
      }
      wFile.close();

      // demand loss
      wFile.open("Network Data/" + caseName + "/demand_loss_Nm1.txt");
      for(int j=0; j<wds->relativeDemandLoss.size(); j++)
      {
         wFile << wds->relativeDemandLoss[j] << '\n';
      }
      wFile.close();

      // network vulnerability
      wFile.open("Network Data/" + caseName + "/network_vulner_Nm1.txt");
      wFile << wds->globalGamma << '\n';
      wFile.close();

      // segment demand
      wFile.open("Network Data/" + caseName + "/segment_demand_Nm1.txt");
      segmentDemand.clear();
      segmentDemand.resize(wds->numberSegment);
      for(int j=0; j<wds->nodes.size(); j++)
      {
         int seg = wds->nodes[j]->segment;
         segmentDemand[seg] += wds->nodes[j]->demand;
      }
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << segmentDemand[j] << '\n';
      }
      wFile.close();

      // absolute pipeline length
      wFile.open("Network Data/" + caseName + "/absolute_segment_length_Nm1.txt");
      for(int j=0; j<wds->numberSegment; j++)
      {
         wFile << wds->absolutePipeLength[j] << '\n';
      }
      wFile.close();

      // topological fitness function
      ns = wds->numberSegment;
      f1 = 1./ns;
      // finding Lmax
      Lmax = 0.;
      for(int j=0; j<ns; j++)
      {  
         double Li = wds->absolutePipeLength[j];
         if(Lmax<Li)
         {
            Lmax = Li;
         }
      }

      // calculating lamda for fitness function
      f2=0.;
      for(int j=0; j<ns; j++)
      {
         double Li = wds->absolutePipeLength[j];
         f2 -= Li/Lmax;
      }
      f2 += ns;

      wFile.open("Network Data/" + caseName + "/Lambda_Nm1.txt");
      wFile << f1 << '\n';
      wFile << f2 << '\n';
      wFile.close();

      // saving relative demand loss
      /*for(unsigned int j=0; j<wds->nodes.size(); j++)
      {  
         int seg = wds->nodes[j]->segment;
         double dem_loss = wds->relativeDemandLoss[seg];
         double vul = wds->localGamma[seg];
         wds->nodes[j]->setProperty("userOutput",vul);
      }
      wds->saveResult("userOutput","Node");*/
   }

   cout << endl << endl;
   return 0;
}

