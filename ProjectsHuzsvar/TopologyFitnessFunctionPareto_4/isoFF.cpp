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
   ofstream wFile("results.txt");


   double nominalISODiameter = 0.2;
   double refA = nominalISODiameter*nominalISODiameter*M_PI/4.;

   if(argc < 2)
   {
      cout << " !ERROR! At least two input arguments is required" << endl;
      exit(-1);
   }

   string caseName = argv[1];
   string runType = argv[2];
   
   Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");

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
   else if(runType == "Valve") // number of iso valves at default
   {
      double n = wds->valveISOIndex.size();
      wFile << n << '\n';
   }
   else if(runType == "GetOriginalGamma")
   {
      wds->buildSegmentGraph();
      wds->calculateVulnerability();
      ofstream wFile2("results/results_"+caseName+"_originalGammaDistribution.txt");
      for (int i = 0; i < wds->getNumberSegment(); ++i)
      {
         wFile2 << wds->localGamma.at(i) << '\n';
         cout << wds->localGamma.at(i) << '\n';
      }
      wFile2.close();  
   }
   // fitness functions //
   // ----------------- //

   if(runType == "Topology" || runType == "Topology_Length" || runType == "Topology_Volume" || runType == "Topology_Demand" || runType == "Gamma" || runType == "Sigma_gamma" || runType == "gamma_distribution" || runType == "get_Topology")
   {
      // deleting every existing 
      vector<int> ISOValvesToDelete;
      for(unsigned int i=0; i<wds->valveISOIndex.size(); i++)
      {
         ISOValvesToDelete.push_back(wds->valveISOIndex[i]);
      }
      wds->deleteISOValves(ISOValvesToDelete);
      // todo string helyett index

      // loading iso valve positions
      string inFileName = argv[3];
      vector<string> fileData = readVectorString(inFileName);
      vector<int> addISOPipe;
      vector<bool> isStart;
      vector<string> addISOName;
      vector<double> refCros;
      for(int i=0; i<fileData.size(); i++)
      {
         stringstream ss(fileData[i]);
         vector<string> sv;
         cout << fileData[i] << endl;
         while(ss.good())
         {
            string substr;
            getline(ss,substr,',');
            sv.push_back(substr);
         }
         cout << sv[0] << "," << sv[1] << endl;
         addISOPipe.push_back(stoi(sv[0]));
         isStart.push_back(stoi(sv[1]));
         addISOName.push_back("ISO_" + to_string(i));
         refCros.push_back(refA);
      }

      // adding the new iso valves
      wds->addNewISOValves(addISOName, addISOPipe, isStart, 1000., refCros, 0.);
      wds->buildSegmentGraph();

      if(runType == "Topology") // lambda = a*f1 + b*f2, f1=1/ns, f2 = ns+sum_i(Li/Lmax)
      {
         double a,b;
         if(strcmp(argv[4],"a") == 0)
         {
            a = stoi(argv[5]);
         }
         else if(strcmp(argv[4],"b") == 0)
         {
            b = stoi(argv[5]);
         }

         if(strcmp(argv[6],"a") == 0)
         {
            a = stoi(argv[7]);
         }
         else if(strcmp(argv[6],"b") == 0)
         {
            b = stoi(argv[7]);
         }

         double ns = wds->numberSegment;
         double f1 = 1./ns;

         // finding Lmax
         double Lmax = 0.;
         for(int i=0; i<ns; i++)
         {  
            double Li = wds->absolutePipeLength[i];
            if(Lmax<Li)
            {
               Lmax = Li;
            }
         }

         // calculating lamda for fitness function
         double f2=0.;
         for(int i=0; i<ns; i++)
         {
            double Li = wds->absolutePipeLength[i];
            f2 += Li/Lmax;
            cout << i << " . szegmens " << f2 << endl;
         }
         f2 = ns - f2;
         double sigma = standardDeviation(wds->absolutePipeLength);

         // write to file
         double lambda = a*f1 + b*f2;
         wFile << sigma << '\n';
         cout << sigma << endl;
      }
      else if(runType == "Topology_Length") // lambda = a*f1 + b*f2, f1=1/ns, f2 = ns+sum_i(Li/Lmax)
      {
         double ns = wds->numberSegment;
         double Lmax = 0.;
         vector<double> DL;
         for(int i=0; i<ns; i++)
         {  
            double Li = wds->absolutePipeLength[i];
            if(Lmax<Li)
            {
               Lmax = Li;
            }
         }

         for(int i=0; i<ns; i++)
         {
            double Li = wds->absolutePipeLength[i];
            DL.push_back(Li/Lmax);
            cout << i << " . szegmens " << DL[i] << endl;
         }
         double sigma = standardDeviation(DL);
         wFile << sigma << '\n';
         cout << sigma << endl;
      }
      else if(runType == "Topology_Volume")
      {
         double ns = wds->numberSegment;
         vector<double> Q_i(wds->numberSegment,0.);
         vector<double> Q_i_rel(wds->numberSegment,0.);
         for(int i=0; i<wds->edges.size(); i++)
         {
            if(wds->edges[i]->typeCode == 1 || wds->edges[i]->typeCode == 0) // pipe, pipeCV
            {
              double L = wds->edges[i]->getDoubleProperty("length");
              double D = wds->edges[i]->getDoubleProperty("diameter");
              Q_i[wds->edges[i]->segment] += L*D;
            }
         }
         double Qmax = 0.;
         for(int i=0; i<ns; i++)
         {            
            if(Qmax<Q_i[i])
            {
               Qmax = Q_i[i];
            }
         }
         for(int i=0; i<ns; i++)
         {
            Q_i_rel[i] = Q_i[i]/Qmax;
            cout << i << " . szegmens " << Q_i_rel[i] << endl;
         }
         double sigma = standardDeviation(Q_i_rel);
         wFile << sigma << '\n';
         cout << sigma << endl;
      }
      else if(runType == "Topology_Demand")
      {
         double ns = wds->numberSegment;
         vector<double> D_i(wds->numberSegment,0.);
         vector<double> D_i_rel(wds->numberSegment,0.);
         vector<bool> Allocated_Demand(wds->nodes.size(),0.);
         vector< vector<int> > SegmentVector = wds->getSegmentVector();
         for(int i=0; i< SegmentVector.size(); i++)
         {
            for (int j = 0; j < SegmentVector[i].size(); ++j)
            {
               if (Allocated_Demand[SegmentVector[i][j]] != true)
               {
                  D_i[i] += wds->nodes.at(SegmentVector[i][j])->getProperty("demand");
               }
            }
         }
         double Dmax = 0.;
         for(int i=0; i<ns; i++)
         {            
            if(Dmax<D_i[i])
            {
               Dmax = D_i[i];
            }
         }
         for(int i=0; i<ns; i++)
         {
            D_i_rel[i] = D_i[i]/Dmax;
            cout << i << " . szegmens " << D_i_rel[i] << endl;
         }
         double sigma = standardDeviation(D_i_rel);
         wFile << sigma << '\n';
         cout << sigma << endl;
      } 
      else if(runType == "Gamma")
      {
         wds->calculateVulnerability();
         wFile << wds->globalGamma << '\n';
      }
      else if(runType == "Sigma_gamma")
      {  
         wds->calculateVulnerability();
         double sigma = standardDeviation(wds->localGamma);
         wFile << sigma << '\n';
      }
      else if(runType == "gamma_distribution")
      {  
         wds->calculateVulnerability();
         double counter = 0;
         ofstream wFile3("results/results_"+caseName+"_modifiedGammaDistribution.txt");
         for (int i = 0; i < wds->getNumberSegment(); ++i)
         {
            wFile3 << wds->localGamma[i] << '\n';
            cout << wds->localGamma[i] << '\n';
         }
         wFile3.close();
      }
      else if(runType == "get_Topology")
      {  
         wds->calculateVulnerability();
         double ns = wds->numberSegment;
         double f1 = 1./ns;
         ofstream wFile3("results/results_"+caseName+"_"+argv[4]+"_TopologyFitnessFunction.txt");
         // finding Lmax
         double Lmax = 0.;
         for(int i=0; i<ns; i++)
         {  
            double Li = wds->absolutePipeLength[i];
            if(Lmax<Li)
            {
               Lmax = Li;
            }
         }

         // calculating lamda for fitness function
         double f2=0.;
         for(int i=0; i<ns; i++)
         {
            double Li = wds->absolutePipeLength[i];
            f2 -= Li/Lmax;
         }
         f2 += ns;

         // write to file
         double lambda = f1 + f2;
         wFile3 << f1 << "," << f2 << "," << wds->globalGamma << "," << argv[4] << '\n';
         wFile3.close();
      }

   }
   wFile.close();   
   cout << endl << endl;
   return 0;
}

