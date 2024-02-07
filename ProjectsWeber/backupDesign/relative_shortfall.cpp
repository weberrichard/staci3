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

   string caseName, runType;
   if(argc != 3)
   {
      cout << "Wrong number of inputs, correct: 2 (caseName, runType[pn,of,of2])" << endl;
      exit(-1);
   }
   else
   {
      caseName = argv[1];
      runType = argv[2];
   }
   
   Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");
   if(runType == "pn") // pipe number
   {
      int pn = wds->pipeIndex.size();
      ofstream wFile;
      wFile.open("pn.txt");
      wFile << pn << endl;
      wFile.close();
   }
   else if(runType == "of") // classic optimisation
   {
      double cost=0., cost_orig=0.;
      double Arp=13.,Brp=29.,Crp=1200.;
      vector<double> diameter = readVectorDouble("diameter.txt");
      for(int i=0; i<wds->pipeIndex.size(); i++)
      {  
         double length=wds->edges[wds->pipeIndex[i]]->getDoubleProperty("length");

         double diam_orig =  wds->edges[wds->pipeIndex[i]]->getDoubleProperty("diameter");
         cost_orig += (Arp + Brp*diam_orig + Crp*diam_orig*diam_orig)*length;

         wds->edges[wds->pipeIndex[i]]->setDoubleProperty("diameter",diameter[i]*1.e-3);
         cost += (Arp + Brp*diameter[i]*1.e-3 + Crp*diameter[i]*diameter[i]*1.e-6)*length;
      }
      double cost_rel = cost/cost_orig;

      wds->isPressureDemand = true;
      wds->solveSystem();
      double rs=0., sumd=0.;
      for(int i=0; i<wds->nodes.size(); i++)
      {
         double d = wds->nodes[i]->getProperty("demand");
         if(d!=0)
         {
            double c = wds->nodes[i]->getProperty("consumption");
            rs += (d-c);
            sumd += d;
         }
      }
      rs /= sumd;

      double a = 1000.;
      double b = 1.;
      double of = a*rs + b*cost_rel;

      // cout << "rs : " << rs << " c: " << cost_rel << endl;

      ofstream wFile;
      wFile.open("of.txt");
      wFile << of << endl;
      wFile.close();
   }
   else if(runType == "of2") // backup optimisation
   {
      double cost=0., cost_orig=0.;
      double Arp=13.,Brp=29.,Crp=1200.;
      vector<double> diameter = readVectorDouble("diameter2.txt");
      for(int i=0; i<wds->pipeIndex.size(); i++)
      {  
         double length=wds->edges[wds->pipeIndex[i]]->getDoubleProperty("length");

         double diam_orig =  wds->edges[wds->pipeIndex[i]]->getDoubleProperty("diameter");
         cost_orig += (Arp + Brp*diam_orig + Crp*diam_orig*diam_orig)*length;

         wds->edges[wds->pipeIndex[i]]->setDoubleProperty("diameter",diameter[i]*1.e-3);
         cost += (Arp + Brp*diameter[i]*1.e-3 + Crp*diameter[i]*diameter[i]*1.e-6)*length;
      }
      double cost_rel = cost/cost_orig;

      // cout << "cost: " << cost << endl;
      // cout << "cost_orig: " << cost_orig << endl;

      wds->calculateVulnerability();
      double rs = wds->backupRelativeShortfall;

      // cout << "rs: " << rs << endl;

      double a = 1000.;
      double b = 1.;
      double of2 = a*rs + b*cost_rel;

      //cout << "of2: " << of2 << endl;

      ofstream wFile;
      wFile.open("of2.txt");
      wFile << of2 << endl;
      wFile.close();
   }
   
   return 0;
}