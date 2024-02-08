#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../HydraulicSolver.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
   // Name of containing folder of staci file
   string caseFolder = "../../Networks/Sopron/";

   string caseName;
   if(argc != 2)
   {
      cout << "Wrong number of inputs, correct: 2 (caseName)" << endl;
      exit(-1);
   }
   else
   {
      caseName = argv[1];
   }
   
   HydraulicSolver *wds = new HydraulicSolver(caseFolder + caseName + ".inp");
   wds->solveSystem();
   
   HydraulicSolver *wds_simp = new HydraulicSolver(caseFolder + caseName + "_simp.inp");
   wds_simp->solveSystem();

   double pres_error=0.;
   for(int i=0; i<wds_simp->numberNodes; i++)
   {
      int idx = wds->nodeIDtoIndex(wds_simp->nodes[i]->name, false);

      if(idx>-1)
      {
         pres_error += abs(wds->nodes[idx]->head - wds_simp->nodes[i]->head)/wds->nodes[idx]->head;
      }
   }
   cout << wds_simp->numberNodes << endl;
   pres_error /= wds_simp->numberNodes;

   double flow_error=0.;
   for(int i=0; i<wds_simp->numberEdges; i++)
   {
      int idx = wds->edgeIDtoIndex(wds_simp->edges[i]->name, false);

      if(idx>-1 && abs(wds->edges[idx]->volumeFlowRate)>0.)
      {
         flow_error += abs(wds->edges[idx]->volumeFlowRate - wds_simp->edges[i]->volumeFlowRate)/abs(wds->edges[idx]->volumeFlowRate);
      }
   }
   flow_error /= wds_simp->numberEdges;

   cout << "[*] " << caseName << endl;
   cout << "pres_error: " << pres_error << " flow_error: " << flow_error << endl;

   return 0;
}