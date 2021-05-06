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

int main(int argc, char* argv[])
{
   vector<string> caseName;
   caseName.push_back("tomalom"); 
   caseName.push_back("tomalom_sigma_opt"); 
   caseName.push_back("szakov"); 
   caseName.push_back("szakov_sigma_opt");

   string caseFolder = "../../Networks/Sopron/";
   ofstream wFile;

   for(int I=0; I<caseName.size(); I++)
   {
      cout << " case: " << caseName[I] << endl;
      vector<double> isoLength;
      HydraulicSolver *wds = new HydraulicSolver(caseFolder + caseName[I] + ".inp");

      for(int i=0; i<wds->valveISOIndex.size(); i++)
      {
         int index = wds->valveISOIndex[i];
         int ns = wds->edges[index]->startNodeIndex;
         int ne = wds->edges[index]->endNodeIndex;

         double xs = wds->nodes[ns]->xPosition;
         double ys = wds->nodes[ns]->yPosition;
         double xe = wds->nodes[ne]->xPosition;
         double ye = wds->nodes[ne]->yPosition;

         double x = (xs+xe)*.5;
         double y = (ys+ye)*.5;

         double lmin = 1e100;
         for(int j=0; j<wds->presIndex.size(); j++)
         {
            int in = wds->presIndex[j];
            int n = wds->edges[in]->startNodeIndex;
            double xp = wds->nodes[n]->xPosition;
            double yp = wds->nodes[n]->yPosition;

            double l = pow((xp-x)*(xp-x) + (yp-y)*(yp-y),.5);
            if(l<lmin)
            {
               lmin=l;
            }
         }
         for(int j=0; j<wds->poolIndex.size(); j++)
         {
            int in = wds->poolIndex[j];
            int n = wds->edges[in]->startNodeIndex;
            double xp = wds->nodes[n]->xPosition;
            double yp = wds->nodes[n]->yPosition;

            double l = pow((xp-x)*(xp-x) + (yp-y)*(yp-y),.5);
            if(l<lmin)
            {
               lmin=l;
            }
         }
         isoLength.push_back(lmin);
      }

      mkdir("Network Data",0777);
      mkdir(("Network Data/" + caseName[I]).c_str(),0777);

      wFile.open("Network Data/" + caseName[I] + "/isoLength.txt");
      for(int i=0; i<isoLength.size(); i++)
      {
         wFile << isoLength[i] << '\n';
      }
      wFile.close();
   }

   cout << endl << endl;
   return 0;
}

