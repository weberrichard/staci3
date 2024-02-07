#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../Staci.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
   // Name of containing folder of staci file
   string caseFolder = "../../Networks/Sopron/";

   string caseName;
   double demandCut;
   if(argc != 3)
   {
      cout << "Wrong number of inputs, correct: 2 (caseName, demandCut)" << endl;
      exit(-1);
   }
   else
   {
      caseName = argv[1];
      demandCut = stod(argv[2],0);
   }
   
   Staci *wds = new Staci(caseFolder + caseName + ".inp");
   wds->simplifySystem(demandCut);
   wds->saveSystem(caseFolder + caseName + "_simp.inp");
   
   return 0;
}