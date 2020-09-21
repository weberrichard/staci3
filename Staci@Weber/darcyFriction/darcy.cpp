#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../HydraulicSolver.h"
#include "../../Shutdown.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){

  // Name of containing folder of staci file
  string caseFolder = "../../Networks/Basic/";
  //string caseFolder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";

  string caseName;
  if(argc == 1){
    //caseName = "C-town.inp";
    //caseName = "Anytown.inp";
    caseName =  "linear_3.inp";
    //caseName =  "ky2.inp";
    //caseName =  "Net1.inp";
  }else if(argc == 2){
    caseName = argv[1];
  }

  cout << endl << "Case: " << caseName << endl;
  srand((unsigned int) time(0));

  HydraulicSolver *wds;
  wds = new HydraulicSolver(caseFolder + caseName);
  wds->printLevel = 3;
  wds->initialization();
  wds->solveSystem();
  cout << endl << "vf: " << endl;
  for(int i=0; i<wds->edges.size(); i++)
    cout << wds->edges[i]->name << "  " << wds->edges[i]->getDoubleProperty("velocity") << endl;
  cout << endl << "head: " << endl;
  for(int i=0; i<wds->nodes.size(); i++)
    cout << wds->nodes[i]->name << "  " << wds->nodes[i]->getProperty("head") << endl;

  cout << endl << endl;
  return 0;
}