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
  string caseFolder = "../../Networks/Sopron/";
  //string caseFolder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";

  string caseName;
  if(argc == 1){
    //caseName = "C-town.inp";
    //caseName = "Anytown.inp";
    //caseName = "m_tv_k12_200324.spr";
    caseName =  "tomalom.inp";
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
  wds->saveResult("head","Node");
  wds->saveResult("volumeFlowRateAbsLPS","Pipe");

  /*Shutdown *sd;
  sd = new Shutdown(caseFolder + caseName);
  sd->buildSegmentGraph();
  sd->saveResult("segment","All");*/

  cout << endl << endl;
  return 0;
}