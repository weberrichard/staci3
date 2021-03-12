#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>

#include "../../DMA.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){

  // Name of containing folder of staci file
  string case_folder = "../../Networks/Sopron/";
  //string case_folder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";

  string case_name;
  if(argc == 1){
    case_name = "balf.inp";
    //case_name = "grid_9.inp";
    //case_name = "villasor.inp";
    //case_name = "Anytown.inp";
    //case_name =  "linear_9.inp";
    //case_name =  "grid_9.inp";
    //case_name =  "ky2.inp";
    //case_name =  "Net1.inp";
  }else if(argc == 2){
    case_name = argv[1];
  }

  cout << endl << "Case: " << case_name << endl;

  DMA *wds;
  wds = new DMA(case_folder + case_name);

  // adding flow meters
  vector<string> flowMeterName; flowMeterName.push_back("FM1");
  vector<string> pipeName; pipeName.push_back("PIPE_1406709");
  vector<bool> isStart; isStart.push_back(true);
  double dens = wds->density;
  vector<double> refCros; refCros.push_back(0.2*0.2*M_PI/4.);

  wds->addNewFlowMeter(flowMeterName, pipeName, isStart, dens, refCros, 0.);

  wds->determineDMAZones();
  for(int i=0; i<wds->numberEdges; i++)
  {
    cout << wds->edges[i]->name << "  " << wds->edges[i]->DMAZone << endl;
  }
  cout << "DMA: " << wds->numberDMAZones << endl;

  cout << endl << endl;
  return 0;
}