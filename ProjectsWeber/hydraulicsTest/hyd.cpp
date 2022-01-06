#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>

#include "../../SeriesHydraulics.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){

  // Name of containing folder of staci file
  //string case_folder = "../../Networks/ky/";
  string case_folder = "../../Networks/Sopron/";

  string case_name;
  if(argc == 1){
    //case_name = "Net1orig.inp";
    case_name = "ky2.inp";
    //case_name = "C-town.inp";
  }else if(argc == 2){
    case_name = argv[1];
  }

  cout << endl << "Case: " << case_name << endl;

  // double demConv = 1000.;
  double demConv = 15850.323141;

  double headConv = 1.4219702063247;

  SeriesHydraulics *wds;
  wds = new SeriesHydraulics(case_folder + case_name);
  wds->isPressureDemand = false;
  //wds->solveSystem();
  wds->seriesSolve();

  cout << wds->nodes[0]->name << ": " << wds->nodes[0]->demand*demConv << ", " << wds->nodes[0]->head*headConv << endl;
  cout << wds->nodes[1]->name << ": " << wds->nodes[1]->demand*demConv << ", " << wds->nodes[1]->head*headConv << endl;

  cout << endl;
  cout << wds->edges[0]->name << ": " << wds->edges[0]->volumeFlowRate*demConv << endl;
  cout << wds->edges[1]->name << ": " << wds->edges[1]->volumeFlowRate*demConv << ", " << (wds->edges[1]->getDoubleProperty("bottomLevel")+wds->edges[1]->getDoubleProperty("waterLevel"))/0.3048 << ", " << wds->edges[1]->getDoubleProperty("waterLevel")*headConv << endl;
  cout << wds->edges[2]->name << ": " << wds->edges[2]->volumeFlowRate*demConv << ", " << (wds->edges[2]->getDoubleProperty("bottomLevel")+wds->edges[2]->getDoubleProperty("waterLevel"))/0.3048 << ", " << wds->edges[2]->getDoubleProperty("waterLevel")*headConv << endl;
  cout << wds->edges[3]->name << ": " << wds->edges[3]->volumeFlowRate*demConv << ", " << (wds->edges[3]->getDoubleProperty("bottomLevel")+wds->edges[3]->getDoubleProperty("waterLevel"))/0.3048 << ", " << wds->edges[3]->getDoubleProperty("waterLevel")*headConv << endl;

  cout << endl;
  cout << wds->edges[4]->name << ": " << wds->edges[4]->volumeFlowRate*demConv << ", " << wds->edges[4]->getDoubleProperty("velocity")/0.3048 << endl;


  cout << endl;
  cout << wds->edges[wds->pumpIndex[0]]->name << ": " << wds->edges[wds->pumpIndex[0]]->volumeFlowRate << ", " << wds->edges[wds->pumpIndex[0]]->status << endl;

  cout << endl;
  cout << "friction model: " << wds->frictionModel << endl;
  cout << wds->edges[10]->name << ": " << wds->edges[10]->getIntProperty("frictionModel") << endl;

  double sumDem=0.;
  for(int i=0; i<wds->nodes.size(); i++)
  	sumDem += wds->nodes[i]->demand;

  cout << "sum dem: " << sumDem*demConv << endl;

  cout << endl << endl;
  return 0;
}