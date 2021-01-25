#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../Sensitivity.h"

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
    //caseName =  "tomalom.inp";
    caseName =  "varisi.inp";
    //caseName =  "ky2.inp";
    //caseName =  "Net1.inp";
  }else if(argc == 2){
    caseName = argv[1];
  }

  cout << endl << "Case: " << caseName << endl;
  srand((unsigned int) time(0));

  Sensitivity *wds;
  wds = new Sensitivity(caseFolder + caseName);
  wds->printLevel = 3;
  wds->initialization();

  wds->calculateSensitivity("roughness");

  vector<double> Snode(wds->nodes.size(),0.);
  vector<double> Sedge(wds->edges.size(),0.);

  for(int i=0; i<wds->nodes.size(); i++)
  {
    for(int j=0; j<wds->edges.size(); j++)
    {
      Snode[i] += wds->pressureSensitivity(i,j);
      Sedge[j] += wds->pressureSensitivity(i,j);
    }
  }

  double SnodeMax=0, SedgeMax=0;
  for(int i=0; i<Snode.size(); i++)
  {
    if(abs(Snode[i]) > SnodeMax)
      SnodeMax = abs(Snode[i]);
  }
  for(int i=0; i<Sedge.size(); i++)
  {
    if(abs(Sedge[i]) > SedgeMax)
      SedgeMax = abs(Sedge[i]);
  }

  for(int i=0; i<wds->nodes.size(); i++)
  {
    wds->nodes[i]->setProperty("userOutput",abs(Snode[i])/SnodeMax);
  }
  for(int i=0; i<wds->edges.size(); i++)
  {
    wds->edges[i]->setEdgeDoubleProperty("userOutput",abs(Sedge[i])/SedgeMax);
  }

  //wds->saveResult("userOutput","Node");
  //wds->saveResult("userOutput","Pipe");
  wds->saveResult("pressure","Node");

  cout << endl << endl;
  return 0;
}