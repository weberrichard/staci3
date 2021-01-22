#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>

#include "../../Sensitivity.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){

  // Name of containing folder of staci file
  string case_folder = "../../Networks/Sopron/";
  //string case_folder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";

  string case_name;
  if(argc == 1){
    //case_name = "C-town.inp";
    case_name = "villasor.inp";
    //case_name = "Anytown.inp";
    //case_name =  "linear_9.inp";
    //case_name =  "grid_9.inp";
    //case_name =  "ky2.inp";
    //case_name =  "Net1.inp";
  }else if(argc == 2){
    case_name = argv[1];
  }

  cout << endl << "Case: " << case_name << endl;
  srand((unsigned int) time(0));

  Sensitivity *wds;
  clock_t ido = clock();
  wds = new Sensitivity(case_folder + case_name);
  cout << "Nodes: " << wds->nodes.size() << endl << "Edges: " << wds->edges.size() << endl;
  cout << endl << "\nKonst:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  ido = clock();
  wds->initialization();
  cout << endl << "\nInit:    " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  ido = clock();
  wds->solveSystem();
  cout << endl << "\nSolver:  " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;
  /*printf("   name   | vf\n");
  for(int i=0; i<wds->edges.size(); i++)
    printf(" %8s | %8.3f\n", wds->edges[i]->name.c_str(), wds->edges[i]->volumeFlowRate);
  printf("   name   | head\n");
  for(int i=0; i<wds->nodes.size(); i++)
    printf(" %8s | %8.3f\n", wds->nodes[i]->name.c_str(), wds->nodes[i]->head);*/

  wds->initialization();
  ido = clock();
  wds->calculateSensitivity("demand");
  cout << endl << "\nSensitivity, demand:  " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  // debug
  //cout << endl << " Sens, mf-dem: " << endl << wds->massFlowRateSensitivity << endl;
  //cout << endl << " Sens, pres-dem: " << endl << wds->pressureSensitivity << endl;
  // debug

  wds->initialization();
  ido = clock();
  wds->calculateSensitivity("roughness");
  cout << endl << "\nSensitivity, fric:  " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  // debug
  //cout << endl << " Sens, mf-dem: " << endl << wds->massFlowRateSensitivity << endl;
  //cout << endl << " Sens, pres-dem: " << endl << wds->pressureSensitivity << endl;
  // debug

  cout << endl << endl;
  return 0;
}