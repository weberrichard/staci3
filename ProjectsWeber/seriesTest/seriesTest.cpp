#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../SeriesHydraulics.h"
#include "../../Shutdown.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){

  srand((unsigned int) time(0));

  // Name of containing folder of staci file
  string caseFolder = "../../Networks/";
  //string caseFolder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";

  string caseName;
  if(argc == 1){
    //caseName = "C-town.inp";
    //caseName = "Anytown.inp";
    caseName = "MUModel.inp";
    //caseName =  "grid_9.inp";
    //caseName =  "ky2.inp";
    //caseName =  "Net1.inp";
  }else if(argc == 2){
    caseName = argv[1];
  }

  cout << endl << " Case name: " << caseName << endl;

  SeriesHydraulics *sh;
  sh = new SeriesHydraulics(caseFolder + caseName);

  sh->printLevel = 3;
  clock_t ido = clock();
  sh->seriesSolve();
  cout << endl << "\n Sol time:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  vector<string> edgesID, nodesID;
  /*int npp = sh->poolIndex.size() + sh->presIndex.size();
  for(int i=sh->numberNodes-npp; i<sh->numberNodes; i++)
  {
    nodesID.push_back(sh->nodes[i]->name);
    edgesID.push_back(sh->nodes[i]->name);
  }*/

  //edgesID.push_back("WTP-TURITEA");
  //edgesID.push_back("1");
  //edgesID.push_back("3");
  //edgesID.push_back("RSV-NGAHERE-INFLOW");
  //for(int i=0; i<sh->valveIndex.size(); i++)
  // edgesID.push_back(sh->edges[sh->valveIndex[i]]->name);

  sh->saveToFile(edgesID,nodesID,"lps","m");

  //vector<string> tableEdge;
  //tableEdge.push_back("T5");
  //sh->timeTableEdge(tableEdge,"lps");

  vector<string> tableNode;
  //tableNode.push_back("node_1");
  //tableNode.push_back("node_2");
  //tableNode.push_back("node_3");
  //tableNode.push_back("12");
  //tableNode.push_back("22");
  //tableNode.push_back("32");
  tableNode.push_back("node_13");
  //tableNode.push_back("23");
  //tableNode.push_back("31");
  //tableNode.push_back("32");
  sh->timeTableNode(tableNode,"m");

  cout << endl << endl;
  return 0;
}
