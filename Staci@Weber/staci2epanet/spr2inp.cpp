/*===================================================================*\
                              staci2epanet  
                            ----------------

  This code converts a STACI file (.spr) to an EPANET file (.inp).
  Exporting nodes, pipelines, valves, pumps, pools, pressure points.

  HOW TO USE: 

      - CASE_NAME is the name of the file without .spr extension

      - SET:  nothing :)

      - OUTPUT: CASE_NAME.inp            File for EPANET

      - RUN: ./staci2epanet CASE_NAME

      - COMPILE: make -f make_staci2epanet.mk

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>

#include "../../HydraulicSolver.h"

using namespace std;

Staci *wdsSPR;
HydraulicSolver *wdsINP;

int main(int argc, char* argv[])
{
  string case_folder = "../../Networks/";
  string case_name;
  if(argc == 1){
    case_name = "linear_3.spr";
  }else if(argc == 2){
    case_name = argv[1];
  }

  // Creating the case, loading the data
  wdsSPR = new Staci(case_folder + case_name);
  cout << endl << "[*] File loaded successfully: " << case_folder + case_name << endl;

  // Adding the elevation to pressure points
  for(int i=0; i<wdsSPR->edges.size(); i++)
  {
    if(wdsSPR->edges[i]->typeCode == -2)
    {
      int idx=-1;
      for(int j=0; j<wdsSPR->nodes.size(); j++)
        if(wdsSPR->nodes[j]->name==wdsSPR->edges[i]->startNodeName)
          idx = j;

    wdsSPR->edges[i]->setDoubleProperty("head", (wdsSPR->edges[i]->getDoubleProperty("head") + wdsSPR->nodes[idx]->getProperty("height"))); // meter to Pascal minus elevation
    }
  }

  // Warning about the isolation valves if the connecting node does not have a rank of 2
  vector<int> nodeRank(wdsSPR->nodes.size(), 0);
  vector<int> nodeRankISO(wdsSPR->nodes.size(), 0);
  for(int i=0; i<wdsSPR->edges.size(); i++)
  {
    if(wdsSPR->edges[i]->typeCode == 1) // normal pipe
    {
      int si=wdsSPR->edges[i]->startNodeIndex;
      int ei=wdsSPR->edges[i]->endNodeIndex;
      nodeRank[si]++;
      nodeRank[ei]++;
      if(wdsSPR->nodes[si]->name == "NODE_11211")
      {
        cout << endl << "edge name: " << wdsSPR->edges[i]->name << "  nodeRank[si]: " << nodeRank[si] << " si: " << si << endl;
      }
      if(wdsSPR->nodes[ei]->name == "NODE_11211")
      {
        cout << endl << "edge name: " << wdsSPR->edges[i]->name << "  nodeRank[ei]: " << nodeRank[ei] << " ei: " << ei << endl;
      }
    }
    if(wdsSPR->edges[i]->typeCode == 9) // iso valve
    {
      int si=wdsSPR->edges[i]->startNodeIndex;
      int ei=wdsSPR->edges[i]->endNodeIndex;
      nodeRankISO[si]++;
      nodeRankISO[ei]++;
    }
  }

  // Checking for double iso valve connections
  for(int i=0; i<wdsSPR->nodes.size(); i++)
  {
    if(nodeRankISO[i]>1)
    {
      cout << endl << "!WARNING! Node " << wdsSPR->nodes[i]->name << " has a rank of " << nodeRankISO[i] << " with ISO VALVES\n Do something about it... " << endl;
      //exit(-1);
    }
  }

  vector<int> isolationNodeStart;
  vector<int> isolationNodeEnd;
  vector<string> nodesToDelete;
  vector<string> edgesToDelete;
  for(int i=0; i<wdsSPR->edges.size(); i++)
  {
    if(wdsSPR->edges[i]->typeCode == 9)
    {
      int si = wdsSPR->edges[i]->startNodeIndex;
      int ei = wdsSPR->edges[i]->endNodeIndex;
      if(nodeRank[ei]>0)
      {
        cout << endl << "!WARNING! End node " << wdsSPR->edges[i]->endNodeName << " of iso valve " << wdsSPR->edges[i]->name << " has a rank of " << nodeRank[ei] << endl << " Do something... " << endl;
        //exit(-1);
      }

      // Deleting the node and the valve if the rank is one
      if(nodeRank[si] == 1)
      {
        cout << "\n Node " << wdsSPR->nodes[si]->name << " has rank of " << nodeRank[si] << ", thus it the end node and the valve will not be added to the INP" << endl << " Node name: " << wdsSPR->nodes[ei]->name << endl << " Edge name: " << wdsSPR->edges[i]->name << endl;
        edgesToDelete.push_back(wdsSPR->edges[i]->name);
        nodesToDelete.push_back(wdsSPR->nodes[ei]->name);
      }
      else if(nodeRank[si] == 2)
      {
        isolationNodeStart.push_back(si);
        isolationNodeEnd.push_back(wdsSPR->edges[i]->endNodeIndex);
      }
      else
      {
        // Warning about outstanding cases
        cout << endl << " !WARNING! WTF??? Node " << wdsSPR->nodes[si]->name << " has a rank of " << nodeRank[si] << "\n Do something about it..." << endl;
      }
    }
  }

  // Reconnecting the pipes with including the isolation valve
  for(int j=0; j<isolationNodeStart.size(); j++)
  {
    for(int i=0; i<wdsSPR->edges.size(); i++)
    {
      if(wdsSPR->edges[i]->typeCode == 1) // normal pipe
      {
        if(wdsSPR->edges[i]->startNodeIndex == isolationNodeStart[j])
        {
          wdsSPR->edges[i]->startNodeName = wdsSPR->nodes[isolationNodeEnd[j]]->name;
          break;
        }
        if(wdsSPR->edges[i]->endNodeIndex == isolationNodeStart[j])
        {
          wdsSPR->edges[i]->endNodeName = wdsSPR->nodes[isolationNodeEnd[j]]->name;
          break;
        }
      }
    }
  }

  // finding the nodes with pressurepoint or pool
  vector<bool> isRealNode(wdsSPR->nodes.size(),true);
  for(int i=0; i<wdsSPR->nodes.size(); i++)
  {
    for(int j=0; j<wdsSPR->edges.size(); j++)
    {
      int typeCode = wdsSPR->edges[j]->typeCode;
      if(typeCode == -1 || typeCode == -2)
      {
        if(wdsSPR->nodes[i]->name == wdsSPR->edges[j]->startNodeName)
        {
          isRealNode[i] = false;
          cout << endl << " Node " << wdsSPR->nodes[i]->name << " is pres/pool points, renaming to " << wdsSPR->edges[j]->name << endl;
          if(wdsSPR->nodes[i]->getProperty("demand") != 0)
          {
            cout << "!WARNING! Pressure point " << wdsSPR->nodes[i]->name << " DOES have non-zero demand " << wdsSPR->nodes[i]->demand << endl;
            exit(-1);
          }
          string node_name = wdsSPR->nodes[i]->name;
          wdsSPR->nodes[i]->name = wdsSPR->edges[j]->name;
          for(int k=0; k<wdsSPR->edges.size(); k++)
          {
            if(node_name == wdsSPR->edges[k]->startNodeName)
              wdsSPR->edges[k]->startNodeName = wdsSPR->nodes[i]->name;
            if(node_name == wdsSPR->edges[k]->endNodeName)
              wdsSPR->edges[k]->endNodeName = wdsSPR->nodes[i]->name;
          }
        }
      }
    }
  }

  // deleting unnecessary nodes
 /*for(int i=0; i<nodesToDelete.size(); i++)
  {
    for(int j=0; j<wdsSPR->nodes.size(); j++)
    {
      if(nodesToDelete[i] == wdsSPR->nodes[j]->name)
      {
        cout << endl << " Deleted node name: " << wdsSPR->nodes[j]->name << endl;
        wdsSPR->nodes.erase(wdsSPR->nodes.begin()+j);
        isRealNode.erase(isRealNode.begin()+j);
      }
    }
  }
  // deleting unnecessary edges
  for(int i=0; i<edgesToDelete.size(); i++)
  {
    for(int j=0; j<wdsSPR->edges.size(); j++)
    {
      if(edgesToDelete[i] == wdsSPR->edges[j]->name)
      {
        cout << endl << " Deleted edge name: " << wdsSPR->edges[j]->name << endl;
        wdsSPR->edges.erase(wdsSPR->edges.begin()+j);
      }
    }
  }*/

  // Renaming everything
  /*for(int i=0; i<wdsSPR->nodes.size(); i++)
  {
    string newName, oldName;
    oldName = wdsSPR->nodes[i]->name;
    if(isRealNode[i])
    {
      newName = "node_"+to_string(i);
    }
    else
    {
      int poolCoutner=0, presCoutner=0;
      for(int j=0; j<wdsSPR->edges.size(); j++)
      {
        int typeCode = wdsSPR->edges[j]->typeCode;
        if(typeCode == -1)
        {
          newName = "pool_"+to_string(poolCoutner);
          wdsSPR->edges[j]->name = newName;
          poolCoutner++;
        }
        else if(typeCode == -2)
        {
          newName = "pres_"+to_string(presCoutner);
          wdsSPR->edges[j]->name = newName;
          presCoutner++;
        }
      }
    }

    wdsSPR->nodes[i]->name = newName;
    for(int j=0; j<wdsSPR->edges.size(); j++)
    {
      if(wdsSPR->edges[j]->startNodeName == oldName)
        wdsSPR->edges[j]->startNodeName = newName;
      if(wdsSPR->edges[j]->endNodeName == oldName)
        wdsSPR->edges[j]->endNodeName = newName;
    }
  }

  int pipeCounter=0, valveCounter=0;
  for(int i=0; i<wdsSPR->edges.size(); i++)
  {
    int typeCode = wdsSPR->edges[i]->typeCode;
    string newName;
    if(typeCode == 1)
    {
      newName = "pipe_"+to_string(pipeCounter);
      wdsSPR->edges[i]->name = newName;
      pipeCounter++;
    }
    else if(typeCode == 9)
    {
      newName = "valve_"+to_string(valveCounter);
      wdsSPR->edges[i]->name = newName;
      valveCounter++;
    }
    //else if(typeCode == -1)
    //{
    //  newName = "pool_"+to_string(poolCoutner);
    //  poolCoutner++;
    //  string oldName = wdsSPR->edges[i]->name;
    //  //for(int j=0; j<wdsSPR->edges.size(); j++)
    //  //{
    //  //  if(wdsSPR->edges[j]->startNodeName == oldName)
    //  //    wdsSPR->edges[j]->startNodeName = newName;
    //  //  if(wdsSPR->edges[j]->endNodeName == oldName)
    //  //    wdsSPR->edges[j]->endNodeName = newName;
    //  //}
    //}
    //else if(typeCode == -2)
    //{
    //  newName = "pres_"+to_string(presCoutner);
    //  presCoutner++;
    //  string oldName = wdsSPR->edges[i]->name;
    //  //for(int j=0; j<wdsSPR->edges.size(); j++)
    //  //{
    //  //  if(wdsSPR->edges[j]->startNodeName == oldName)
    //  //    wdsSPR->edges[j]->startNodeName = newName;
    //  //  if(wdsSPR->edges[j]->endNodeName == oldName)
    //  //    wdsSPR->edges[j]->endNodeName = newName;
    //  //}
    //}
    else if(typeCode != -1 && typeCode != -2)
    {
      cout << endl << "!WARNING! Unhandeld edge type: " << typeCode << "  name: " << wdsSPR->edges[i]->name << endl;
      exit(-1);
    }
  }*/

  // Start of writing the file
  string file_name = case_folder + case_name.substr(0,case_name.length()-4) + ".inp";
  FILE *wfile;
  wfile = fopen(file_name.c_str(),"w");

  fprintf(wfile, "[TITLE]\n");
  fprintf(wfile, "%s\n", ("Converted from Epanet to Staci, case name: " + case_name + "\n").c_str());

  fprintf(wfile, "\n[JUNCTIONS]\n");
  fprintf(wfile, ";ID              \tElev        \tDemand      \tPattern         \n");
  for(int i=0; i<wdsSPR->nodes.size(); i++)
  {
    if(isRealNode[i])
      fprintf(wfile, " %-16s\t%-12.4f\t%-12.6f\t%-16s\t;\n", wdsSPR->nodes[i]->name.c_str(),wdsSPR->nodes[i]->getProperty("geodeticHeight"),wdsSPR->nodes[i]->getProperty("demand"),"");
  }

  fprintf(wfile, "\n[RESERVOIRS]\n");
  fprintf(wfile, ";ID              \tHead        \tPattern         \n");
  for(int i=0; i<wdsSPR->edges.size(); i++){
    if(wdsSPR->edges[i]->type=="PressurePoint"){
      fprintf(wfile," %-16s\t%-12.4f\t%-16s\n",wdsSPR->edges[i]->getEdgeStringProperty("name").c_str(),wdsSPR->edges[i]->getDoubleProperty("head"),"");
    }
  }

  // TODO: min/max level
  fprintf(wfile, "\n[TANKS]\n");
  fprintf(wfile, ";ID             \tElevation   \tInitLevel   \tMinLevel    \tMaxLevel    \tDiameter    \tMinVol      \tVolCurve\n");
  for(int i=0; i<wdsSPR->edges.size(); i++){
    if(wdsSPR->edges[i]->type=="Pool")
      fprintf(wfile, " %-16s\t%-12.4f\t%-12.4f\t%-12.4f\t%-12.4f\t%-12.4f\t%-12.4f\t%-16s\t\n", wdsSPR->edges[i]->getEdgeStringProperty("name").c_str(),wdsSPR->edges[i]->getDoubleProperty("bottomLevel"),wdsSPR->edges[i]->getDoubleProperty("waterLevel"),0.0,1000.0,pow(wdsSPR->edges[i]->getDoubleProperty("referenceCrossSection"),.5)*4./M_PI,0.,"");
  }

  fprintf(wfile, "\n[PIPES]\n"); // diameter is in mm in INP, while it is in m in SPR
  fprintf(wfile, ";ID             \tNode1           \tNode2           \tLength      \tDiameter    \tRoughness   \tMinorLoss   \tStatus\n");
  for(int i=0; i<wdsSPR->edges.size(); i++){
    if(wdsSPR->edges[i]->typeCode == 1)
      fprintf(wfile, " %-16s\t%-16s\t%-16s\t%-12.4f\t%-12.4f\t%-11.4f\t%-12.4f\t%-6s;\n",wdsSPR->edges[i]->getEdgeStringProperty("name").c_str(),wdsSPR->edges[i]->startNodeName.c_str(),wdsSPR->edges[i]->endNodeName.c_str(),wdsSPR->edges[i]->getDoubleProperty("length"),wdsSPR->edges[i]->getDoubleProperty("diameter")*1000.,wdsSPR->edges[i]->getDoubleProperty("roughness"),0.0,"Open");
    //if(wdsSPR->edges[i]->type == "PressurePoint")
    //  fprintf(wfile, " %-16s\t%-16s\t%-16s\t%-12.4f\t%-12.4f\t%-12.4f\t%-12.4f\t%-6s;\n",("PIPE_" + wdsSPR->edges[i]->getEdgeStringProperty("name")).c_str(),wdsSPR->edges[i]->getEdgeStringProperty("name").c_str(),wdsSPR->edges[i]->startNodeName.c_str(),0.0,1.0,0.0,0.0,"Open");
  }
  // TODO: Pump curves

  fprintf(wfile, "\n[PUMPS]\n");
  fprintf(wfile, ";ID             \tNode1           \tNode2           \tParameters\n");
  int pump_counter=0;
  for(int i=0; i<wdsSPR->edges.size(); i++){
    if(wdsSPR->edges[i]->type == "Pump"){
      fprintf(wfile, " %-16s\t%-16s\t%-16s\t%-6s;\n", wdsSPR->edges[i]->getEdgeStringProperty("name").c_str(), wdsSPR->edges[i]->startNodeName.c_str(), wdsSPR->edges[i]->endNodeName.c_str(), ("HEAD "+to_string(pump_counter)).c_str());
      pump_counter++;
    }
  }

  // TODO: Different type of valves
  fprintf(wfile, "\n[VALVES]\n");
  fprintf(wfile, ";ID             \tNode1           \tNode2           \tDiameter    \tType\tSetting     \tMinorLoss   \n");
  for(int i=0; i<wdsSPR->edges.size(); i++){
    if(wdsSPR->edges[i]->typeCode == 9)
      fprintf(wfile, " %-16s\t%-16s\t%-16s\t%-12.4f\t%-4s\t%-12.4f\t%-12.4f\n", wdsSPR->edges[i]->getEdgeStringProperty("name").c_str(), wdsSPR->edges[i]->startNodeName.c_str(), wdsSPR->edges[i]->endNodeName.c_str(), 100., "ISO",0.0,0.0);
  }

  fprintf(wfile, "\n[TAGS]\n");

  fprintf(wfile, "\n[DEMANDS]\n");
  fprintf(wfile, ";Junction       \tDemand      \tPattern         \tCategory\n");

  // TODO: ADD STATUS i.e. pumps/valves are closed or open
  fprintf(wfile, "\n[STATUS]\n");
  fprintf(wfile, ";ID              \tStatus/Setting\n");

  // TODO: ADD PATTERNS
  fprintf(wfile, "\n[PATTERNS]\n");
  fprintf(wfile, ";ID             \tMultipliers\n");

  // TODO: PUMP CURVES
  fprintf(wfile, "\n[CURVES]\n");
  fprintf(wfile, ";ID             \tX-Value     \tY-Value\n");
  pump_counter = 0;
  for(int i=0; i<wdsSPR->edges.size(); i++){
    if(wdsSPR->edges[i]->type == "Pump"){
      vector <double> pumpCurveVolumeFlowRate = wdsSPR->edges[i]->getVectorProperty("volumeFlowRate");
      vector <double> pumpCurveHead = wdsSPR->edges[i]->getVectorProperty("head");
      for(int j=0; j<pumpCurveVolumeFlowRate.size(); j++)
        fprintf(wfile, " %-16i\t%-12.4f\t%-12.4f\n", pump_counter, pumpCurveVolumeFlowRate[j], pumpCurveHead[j]);
      pump_counter++;
    }
  }

  fprintf(wfile, "\n[CONTROLS]\n");

  fprintf(wfile, "\n[RULES]\n");

  fprintf(wfile, "\n[ENERGY]\n");

  fprintf(wfile, "\n[EMITTERS]\n");

  fprintf(wfile, "\n[QUALITY]\n");

  fprintf(wfile, "\n[SOURCES]\n");

  fprintf(wfile, "\n[REACTIONS]\n");

  fprintf(wfile, "\n[MIXING]\n");

  // TODO: for extended calculations
  fprintf(wfile, "\n[TIMES]\n");

  fprintf(wfile, "\n[REPORT]\n");

  fprintf(wfile, "\n[OPTIONS]\n");
  fprintf(wfile, " Units              \tLPS\n");
  fprintf(wfile, " Headloss           \tC-F\n");
  fprintf(wfile, " Specific Gravity   \t1.0\n");
  fprintf(wfile, " Viscosity          \t1.0\n");
  fprintf(wfile, " Trials             \t40\n");
  fprintf(wfile, " Accuracy           \t0.001\n");
  fprintf(wfile, " CHECKFREQ          \t2\n");
  fprintf(wfile, " MAXCHECK           \t10\n");
  fprintf(wfile, " DAMPLIMIT          \t0\n");
  fprintf(wfile, " Unbalanced         \tContinue 10\n");
  fprintf(wfile, " Pattern            \t1\n");
  fprintf(wfile, " Demand Multiplier  \t1.0\n");
  fprintf(wfile, " Emitter Exponent   \t0.5\n");
  fprintf(wfile, " Quality            \tChlorine mg/L\n");
  fprintf(wfile, " Diffusivity        \t1.0\n");
  fprintf(wfile, " Tolerance          \t0.01\n");

  // TODO: xpos or xcoord??
  fprintf(wfile, "\n[COORDINATES]\n");
  fprintf(wfile, ";Node           \tX-Coord         \tY-Coord\n");
  for(int i=0; i<wdsSPR->nodes.size(); i++)
  {
    fprintf(wfile, " %-16s\t%-16.4f\t%-16.4f\n", wdsSPR->nodes[i]->name.c_str(), wdsSPR->nodes[i]->xPosition, wdsSPR->nodes[i]->yPosition);
  }

  fprintf(wfile, "\n[VERTICES]\n");
  fprintf(wfile, "\n[LABELS]\n");
  fprintf(wfile, "\n[BACKDROP]\n");

  // End of writing the file
  fprintf(wfile, "\n[END]\n");
  fclose(wfile);

  cout << endl << "[*] File written successfully: " << file_name << endl;

  string case_name_check = case_name.substr(0,case_name.length()-4) + ".inp";
  cout << endl << "[!] Checking the case: " << case_name_check << endl;

  wdsINP = new HydraulicSolver(case_folder + case_name_check);

  bool hydOk = wdsINP->solveSystem();
    
  if(hydOk)
  {
    cout << "[!] Hydraulics ok :)";
  }
  else
  {
    cout << "[!] Hydraulics NOT ok :(";
  }
  int warningCounter = 0;

  // checking the rank of the wdn at iso valve nodes
  vector<int> checkRank(wdsINP->nodes.size(),0);
  for(int i=0; i<wdsINP->edges.size(); i++)
  {
    checkRank[wdsINP->edges[i]->startNodeIndex]++;
    if(wdsINP->edges[i]->numberNode==2)
      checkRank[wdsINP->edges[i]->endNodeIndex]++;
  }
  for(int i=0; i<wdsINP->nodes.size(); i++)
  {
    bool isISO=false;
    for(int j=0; j<wdsINP->nodes[i]->edgeIn.size(); j++)
    {
      if(wdsINP->edges[wdsINP->nodes[i]->edgeIn[j]]->typeCode == 9)
        isISO = true;
    }
    for(int j=0; j<wdsINP->nodes[i]->edgeOut.size(); j++)
    {
      if(wdsINP->edges[wdsINP->nodes[i]->edgeOut[j]]->typeCode == 9)
        isISO = true;
    }
    if(isISO && checkRank[i]!=2)
    {
      cout << endl << "!WARNING! Node " << wdsINP->nodes[i]->name << " has a rank of "  << checkRank[i] << endl;
      warningCounter++;
    }
  }

  // checking whether there are two iso valve directly parallel
  for(int i=0; i<wdsINP->edges.size(); i++)
  {
    if(wdsINP->edges[i]->typeCode == 9)
    {
      bool gotProblem = false;
      int si = wdsINP->edges[i]->startNodeIndex;
      int ei = wdsINP->edges[i]->endNodeIndex;
      for(int j=0; j<wdsINP->nodes[si]->edgeIn.size(); j++)
      { 
        int idx = wdsINP->nodes[si]->edgeIn[j];
        int tc = wdsINP->edges[idx]->typeCode;
        if((tc == -1 || tc == -2 || tc == 9) && idx != i)
        {
          gotProblem = true;
          break;
        }
      }
      for(int j=0; j<wdsINP->nodes[si]->edgeOut.size(); j++)
      { 
        int idx = wdsINP->nodes[si]->edgeOut[j];
        int tc = wdsINP->edges[idx]->typeCode;
        if((tc == -1 || tc == -2 || tc == 9) && idx != i)
        {
          gotProblem = true;
          break;
        }
      }
      for(int j=0; j<wdsINP->nodes[ei]->edgeIn.size(); j++)
      { 
        int idx = wdsINP->nodes[ei]->edgeIn[j];
        int tc = wdsINP->edges[idx]->typeCode;
        if((tc == -1 || tc == -2 || tc == 9) && idx != i)
        {
          gotProblem = true;
          break;
        }
      }
      for(int j=0; j<wdsINP->nodes[ei]->edgeOut.size(); j++)
      { 
        int idx = wdsINP->nodes[si]->edgeOut[j];
        int tc = wdsINP->edges[idx]->typeCode;
        if((tc == -1 || tc == -2 || tc == 9) && idx != i)
        {
          gotProblem = true;
          break;
        }
      }
      if(gotProblem)
      {
        cout << endl << " !WARNING! ISO valve " << wdsINP->edges[i]->name << " is in parallel with pres/pool/iso edges" << endl << " Do something about it..." << endl;
        warningCounter++;
      }
    }
  }

  if(warningCounter == 0)
  {
    cout << endl << "[!] Everything is ok :)" << endl;
  }
  else
  {
    cout << endl << "[!] There was " << warningCounter << " number of warning! :(" << endl;
  }

  cout << endl << endl;
  return 0;
}