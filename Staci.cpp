#include "Staci.h"

using namespace Eigen;

Staci::Staci(string fileName) {
  definitionFile = fileName;

  // getting rid of global path
  caseName = definitionFile.substr(definitionFile.rfind('/')+1);
  // getting rid of extension
  caseName = caseName.substr(0,caseName.length()-4);

  string fileFormat = definitionFile.substr(definitionFile.length()-3,3); // SPR or INP
  if(fileFormat == "inp") // Standard EPANET inp format
  {
    loadSystem();
  }
  else
  {
    cout << endl << "Unkown file format: " << fileFormat << endl << "Available file formats are: inp" << endl;
    exit(-1);
  }

  // calculating the number of nodes and edges
  numberEdges = edges.size();
  numberNodes = nodes.size();
  
  // Finding the indicies of nodes for the edges and vica versa
  buildSystem();
}

//--------------------------------------------------------------
Staci::~Staci() { 
}

//--------------------------------------------------------------
void Staci::buildSystem(){
  bool startGotIt = false, endGotIt = false;
  int j = 0;
  int startNode = -1, endNode = -1;

  // Clearing the in/out going edges from nodes
  for(int i=0; i<numberNodes; i++){
    nodes[i]->edgeIn.clear();
    nodes[i]->edgeOut.clear();
  }

  for(int i=0; i<numberEdges; i++)
  {
    startGotIt = false;
    j = 0;
    while((j < numberNodes) && (!startGotIt)){
      if((edges[i]->getEdgeStringProperty("startNodeName")).compare(nodes[j]->getName()) == 0){
        startGotIt = true;
        startNode = j;
        nodes[j]->edgeOut.push_back(i);
      }
      j++;
    }
    if(!startGotIt){
      cout << "\n!!! ERROR !!! Edge name: " << edges[i]->getEdgeStringProperty("name").c_str() << ", startnode not found !!!";
      cout << endl << "startNode: " << startNode << endl << "Exiting..." << endl;
      exit(-1);
    } 

    if(edges[i]->getEdgeIntProperty("numberNode") == 2){
      endGotIt = false;
      j = 0;
      while ((j < numberNodes) && (!endGotIt)) {
        if ((edges[i]->getEdgeStringProperty("endNodeName")).compare(nodes[j]->getName()) == 0) {
          endGotIt = true;
          endNode = j;
          nodes[j]->edgeIn.push_back(i);
        }
        j++;
      }
      if(!endGotIt) {
        cout << "\n!!! ERROR !!! Edge name: " << edges[i]->getEdgeStringProperty("name").c_str() << ", startnode not found !!!";
        cout << endl << "startNode: " << startNode << endl << "Exiting..." << endl;
        exit(-1);
      }
    }

    if(edges[i]->getEdgeIntProperty("numberNode") == 2){
      edges[i]->setEdgeIntProperty("startNodeIndex", startNode);
      edges[i]->setEdgeIntProperty("endNodeIndex", endNode);
    }
    else{
      edges[i]->setEdgeIntProperty("startNodeIndex", startNode);
      edges[i]->setEdgeIntProperty("endNodeIndex", -1);
    }
  }
}

//--------------------------------------------------------------
void Staci::listSystem() {
  cout << "\n\n Nodes:\n--------------------------";
  for (int i = 0; i < numberNodes; i++)
    cout << nodes[i]->info(true);
  cout << "\n\n Edges:\n--------------------------";
  for (int i = 0; i < numberEdges; i++)
    cout << edges[i]->info();
}


/*! WR: find the the given node IDs and gives back the Indicies
  *id:  vector containing the IDs*/
vector<int> Staci::ID2Index(const vector<string> &id){
  int n_id = id.size();
  bool gotIt = false;
  vector<int> idx(n_id);
  for(int j=0;j<n_id;j++){
    int i=0;
    gotIt = false;
    while(!gotIt && i<numberNodes){
      if(id[j] == nodes[i]->getName()){
        idx[j] = i;
        gotIt = true;
      }
      i++;
    }
    if(gotIt == false)
      cout << "\n!!!WARNING!!!\nStaci:ID2Indicies function\nNode is not existing, id: " << id[j] << endl<< "\nContinouing...";
  }
  return idx;
}

//--------------------------------------------------------------
void Staci::checkSystem()
{
  ostringstream msg1;
  bool stop = false;

  cout << endl << " [*] Checking System  ";

  cout << "\n [*] Looking for identical IDs  ";
  string name1, name2;
  for(int i = 0; i < numberEdges; i++){
    name1 = edges.at(i)->getEdgeStringProperty("name");
    for(int j = 0; j < numberEdges; j++){
      name2 = edges.at(j)->getEdgeStringProperty("name");
      if(i != j){
        if(name1 == name2){
          cout << "\n !!!ERROR!!! edge #" << i << ": " << name1 << " and edge #" << j << ": " << name2 << " with same ID!" << endl;
          stop = true;
        }
      }
    }

    for(int j = 0; j < numberNodes; j++){
      name2 = nodes.at(j)->getName();
      if(i != j){
        if(name1 == name2){
          cout << "\n !!!ERROR!!! edge #" << i << ": " << name1 << " and node #" << j << ": " << name2 << " with same ID!" << endl;
          stop = true;
        }
      }
    }
  }

  if (stop)
    exit(-1);
  else
    cout << endl << " [*] Checking System:  OK";
}

//--------------------------------------------------------------
void Staci::saveResult(string property, string element){

  vector<string> allElement{"All","Node","Pipe","Pump","PressurePoint","Pool"};

  bool elementExist = false;
  int i=0;
  while(i<allElement.size() && !elementExist){
    if(element == allElement[i])
      elementExist = true;
    i++;
  }

  if(elementExist)
  { 
    mkdir(caseName.c_str(),0777);
    if(element == "All")
    {
      remove((caseName + "/Node.txt").c_str());
      remove((caseName + "/Pipe.txt").c_str());
      remove((caseName + "/Pump.txt").c_str());
      remove((caseName + "/Pool.txt").c_str());
      remove((caseName + "/Pres.txt").c_str());
    }

    ofstream wfile;
    if(element == "Node" || element == "All")
    { 
      remove((caseName + "/Node.txt").c_str());
      wfile.open(caseName + "/Node.txt");
      for(int i=0; i<numberNodes; i++)
        wfile << nodes[i]->getProperty(property) << endl;
      wfile.close();
    }

    if(element == "Pipe" || element == "All")
    { 
      remove((caseName + "/Pipe.txt").c_str());
      wfile.open(caseName + "/Pipe.txt");
      for(int i=0; i<numberEdges; i++){
        if(edges[i]->type == "Pipe"){
          wfile << edges[i]->getDoubleProperty(property) << endl;
        }
      }
      wfile.close();
    }

    if(element == "Pump" || element == "All")
    {
      remove((caseName + "/Pump.txt").c_str());
      wfile.open(caseName + "/Pump.txt");
      for(int i=0; i<numberEdges; i++){
        if(edges[i]->type == "Pump"){
          wfile << edges[i]->getEdgeDoubleProperty(property) << endl;
        }
      }  wfile.close();
    }

    if(element == "PressurePoint" || element == "All")
    {
      remove((caseName + "/Pres.txt").c_str());
      wfile.open(caseName + "/Pres.txt");
      for(int i=0; i<numberEdges; i++){
        if(edges[i]->type == "PressurePoint"){
          wfile << edges[i]->getEdgeDoubleProperty(property) << endl;
        }
      }
      wfile.close();
    }

    if(element == "Pool" || element == "All")
    { 
      remove((caseName + "/Pool.txt").c_str());
      wfile.open(caseName + "/Pool.txt");
      for(int i=0; i<numberEdges; i++){
        if(edges[i]->type == "Pool"){
          wfile << edges[i]->getEdgeDoubleProperty(property) << endl;
        }
      }
      wfile.close();
    }
  }
  else
  {
    cout << endl << "Elemenet: " << element << " does not exist in Staci::saveResult() function" << endl << "Possible elements: ";
    for(int i=0; i<allElement.size(); i++)
      cout << allElement[i] << ", ";
    cout << endl;
  }
}

//--------------------------------------------------------------
/*void Staci::changeEdgeStatus(int idx, bool state){

  if(state == false && !edges[idx]->isClosed) // Closing an edge
  {
    edges[idx]->isClosed = true;
    int nodeFrom = edges[idx]->getEdgeIntProperty("startNodeIndex");
    for(int i=0; i<nodes[nodeFrom]->edgeOut.size(); i++)
      if(nodes[nodeFrom]->edgeOut[i] == idx)
        nodes[nodeFrom]->edgeOut.erase(nodes[nodeFrom]->edgeOut.begin() + i);

    if(nodes[nodeFrom]->edgeOut.size() + nodes[nodeFrom]->edgeIn.size() == 0)
      nodes[nodeFrom]->isClosed = true;

    if(edges[idx]->getEdgeIntProperty("numberNode") == 2){
      int nodeTo = edges[idx]->getEdgeIntProperty("endNodeIndex");
      for(int i=0; i<nodes[nodeTo]->edgeIn.size(); i++)
        if(nodes[nodeTo]->edgeIn[i] == idx)
          nodes[nodeTo]->edgeIn.erase(nodes[nodeTo]->edgeIn.begin() + i);

      if(nodes[nodeTo]->edgeOut.size() + nodes[nodeTo]->edgeIn.size() == 0)
        nodes[nodeTo]->isClosed = true;
    }
  }
  else if(state == true && edges[idx]->isClosed)
  {
    edges[idx]->isClosed = false;
    int nodeFrom = edges[idx]->getEdgeIntProperty("startNodeIndex");
    nodes[nodeFrom]->edgeOut.push_back(idx);
    if(nodes[nodeFrom]->edgeOut.size() + nodes[nodeFrom]->edgeIn.size() == 1)
      nodes[nodeFrom]->isClosed = false;

    if(edges[idx]->getEdgeIntProperty("numberNode") == 2){
      int nodeTo = edges[idx]->getEdgeIntProperty("endNodeIndex");
      nodes[nodeTo]->edgeIn.push_back(idx);
      if(nodes[nodeTo]->edgeOut.size() + nodes[nodeTo]->edgeIn.size() == 1)
        nodes[nodeTo]->isClosed = false;
    }
  }
}*/

//--------------------------------------------------------------
/*void Staci::changeEdgeStatus(string ID, bool state){

  int i=0, idx=-1;
  bool gotIt=false;
  while(i<numberEdges && !gotIt)
  {
    if(ID.compare(edges[i]->getEdgeStringProperty("name")) == 0)
    {
      gotIt = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n!!!WARNING!!!\nStaci:changeEdgeStatus function\nEdge is not existing, id: " << ID << endl<< "\nContinouing...";
  }
  else
  {
    changeEdgeStatus(idx,state);
  }
}*/

//--------------------------------------------------------------
/*int Staci::getVectorIndex(const vector<int> &v, int value){
  int idx = -1;
  for(int i=0; i<v.size(); i++)
    if(v[i] == value){
      idx = i;
      break;
    }
  if(idx == -1)
    cout << endl << "!!! WARNING !!! Staci::getVectorIndex did not found value " << value << endl;

  return idx;
}*/

//--------------------------------------------------------------
int Staci::nodeIDtoIndex(string ID){
  int i=0, idx=-1;
  bool gotIt=false;
  while(i<numberNodes && !gotIt)
  {
    if(ID.compare(nodes[i]->getName()) == 0)
    {
      gotIt = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n!!!WARNING!!!\nStaci:nodeIDtoIndex function\nNode is not existing, ID: " << ID << endl<< "\nContinouing...";
  }
  return idx;
}

//--------------------------------------------------------------
int Staci::edgeIDtoIndex(string ID){
  int i=0, idx=-1;
  bool gotIt=false;
  while(i<numberEdges && !gotIt)
  {
    if(ID.compare(edges[i]->getEdgeStringProperty("name")) == 0)
    {
      gotIt = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n!!!WARNING!!!\nStaci:edgeIDtoIndex function\nNode is not existing, ID: " << ID << endl<< "\nContinouing...";
  }
  return idx;
}