#include "Staci.h"

using namespace Eigen;

Staci::Staci(string fileName)
{
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
  else if(fileFormat == "spr") // Old STACI spr format
  {
    IOxml IOxmlObj(definitionFile.c_str());
    IOxmlObj.loadSystem(nodes, edges);
    frictionModel = "C-F";
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

  // Fill up the index vectors to make the code more efficient
  buildIndexing();
}

//--------------------------------------------------------------
Staci::~Staci(){}

//--------------------------------------------------------------
void Staci::buildSystem()
{
  bool startGotIt = false, endGotIt = false;
  int j = 0;
  int startNode = -1, endNode = -1;

  // Clearing the in/out going edges from nodes
  for(int i=0; i<numberNodes; i++)
  {
    nodes[i]->edgeIn.clear();
    nodes[i]->edgeOut.clear();
  }

  for(int i=0; i<numberEdges; i++)
  {
    startGotIt = false;
    j = 0;
    while((j < numberNodes) && (!startGotIt)){
      if((edges[i]->startNodeName).compare(nodes[j]->name) == 0){
        startGotIt = true;
        startNode = j;
        nodes[j]->edgeOut.push_back(i);
      }
      j++;
    }
    if(!startGotIt){
      cout << "\n!!! ERROR !!! Edge name: " << edges[i]->name.c_str() << ", startnode not found !!!";
      cout << endl << "startNode: " << edges[i]->startNodeName.c_str() << endl << "Exiting..." << endl;
      exit(-1);
    } 

    if(edges[i]->numberNode == 2){
      endGotIt = false;
      j = 0;
      while ((j < numberNodes) && (!endGotIt)) {
        if ((edges[i]->endNodeName).compare(nodes[j]->name) == 0) {
          endGotIt = true;
          endNode = j;
          nodes[j]->edgeIn.push_back(i);
        }
        j++;
      }
      if(!endGotIt) {
        cout << "\n!!! ERROR !!! Edge name: " << edges[i]->name.c_str() << ", startnode not found !!!";
        cout << endl << "startNode: " << startNode << endl << "Exiting..." << endl;
        exit(-1);
      }
    }

    if(edges[i]->numberNode == 2){
      edges[i]->startNodeIndex = startNode;
      edges[i]->endNodeIndex = endNode;
    }
    else{
      edges[i]->startNodeIndex = startNode;
      edges[i]->endNodeIndex = -1;
    }
  }
}

//--------------------------------------------------------------
void Staci::buildIndexing()
{
  poolIndex.clear(); 
  presIndex.clear(); 
  pumpIndex.clear(); 
  valveIndex.clear();
  valveISOIndex.clear();
  pipeIndex.clear();
  pipeCVIndex.clear();
  flowMeterIndex.clear();

  for(int i=0; i<edges.size(); i++)
  { 
    int typeCode = edges[i]->typeCode;
    if(typeCode == 1) // normal pipe
    {
      pipeIndex.push_back(i);
    }
    else if(typeCode == 0) // pipe with cv
    {
      pipeIndex.push_back(i);
      pipeCVIndex.push_back(i);
    }
    else if(typeCode == 2) // pump
    {
      pumpIndex.push_back(i);
    }
    else if(typeCode == 3 || typeCode == 4 || typeCode == 5 || typeCode == 6 || typeCode == 7 || typeCode == 8) // active valves
    {
      valveIndex.push_back(i);
    }
    else if(typeCode == 9) // ISO valves
    {
      valveIndex.push_back(i);
      valveISOIndex.push_back(i);
    }
    else if(typeCode == 10) // Flow Meter
    {
      flowMeterIndex.push_back(i);
    }
    else if(typeCode == -1) // pool
    {
      poolIndex.push_back(i);
    }
    else if(typeCode == -2) // pressure point
    {
      presIndex.push_back(i);
    }
    else
    {
      cout << endl << "!WARNING! Unkown typeCode: " << typeCode << ", type: " << edges[i]->type << ", name: " << edges[i]->name << endl;
    }
  }
}


//--------------------------------------------------------------
void Staci::listSystem()
{
  cout << "\n\n Nodes:\n--------------------------";
  for (int i = 0; i < numberNodes; i++)
  {
    cout << "\n Node index      : " << i;
    cout << nodes[i]->info(true);
  }
  cout << "\n\n Edges:\n--------------------------";
  for (int i = 0; i < numberEdges; i++)
  {
    cout << "\n\n Edge index            : " << i;
    cout << edges[i]->info();
  }
}

//--------------------------------------------------------------
void Staci::listSystemShort()
{
   cout << "\n\n Nodes:\n--------------------------\n";
   for(int i=0; i<numberNodes; i++)
   {
      printf("%10s\n",nodes[i]->name.c_str());
   }
   cout << "\n\n Edges:\n--------------------------\n";
   for(int i=0; i<numberEdges; i++)
   {
      int is = edges[i]->startNodeIndex;
      int ie = edges[i]->endNodeIndex;
      if(ie == -1)
      {
         printf("%10s, %10s\n", edges[i]->name.c_str(), nodes[is]->name.c_str());
      }
      else
      {
         printf("%10s, %10s-%10s\n", edges[i]->name.c_str(), nodes[is]->name.c_str(), nodes[ie]->name.c_str());
      }
   }
   cout << endl;
}


//--------------------------------------------------------------
/*vector<int> Staci::ID2Index(const vector<string> &id){
  int n_id = id.size();
  bool gotIt = false;
  vector<int> idx(n_id);
  for(int j=0;j<n_id;j++){
    int i=0;
    gotIt = false;
    while(!gotIt && i<numberNodes){
      if(id[j] == nodes[i]->name){
        idx[j] = i;
        gotIt = true;
      }
      i++;
    }
    if(gotIt == false)
      cout << "\n!!!WARNING!!!\nStaci:ID2Indicies function\nNode is not existing, id: " << id[j] << endl << "Continouing..." << endl;
  }
  return idx;
}*/

//--------------------------------------------------------------
void Staci::checkSystem()
{
  ostringstream msg1;
  bool stop = false;

  cout << endl << " [*] Checking System  ";

  cout << "\n [*] Looking for identical IDs  ";
  string name1, name2;
  for(int i = 0; i < numberEdges; i++){
    name1 = edges.at(i)->name;
    for(int j = 0; j < numberEdges; j++){
      name2 = edges.at(j)->name;
      if(i != j){
        if(name1 == name2){
          cout << "\n !!!ERROR!!! edge #" << i << ": " << name1 << " and edge #" << j << ": " << name2 << " with same ID!" << endl;
          stop = true;
        }
      }
    }

    for(int j = 0; j < numberNodes; j++){
      name2 = nodes.at(j)->name;
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
void Staci::saveResult(string property, string element)
{
  vector<string> allElement{"All","Node","Pipe","Pump","PressurePoint","Pool","Valve"};

  bool elementExist = false;
  int i=0;
  while(i<allElement.size() && !elementExist)
  {
    if(element == allElement[i])
      elementExist = true;
    i++;
  }

  if(elementExist)
  {
    // LINUX
    //mkdir("Network Data",0777);
    //mkdir(("Network Data/" + caseName).c_str(),0777);

    // FOR WINDOWS
    //mkdir("Network Data");
    //mkdir(("Network Data/" + caseName).c_str());

    makeDirectory("Network Data");
    makeDirectory(("Network Data/" + caseName).c_str());

    if(element == "All")
    {
      remove(("/Network Data/" + caseName + "/Node.txt").c_str());
      remove(("/Network Data/" + caseName + "/Pipe.txt").c_str());
      remove(("/Network Data/" + caseName + "/Pump.txt").c_str());
      remove(("/Network Data/" + caseName + "/Pool.txt").c_str());
      remove(("/Network Data/" + caseName + "/Pres.txt").c_str());
      remove(("/Network Data/" + caseName + "/Valve.txt").c_str());
    }

    ofstream wfile;
    if(element == "Node" || element == "All")
    { 
      remove(("Network Data/" + caseName + "/Node.txt").c_str());
      wfile.open("Network Data/" + caseName + "/Node.txt");
      for(int i=0; i<numberNodes; i++)
        wfile << nodes[i]->getProperty(property) << endl;
      wfile.close();
    }

    if(element == "Pipe" || element == "All")
    { 
      remove(("Network Data/" + caseName + "/Pipe.txt").c_str());
      wfile.open("Network Data/" + caseName + "/Pipe.txt");
      for(int i=0; i<numberEdges; i++)
      {
        if(edges[i]->typeCode == 0 || edges[i]->typeCode == 1) // Pipe or PipeCV
        {
          wfile << edges[i]->getDoubleProperty(property) << endl;
        }
      }
      wfile.close();
    }

    if(element == "Pump" || element == "All")
    {
      remove(("Network Data/" + caseName + "/Pump.txt").c_str());
      wfile.open("Network Data/" + caseName + "/Pump.txt");
      for(int i=0; i<numberEdges; i++)
      {
        if(edges[i]->typeCode == 2) // Pump
        {
          wfile << edges[i]->getEdgeDoubleProperty(property) << endl;
        }
      }  wfile.close();
    }

    if(element == "Valve" || element == "All")
    {
      remove(("Network Data/" + caseName + "/Valve.txt").c_str());
      wfile.open("Network Data/" + caseName + "/Valve.txt");
      for(int i=0; i<numberEdges; i++)
      {
        if(edges[i]->typeCode == 9) // ValveISO TODO more valve type
        {
          wfile << edges[i]->getEdgeDoubleProperty(property) << endl;
        }
      }  wfile.close();
    }

    if(element == "PressurePoint" || element == "All")
    {
      remove(("Network Data/" + caseName + "/Pres.txt").c_str());
      wfile.open("Network Data/" + caseName + "/Pres.txt");
      for(int i=0; i<numberEdges; i++)
      {
        if(edges[i]->typeCode == -2)
        {
          wfile << edges[i]->getEdgeDoubleProperty(property) << endl;
        }
      }
      wfile.close();
    }

    if(element == "Pool" || element == "All")
    { 
      remove(("Network Data/" + caseName + "/Pool.txt").c_str());
      wfile.open("Network Data/" + caseName + "/Pool.txt");
      for(int i=0; i<numberEdges; i++)
      {
        if(edges[i]->typeCode == -1)
        {
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
int Staci::nodeIDtoIndex(string ID)
{
  int i=0, idx=-1;
  bool gotIt=false;
  while(i<numberNodes && !gotIt)
  {
    if(ID.compare(nodes[i]->name) == 0)
    {
      gotIt = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n!!!WARNING!!!\nStaci:nodeIDtoIndex function\nNode is not existing, ID: " << ID << "\nContinouing..." << endl;
  }
  return idx;
}

//--------------------------------------------------------------
int Staci::edgeIDtoIndex(string ID)
{
  int i=0, idx=-1;
  bool gotIt=false;
  while(i<numberEdges && !gotIt)
  {
    if(ID.compare(edges[i]->name) == 0)
    {
      gotIt = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n!!!WARNING!!!\nStaci:edgeIDtoIndex function\nNode is not existing, ID: " << ID << "\nContinouing..." << endl;
  }
  return idx;
}
