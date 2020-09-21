#include "Edge.h"

Edge::Edge(const string a_name, const double a_referenceCrossSection, const double a_volumeFlowRate, const double a_density){
  volumeFlowRate = a_volumeFlowRate; // l/s
  density = a_density; // kg/m3
  name = a_name;
  referenceCrossSection = a_referenceCrossSection;
  startNodeIndex = -1;
  endNodeIndex = -1;
  string startNodeName = "startNothingYet";
  endNodeName = "endNothingYet";

  if(abs(a_density) < 1.0e-3){
    ostringstream strstrm;
    strstrm.str("");
    strstrm << "\n\n******** ERROR *********";
    strstrm << "\n\telement name: " << name;
    strstrm << "\n\tmethod      : " << "Edge constructor";
    strstrm << "\n\tmessage     : " << "Density is to low." << "\n\n";
    cout << strstrm.str();
  exit(0);
  }
}

//--------------------------------------------------------------
Edge::~Edge(){}

//--------------------------------------------------------------
string Edge::info(){
  ostringstream strstrm;
  strstrm << "\n Edge name             : " << name;
  strstrm << "\n status                : " << status;
  strstrm << "\n type                  : " << type;
  strstrm << "\n open/closed           : ";
  if(status == 0)
    strstrm << "Closed";
  else if(status == 1)
    strstrm << "Open";
  else if(status == 2)
    strstrm << "Active";
  else
    strstrm << "Status: " << status;
  strstrm << "\n density               : " << density << " [kg/m3]";
  strstrm << "\n referenceCrossSection : " << referenceCrossSection << " [m^2]";
  strstrm << "\n volume flow rate      : " << volumeFlowRate << " [l/s]";
  strstrm << "\n segment               : " << segment;
  return strstrm.str();
}

//--------------------------------------------------------------
double Edge::getEdgeDoubleProperty(string prop)
{
  double out = 0.;
  if(prop == "volumeFlowRate")
    out = volumeFlowRate;
  else if(prop == "volumeFlowRateLPS")
    out = volumeFlowRate*1000.;
  else if(prop == "volumeFlowRateAbs")
    out = abs(volumeFlowRate);
  else if(prop == "volumeFlowRateAbsLPS")
    out = abs(volumeFlowRate)*1000.;
  else if(prop == "segment")
    out = (double)segment;
  else if(prop == "status")
    out = (double)status;
  else if(prop == "userOutput")
    out = userOutput;
  else if(prop == "massFlowRate") // kg/s
    out = volumeFlowRate * density;
  else if(prop == "velocity") // m/s
    out = volumeFlowRate  / referenceCrossSection;
  else if(prop == "density") // kg/m3
    out = density;
  else if(prop == "referenceCrossSection" || prop == "Aref")
    out = referenceCrossSection;
  else
  {
    cout << endl << endl << "DOUBLE Edge::getEdgeDoubleProperty() wrong argument:" << prop;
    cout << ", right values: volumeFlowRate | velocity | density | referenceCrossSection | segment" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
int Edge::getEdgeIntProperty(string prop){
  int out = 0;
  if(prop == "startNodeIndex" || prop == "start_node_index")
    out = startNodeIndex;
  else if(prop == "endNodeIndex" || prop == "end_node_index")
    out = endNodeIndex;
  else if(prop == "numberNode")
    out = numberNode;
  else if(prop == "segment")
    out = segment;
  else
  {
    cout << endl << endl << "INT Edge::getEdgeIntProperty() wrong argument:" << prop;
    cout << ", right values: startNodeIndex | endNodeIndex | numberNode" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
string Edge::getEdgeStringProperty(string prop)
{
  string out = "";
  if(prop == "startNodeName" || prop == "start_node_name")
    out = startNodeName;
  else if(prop == "endNodeName" || prop == "end_node_name")
    out = endNodeName;
  else if(prop == "name")
    out = name;
  else if(prop == "type")
    out = type;
  else
  {
    cout << endl << endl << "STRING Edge::getEdgeStringProperty() wrong argument:" << prop;
    cout << ", right values: startNodeName | endNodeName | name | type" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void Edge::setEdgeDoubleProperty(string prop, double value)
{
  if(prop == "volumeFlowRate")
    volumeFlowRate = value;
  else if(prop == "density")
    density = value;
  else if(prop == "userOutput")
    userOutput = value;
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section")
    referenceCrossSection = value;
  else
  {  
    cout << endl << endl << "Edge::setEdgeProperty( DOUBLE ) wrong argument:" << prop;
    cout << ", right values: volumeFlowRate | density | referenceCrossSection | userOutput" << endl << endl;
  }
}

//--------------------------------------------------------------
void Edge::setEdgeIntProperty(string prop, int value){
  if(prop == "startNodeIndex" || prop == "start_node_index")
    startNodeIndex = value;
  else if(prop == "endNodeIndex" || prop == "end_node_index")
    endNodeIndex = value;
  else if(prop == "numberNode")
    numberNode = value;
  else if(prop == "segment")
    segment = value;
  else
  {  
    cout << endl << endl << "Edge::setEdgeProperty( INT ) wrong argument:" << prop;
    cout << ", right values: startNodeIndex | endNodeIndex | numberNode" << endl << endl;
  }
}

//--------------------------------------------------------------
void Edge::setEdgeStringProperty(string prop, string value){
  if(prop == "startNodeName" || prop == "start_node_name")
    startNodeName = value;
  else if(prop == "endNodeName" || prop == "end_node_name")
    endNodeName = value;
  else if(prop == "name")
    name = value;
  else if(prop == "type")
    type = value;
  else
  {  
    cout << endl << endl << "Edge::setEdgeProperty( STRING ) wrong argument:" << prop;
    cout << ", right values: startNodeName | endNodeName | name | type" << endl << endl;
  }
}
