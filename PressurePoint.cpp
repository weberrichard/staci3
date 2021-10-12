#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include "Edge.h"
#include "PressurePoint.h"

PressurePoint::PressurePoint(const string a_name, const double a_referenceCrossSection, const string a_startNodeName, const double a_density, const double a_head, const double a_totalHead, const double a_volumeFlowRate) :
    Edge(a_name, a_referenceCrossSection, a_volumeFlowRate, a_density){
    type = "PressurePoint";
    numberNode = 1;

    endNodeName = "<none>";
    startNodeName = a_startNodeName;
    head = a_head;
    headNominal = a_totalHead;
    status = 1;
    typeCode = -2;
}

//--------------------------------------------------------------
PressurePoint::~PressurePoint(){
}

//--------------------------------------------------------------
string PressurePoint::info(){
  ostringstream strstrm;
  strstrm << Edge::info();
  strstrm << "\n connection            : " << startNodeName << "(index:" << startNodeIndex << ")";
  strstrm << "\n head                   : " << head << " [m]" << endl;
  return strstrm.str();
}

//--------------------------------------------------------------
double PressurePoint::function(const VectorXd &ppq, VectorXd &fDer)
{ 
  double out; 
  if(status == 1) // OPEN
  {
    out = -ppq(0) + head; 
    fDer(0) = -1.0;    
  }
  else // CLOSED, status is 0 or -1
  {
    //out = ppq(2);
    //fDer(2) = 1.0;

    out = -ppq(0) + head + 1e8 * ppq(2);
    fDer(0) = -1.0;
    fDer(2) =  1e8;
  } 
  return out;
}

//--------------------------------------------------------------
void PressurePoint::initialization(int mode, double value){
  if (mode == 0)
    volumeFlowRate = 0.01;
  else
    volumeFlowRate = value; 
}


//--------------------------------------------------------------
double PressurePoint::getDoubleProperty(string prop){
  double out = 0.;
  if(prop == "head")
    out = head;
  else if(prop == "headNominal")
    out = headNominal;
  else if(prop == "pressure")
    out = head * gravity * density; 
  else if(prop == "massFlowRate") 
    out = volumeFlowRate * density/1000.;
  else if(prop == "volumeFlowRate")
    out = volumeFlowRate;
  else if(prop == "velocity")
    out = volumeFlowRate / 1000. / referenceCrossSection;
  else if(prop == "density")
    out = density; 
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section" || prop == "Aref")
    out = referenceCrossSection;
  else if(prop == "startHeight")
    out = startHeight;
  else if(prop == "waterAge")
    out = waterAge;
  else if(prop == "addedChlorine")
    out = addedChlorine;
  else if(prop == "microbesWater" || prop == "Cf")
    out = Cf;
  else if(prop == "microbesWall" || prop == "Cb")
    out = Cb;
  else if(prop == "substratWater" || prop == "Sf")
    out = Sf;
  else if(prop == "substratWall" || prop == "Sb") 
    out = Sb;
  else 
  {
    cout << endl << endl << "DOUBLE PressurePoint::getDoubleProperty() wrong argument:" << prop;
    cout << ", right values: head | pressure | volumeFlowRate | velocity | density | referenceCrossSection | user1 | user2" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void PressurePoint::setDoubleProperty(string prop, double value){
  if(prop == "head")
    head = value;
  else if(prop == "headNominal")
    headNominal = value;
  else if(prop == "pressure")
    head = value / density / gravity;
  else if(prop == "volumeFlowRate") 
    volumeFlowRate = value;
  else if(prop == "density")
    density = value;
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section")
    referenceCrossSection = value;
  else if(prop == "startHeight")
    startHeight = value;
  else if(prop == "addedChlorine" || prop == "addedAmmountOfChlorine")
    addedChlorine = value;
  else if(prop == "microbesWater" || prop == "Cf")
    Cf = value;
  else if(prop == "microbesWall" || prop == "Cb")
    Cb = value;
  else if(prop == "substratWater" || prop == "Sf")
    Sf = value;
  else if(prop == "substratWall" || prop == "Sb")
    Sb = value;
  else
  {  
    cout << endl << endl << "PressurePoint::setDoubleProperty( DOUBLE ) wrong argument:" << prop;
    cout << ", right values: head | volumeFlowRate | velocity | density | referenceCrossSection" << endl << endl;
  }
}

//--------------------------------------------------------------
string PressurePoint::getStringProperty(string prop){
  string out = "";
  if(prop == "patternID")
    out = patternID;
  else
  {
    cout << endl << endl << "STRING PressurePoint::getStringProperty() wrong argument:" << prop;
    cout << ", right values: patternID" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void PressurePoint::setStringProperty(string prop, string value){
  if(prop == "patternID")
    patternID = value;
  else
  {  
    cout << endl << endl << "PressurePoint::setProperty( STRING ) wrong argument:" << prop;
    cout << ", right values: patternID" << endl << endl;
  }
}

//--------------------------------------------------------------
int PressurePoint::getIntProperty(string prop){
  int out = 0;
  if(prop == "patternIndex")
    out = patternIndex;
  else
  {
    cout << endl << endl << "INT PressurePoint::getIntProperty() wrong argument:" << prop;
    cout << ", right values: patternIndex" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void PressurePoint::setIntProperty(string prop, int value){
  if(prop == "patternIndex")
    patternIndex = value;
  else
  {  
    cout << endl << endl << "PressurePoint::setProperty( STRING ) wrong argument:" << prop;
    cout << ", right values: patternIndex" << endl << endl;
  }
}