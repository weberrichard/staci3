#include "ValveISO.h"

ValveISO::ValveISO(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_referenceCrossSection, const double a_volumeFlowRate) : Valve(a_name, a_startNodeName, a_endNodeName, a_density,a_referenceCrossSection, a_volumeFlowRate) {
  type = "ValveISO";
  status = 1;
  typeCode = 9;
  setting = 0.0 * density/2./(referenceCrossSection*referenceCrossSection); // loss coefficient
}

//--------------------------------------------------------------
ValveISO::~ValveISO() {
}

//--------------------------------------------------------------
string ValveISO::info() {
  ostringstream ss;
  ss << Edge::info();
  ss << "\n segment from          : " << startSegment;
  ss << "\n segment to            : " << endSegment;
  ss << "\n connection            : " << startNodeName << " (index:" << startNodeIndex << ") --> " << endNodeName << " (index:" << endNodeIndex << ")";
  return ss.str();
}

//--------------------------------------------------------------
double ValveISO::function(const VectorXd &ppq, VectorXd &fDer)
{
  double out;
  if(status == 1) // OPEN
  {
    out = ppq(1) - ppq(0) + (endHeight-startHeight) + setting * ppq(2) * abs(ppq(2));
    fDer(0) = -1.0;
    fDer(1) =  1.0;
    fDer(2) = setting * abs(ppq(2));
  }
  else // CLOSED, status is 0 or -1
  {
    //out = ppq(2);
    //fDer(2) =  1.0;

    double k = 10.;
    out = ppq(1) - ppq(0) + (endHeight-startHeight) + k * ppq(2) * abs(ppq(2)) + 1e8 * ppq(2);
    fDer(0) = -1.0;
    fDer(1) =  1.0;
    fDer(2) = 2 * k * abs(ppq(2)) + 1e8;
  }
  return out;
}

//--------------------------------------------------------------
void ValveISO::initialization(int mode, double value) {
  if (mode == 0)
    volumeFlowRate = 1.;
  else
    volumeFlowRate = value;
}

//--------------------------------------------------------------
void ValveISO::setDoubleProperty(string property, double value) {
  if(property == "startHeight")
      startHeight = value;
  else if(property == "endHeight")
      endHeight = value;
  else {
    cout << endl << endl << "ERROR! ValveISO::setDoubleProperty(property), unkown property: property=" << property << endl << endl;
  }
}

//--------------------------------------------------------------
void ValveISO::setIntProperty(string property, int value) {
  if(property == "startSegment")
    startSegment = value;
  else if(property == "endSegment")
    endSegment = value;
  else
    cout << endl << "ERROR! ValveISO::setIntProperty(property), unkown property: property=" << property << endl << endl;
}

//--------------------------------------------------------------
double ValveISO::getDoubleProperty(string property) {
  double out = 0.0;
  if (property == "volumeFlowRate")
    out = volumeFlowRate;
  else if ((property == "length") || (property == "L"))
    out = 0.5;
  else if (property == "velocity")
    out = volumeFlowRate / referenceCrossSection;
  else if (property == "cross_section")
    out = referenceCrossSection;
  else if (property == "startHeight")
    out = startHeight;
  else if (property == "endHeight")
    out = endHeight;
  else
    cout << endl << "ERROR! ValveISO::getDoubleProperty(property), unkown property: property=" << property << endl << endl;

  return out;
}

//--------------------------------------------------------------
string ValveISO::getStringProperty(string property) {
  string out = "";
  if (property == "type")
    out = type;
  else
    cout << endl << "ERROR! ValveISO::getStringProperty(property), unkown property: property=" << property << endl << endl;
  return out;
}

//--------------------------------------------------------------
int ValveISO::getIntProperty(string property) {
  int out = 0;
  if(property == "startSegment")
    out = startSegment;
  else if(property == "endSegment")
    out = endSegment;
  else
    cout << endl << "ERROR! ValveISO::getIntProperty(property), unkown property: property=" << property << endl << endl;
  return out;
}