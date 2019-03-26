#include "Valve.h"

Valve::Valve(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_referenceCrossSection, vector<double> a_charX, vector<double> a_charY, double a_position, const double a_massFlowRate) :
        Edge(a_name, a_referenceCrossSection, a_massFlowRate, a_density) {
    type = "Valve";
    numberNode = 2;
    startNodeName = a_startNodeName;
    endNodeName = a_endNodeName;
    charX = a_charX;
    charY = a_charY;
    position = a_position;
    updateLoss();

}

//--------------------------------------------------------------
Valve::~Valve() {
}

//--------------------------------------------------------------
string Valve::info() {
  ostringstream strstrm;
  strstrm << Edge::info();
  strstrm << endl << "  connection : " << startNodeName << "(index:" << startNodeIndex << ") --> " << endNodeName << "(index:" << endNodeIndex << ")\n";
  cout << setprecision(3);
  vector<double>::iterator it;
  strstrm << "       data : charX [%]      = ";
  for (it = charX.begin(); it != charX.end(); it++)
    strstrm << *(it) << "  ";
  strstrm << "\n";
  strstrm << "                charY [-]   = ";
  for (it = charY.begin(); it != charY.end(); it++)
    strstrm << *(it) << "  ";
  strstrm << "\n";
  strstrm << endl << "       The actual position: " << fixed << position << endl<< "The actual loss coefficient: " << loss << endl;

  return strstrm.str();
}

//--------------------------------------------------------------
double Valve::function(vector<double> x) {
  double result;
  double pe = x[0] * density * gravity;
  double pv = x[1] * density * gravity;
  double he = x[2];
  double hv = x[3];

  result = (pv - pe) / density / gravity + (hv - he) + loss * massFlowRate * fabs(massFlowRate);

  headLoss = loss * massFlowRate * abs(massFlowRate);

  return result;
}

//--------------------------------------------------------------
vector<double> Valve::functionDerivative(vector<double> x) {
  vector<double> result;
  result.push_back(-1.0);
  result.push_back(+1.0);
  result.push_back(+2 * loss * abs(massFlowRate));
  result.push_back(0.0);

  return result;
}

//--------------------------------------------------------------
void Valve::initialization(int mode, double value) {
  if (mode == 0)
    massFlowRate = 1.;
  else
    massFlowRate = value;
}

//--------------------------------------------------------------
void Valve::setDoubleProperty(string property, double value) {
  if (property == "position") {
    position = value;
    updateLoss();
  } else {
    cout << endl << endl << "ERROR! Valve::setDoubleProperty(property), unkown property: property=" << property << endl << endl;
  }
}

//--------------------------------------------------------------
double Valve::getDoubleProperty(string property) {

  double out = 0.0;
  if (property == "position")
    out = position;
  else if (property == "loss")
    out = loss;
  else if (property == "massFlowRate")
    out = massFlowRate;
  else if (property == "headLoss")
    out = headLoss;
  else if (property == "headLoss_per_unit_length")
    out = headLoss;
  else if ((property == "length") || (property == "L"))
    out = 0.5;
  else if (property == "cross_section")
    out = referenceCrossSection;
  else
    cout << endl << "ERROR! Valve::getDoubleProperty(property), unkown property: property=" << property << endl << endl;

  return out;
}

//--------------------------------------------------------------
string Valve::getStringProperty(string property) {

  string out = "";
  if (property == "type")
    out = type;
  else
    cout << endl << "ERROR! Valve::getProperty(property), unkown property: property=" << property << endl << endl;
  return out;
}

//--------------------------------------------------------------
void Valve::updateLoss() {
  loss = interpolate(charX, charY, position);
  if (loss < 0.){
    loss = 0.;
    cout << endl << " !!!WARNING!!! Valve: loss coefficient is smaller than zero.\n Overriding with zero, and contionouing..." << endl;
  }
}