#include "Pool.h"

Pool::Pool(const string a_name, const string a_startNodeName, const double a_density, const double a_referenceCrossSection, const double a_bottomLevel, const double a_waterLevel, const double a_volumeFlowRate) :
    Edge(a_name, a_referenceCrossSection, a_volumeFlowRate, a_density) {
    type = "Pool";
    numberNode = 1;
    startNodeName = a_startNodeName;
    endNodeName = "<none>";
    bottomLevel = a_bottomLevel;
    waterLevel = a_waterLevel;
    initLevel = a_waterLevel;
    status = 1;
    typeCode = -1;
}

//--------------------------------------------------------------
Pool::~Pool() {
}

//--------------------------------------------------------------
string Pool::info() {
  ostringstream strstrm;
  strstrm << Edge::info();
  strstrm << "\n type                  : " << type;
  strstrm << "\n connection            : " << startNodeName << "(index:" << startNodeIndex << ")";
  strstrm << "\n bottom level          : " << bottomLevel << " [m]";
  strstrm << "\n water level           : " << waterLevel << " [m]" << endl;
  return strstrm.str();
}

//--------------------------------------------------------------
double Pool::function(const VectorXd &ppq, VectorXd &fDer)
{
  double out;
  if(status == 1) // open
  {
    out = -ppq(0) - startHeight + (bottomLevel + waterLevel);
    fDer(0) = -1.0;
  }
  else // closed, status is 0 or -1
  {
    //out = ppq(2);
    //fDer(2) = 1.0;

    out = -ppq(0) - startHeight + (bottomLevel + waterLevel) + 1e8 * ppq(2);
    fDer(0) = -1.0;
    fDer(2) =  1e8;
  }
  return out;
}

//--------------------------------------------------------------
void Pool::initialization(int mode, double value)
{
  if (mode == 0)
    volumeFlowRate = 1;
  else
    volumeFlowRate = value;
}

//--------------------------------------------------------------
double Pool::getDoubleProperty(string prop)
{
  double out = 0.;
  if(prop == "waterLevel")
    out = waterLevel;
  else if(prop == "bottomLevel")
    out = bottomLevel;
  else if(prop == "minLevel")
    out = minLevel;
  else if(prop == "maxLevel")
    out = maxLevel;
  else if(prop == "initLevel")
    out = initLevel;
  else if(prop == "head")
    out = (bottomLevel + waterLevel);
  else if(prop == "pressure")
    out = (bottomLevel + waterLevel) * density * gravity;
  else if(prop == "massFlowRate")
    out = volumeFlowRate * density/1000.;
  else if(prop == "volumeFlowRate")
    out = volumeFlowRate;
  else if(prop == "velocity")
    out = volumeFlowRate / 1000. / referenceCrossSection;
  else if(prop == "density")
    out = density;
  else if(prop == "referenceCrossSection" || prop == "Aref")
    out = referenceCrossSection;
  else if(prop == "waterAge")
    out = waterAge;
  else if(prop == "chlorine")
    out = chlorine;
  else
  {
    cout << endl << endl << "DOUBLE Pool::getDoubleProperty() wrong argument:" << prop;
    cout << ", right values: waterLevel | bottomLevel | head | pressure | massFlowRate | velocity | density | referenceCrossSection | user1 | user2" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void Pool::setDoubleProperty(string prop, double value)
{
  if(prop == "waterLevel")
    waterLevel = value;
  else if(prop == "bottomLevel")
    bottomLevel = value;
  else if(prop == "minLevel")
    minLevel = value;
  else if(prop == "maxLevel") 
    maxLevel = value;
  else if(prop == "initLevel")
    initLevel = value;
  else if(prop == "volumeFlowRate")
    volumeFlowRate = value;
  else if(prop == "density")
    density = value;
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section")
    referenceCrossSection = value;
  else if(prop == "startHeight")
    startHeight = value;
  else if(prop == "waterAge")
    waterAge = value;
  else if(prop == "chlorine")
    chlorine = value;
  else
  {  
    cout << endl << endl << "Pool::setProperty( DOUBLE ) wrong argument:" << prop;
    cout << ", right values: head | volumeFlowRate | velocity | density | referenceCrossSection | user1 | user2" << endl << endl;
  }
}

//--------------------------------------------------------------
bool Pool::getBoolProperty(string prop)
{
  double out = false;
  if(prop == "doOverflow")
  {
    out = doOverflow;
  }
  else
  {  
    cout << endl << endl << "Pool::setBoolProperty( BOOL ) wrong argument:" << prop;
    cout << ", right values: doOverflow" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void Pool::setBoolProperty(string prop, bool value)
{
  if(prop == "doOverflow")
  {
    doOverflow = value;
  }
  else
  {
    cout << endl << endl << "BOOL Pool::getBoolProperty() wrong argument:" << prop;
    cout << ", right values: doOverflow" << endl << endl;
  }
}

