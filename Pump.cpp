#include "Pump.h"

using namespace std;
using namespace Eigen;

Pump::Pump(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_referenceCrossSection, vector<double> a_qCurve, vector<double> a_hCurve, const double a_pumpPerformance, const double a_volumeFlowRate, int a_pumpType) : Edge(a_name, a_referenceCrossSection, a_volumeFlowRate, a_density)
{
  type = "Pump";
  typeCode = 2;

  numberNode = 2;
  startNodeName = a_startNodeName;
  endNodeName = a_endNodeName;

  // for HEAD (Q-H curve)
  qCurve = a_qCurve;
  hCurve = a_hCurve;

  // for POWER (const performance)
  constantPerformance = a_pumpPerformance;

  // -1: epanet POWER type with constant performance
  //  0: epanet HEAD type with Q-H curve
  //  1: linear
  //  2: parabolic (least squares method)
  pumpType = a_pumpType;

  status = 1; // OPEN

  int numberPoint = qCurve.size();
  if(numberPoint != hCurve.size())
  {
    cout << endl << "!!!ERROR!!! Pump(" << name << ") characteristic curve: number of Q (" << numberPoint << ") is NOT equal to number of H (" << qCurve.size() << ") !" << endl << "Pump ID: " << a_name << endl << "Exiting..." << endl;
    exit(0);
  }

  if(pumpType == -1)
  {
    curveType = -1; // H*Q*rho*g = CONST

    constantPerformance /= (density*gravity); // H*Q = CONST2
  }
  else if(pumpType == 0) // EPANET HEAD STYLE
  {
    if(numberPoint == 1)
    {
      curveType = 0; // h(q) = a + bq^2

      double a = 4./3.*hCurve[0];
      double b = (hCurve[0] - a) / (qCurve[0]*qCurve[0]);
      coeffCurve.push_back(a);
      coeffCurve.push_back(b);
    }
    else if(numberPoint == 2) 
    {
      curveType = 1; // h(q) = aq + b
    }
    else if(numberPoint == 3)
    {
      curveType = 3; // h(q) = a + bq^c

      double c = log((hCurve[2]-hCurve[0])/(hCurve[1]-hCurve[0])) / log((qCurve[2]-qCurve[0])/(qCurve[1]-qCurve[0]));
      double b = (hCurve[1]-hCurve[0]) / pow(qCurve[1]-qCurve[0], c);
      double a = hCurve[0] - b*pow(qCurve[0], c);
      coeffCurve.push_back(a);
      coeffCurve.push_back(b);
      coeffCurve.push_back(c);
    }
    else if(numberPoint > 3) // changing to linear interpolation
    {
      curveType = 1;
    }
  }
  else if(pumpType == 1) // only linear interpolation
  {
    curveType = 1;
  }
  else if(pumpType == 2) // parabolic function fitting with least squares method
  {
    curveType = 2;

    // Copying vector<double> to Eigen::VectorXd
    VectorXd Q(numberPoint), H(numberPoint);
    for(int i=0; i<numberPoint; i++)
    {
      Q(i) = qCurve[i];
      H(i) = hCurve[i];
    }

    // Fitting second order polynomial to Q-H points with least squares method
    int curveOrder = 2;
    VectorXd C = leastSquaresPolynomial(Q,H,curveOrder);

    // Copying back to vector<double>
    for(int i=0; i<C.rows(); i++)
      coeffCurve.push_back(C(i));
  }
  else
  {
    cout << endl << "!ERROR! curve type: " << curveType << " unkown at pump " << name << endl << "Exiting..." << endl;
    exit(-1);
  }
}

//--------------------------------------------------------------
Pump::~Pump(){}

//--------------------------------------------------------------
string Pump::info()
{
  ostringstream strstrm;
  strstrm << Edge::info();
  strstrm << "\n head                  : " << getDoubleProperty("head") << " [m]";
  strstrm << "\n type                  : " << type;
  strstrm << "\n connection            : " << startNodeName << "(index:" << startNodeIndex << ") --> " << endNodeName << "(index:" << endNodeIndex << ")\n";
  strstrm << "\n pump type             : " << pumpType;
  strstrm << "\n curve type            : " << curveType;

  strstrm << endl;
  return strstrm.str();
}

//--------------------------------------------------------------
double Pump::function(const VectorXd &ppq, VectorXd &fDer)// ppq = [Pstart, Pend, VolumeFlowRate]
{
  double out;
  if(status == 1) // OPEN
  {
    double dp = characteristicCurve(ppq(2));
    out = (ppq(1) - ppq(0)) - dp + (endHeight - startHeight);

    double derivative;
    if(ppq(2)<0.){
      out = 1e8*ppq(2);
      derivative = 1e8;
    }
    else{
      double dx;
      if(abs(ppq(2))>0.1)
        dx = ppq(2)*0.001;
      else
        dx = 1e-4;

      double dp2 = characteristicCurve(ppq(2)+dx);
      derivative = (dp-dp2)/dx;
    }

    fDer(0) = -1.0;
    fDer(1) = 1.0;
    fDer(2) = derivative;
  }
  else // CLOSED, status is 0 or -1
  {
    // slightly unstable
    //out = ppq(2); 
    //fDer(2) = 1.0;

    // alternative solution
    out = ppq(1) - ppq(0) + (endHeight-startHeight) + 1e8 * ppq(2);
    fDer(0) = -1.0;
    fDer(1) =  1.0;
    fDer(2) =  1e8;
  }
  return out;
}

//--------------------------------------------------------------
void Pump::initialization(int mode, double value)
{
  if(mode == 0)
    volumeFlowRate = 1.;
  else
    volumeFlowRate = value;
}

//--------------------------------------------------------------
double Pump::characteristicCurve(double Q)
{
  double H = 0.;

  if(curveType == -1) // P0/Q i.e. POWER is constant
  {
    if(Q<2.8e-5)
      Q = 2.8e-5;
    H = constantPerformance / Q;
  }
  else if(curveType == 0) // a+bq^2
  {
    H = coeffCurve[0]*pow(revolutionNumber,2) + coeffCurve[1] * Q*Q;
  }
  else if(curveType == 1) // linear interpolation
  {
    H = linearInterpolation(qCurve, hCurve, Q);
  }
  else if(curveType == 2) // second order with least squares method
  {
    for(int i=0; i<coeffCurve.size(); i++)
      H += coeffCurve[i] * pow(revolutionNumber,2-i) * pow(Q,i);
  }
  else if(curveType == 3) // a+bq^c
  {
    if(Q<0.0) // avoiding complex numbers
      Q = -Q;

    H = coeffCurve[0] + coeffCurve[1] * pow(Q, coeffCurve[2]);
  }

  return H;
}

//--------------------------------------------------------------
void Pump::checkPump()
{
  if(pumpType >= 0)
  {
    if(volumeFlowRate < min(qCurve))
      cout << endl << "!WARNING! Volume flow rate (" << volumeFlowRate << ") is smaller than the minimum (" << min(qCurve) << ") at Pump: " << name << endl;
    if(volumeFlowRate > max(qCurve))
      cout << endl << "!WARNING! Volume flow rate (" << volumeFlowRate << ") is larger than the maximum (" << max(qCurve) << ") at Pump: " << name << endl;
  }
}

//--------------------------------------------------------------
double Pump::getDoubleProperty(string prop)
{
  double out = 0.0;
  if(prop == "referenceCrossSection" || prop == "reference_cross_section")
    out = referenceCrossSection;
  else if(prop == "head")
    out = characteristicCurve(volumeFlowRate);
  else if(prop == "revolutionNumber")
    out = revolutionNumber;
  else if(prop == "efficiencyNominal")
    out = efficiencyNominal;
  else if(prop == "startHeight")
    out = startHeight;
  else if(prop == "endHeight")
    out = endHeight;
  else if(prop == "massFlowRate")
    out = volumeFlowRate * density/1000.;
  else if(prop == "volumeFlowRate")
    out = volumeFlowRate;
  else if(prop == "velocity")
    out = volumeFlowRate / 1000. / referenceCrossSection;
  else if(prop == "headLoss" || prop == "headloss")
    out = abs(characteristicCurve(volumeFlowRate));
  else{
    cout << endl << endl << "DOUBLE Pump::getDoubleProperty() wrong argument:" << prop;
    cout << ", right values: massFlowRate | head | revolutionNumber | velocity | density | referenceCrossSection | user1 | user2" << endl << endl;
    out = 0.0;
  }
  return out;
}

//--------------------------------------------------------------
void Pump::setDoubleProperty(string prop, double value)
{
  if(prop == "volumeFlowRate")
    volumeFlowRate = value;
  else if(prop == "density")
    density = value;
  else if(prop == "revolutionNumber")
    revolutionNumber = value;
  else if(prop == "efficiencyNominal")
    efficiencyNominal = value;
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section")
    referenceCrossSection = value;
  else if(prop == "startHeight")
    startHeight = value;
  else if(prop == "endHeight")
    endHeight = value;
  else
  {  
    cout << endl << endl << "Pump::setDoubleProperty( DOUBLE ) wrong argument:" << prop;
    cout << ", right values: volumeFlowRate | density | revolutionNumber | referenceCrossSection | user1 | user2" << endl << endl;
  }
}

//--------------------------------------------------------------
void Pump::setIntProperty(string prop, int value)
{
  if(prop == "curveType")
    curveType = value;
  else
  {  
    cout << endl << endl << "void Pump::setStringProperty( INT ) wrong argument:" << prop;
    cout << ", right values: curveType " << endl << endl;
  }
}

//--------------------------------------------------------------
int Pump::getIntProperty(string prop)
{
  int out = 0;
  if(prop == "pumpType")
    out = pumpType;
  else if(prop == "curveType")
    out = curveType;
  else
  {
    cout << endl << endl << "INT Pipe::getIntProperty() wrong argument:" << prop;
    cout << ", right values: frictionModel" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
vector<double> Pump::getVectorProperty(string prop)
{
  vector<double> out;
  if(prop == "qCurve")
    out = qCurve;
  else if(prop == "hCurve")
    out = hCurve;
  else
  {  
    cout << endl << endl << "vector <double> Pump::getVectorProperty( string prop ) wrong argument:" << prop;
    cout << ", right values: qCurve | hCurve " << endl << endl;
  }
  return out;
}