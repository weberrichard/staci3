#include "Node.h"

Node::Node(const string a_name, const double a_xPosition, const double a_yPosition, const double a_geodeticHeight, const double a_demand, const double a_head, const double a_density)
{
  density = a_density;
  demand = a_demand;
  geodeticHeight = a_geodeticHeight;
  xPosition = a_xPosition;
  yPosition = a_yPosition;
  name = a_name;
  head = a_head;
}

//--------------------------------------------------------------
Node::~Node(){}

//--------------------------------------------------------------
void Node::initialization(int mode, double value)
{
  if(mode == 0)
    head = 50.; 
  else
    head = value - geodeticHeight;
} 

//--------------------------------------------------------------
double Node::function(const VectorXd &pq, bool isPressureDemand, VectorXd &fDer)
// pq = [p, Qin1, Qin2, ..., Qout1, Qout2, ...]
// fDer = [df/dp, df/dQin1, df/dQin2, ..., df/dQout1, df/dQou2, ...]
{
  double out = 0.0;
  if(status) // if the node is active
  {
    if(isPressureDemand) // if the demands are depending the nodal pressure
    {
      if(pq(0)>=pdDesiredPressure)
      {
        out -= demand;
        consumption = demand;
        consumptionPercent = 100.0;
        fDer(0) = 0.0;
      }
      else if(pq(0)<pdDesiredPressure && pq(0)>pdMinPressure)
      {
        consumption = getConsumption(pq(0));
        consumptionPercent = consumption / demand;
        out -= consumption;
        fDer(0) = -consumption/(pdExponent*(pq(0)-pdMinPressure));
      }
      else
      {
        out -= 0.0;
        consumption = 0.0;
        consumptionPercent = 0.0;
        fDer(0) = 0.0;
      }
    }
    else // for constant demands
    {
      out -= demand;
      consumption = demand;
      consumptionPercent = 100.0;
      fDer(0) = 0.0;
    }

    for(int i=0; i<edgeIn.size(); i++)
    {
      out += pq(1 + i);
      fDer(1 + i) = 1.0;
    }

    for(int i=0; i<edgeOut.size(); i++)
    {
      out -= pq(1 + edgeIn.size() + i);
      fDer(1 + edgeIn.size() + i) = -1.0;
    }
  }
  else // if the node is NOT active
  {
    out = pq(0);
    consumption = 0.0;
    consumptionPercent = 0.0;
    fDer(0) = 1.0;
  }
  return out;
}
//--------------Overloaded, the evade the false forest in every other cases-----------------------//
double Node::function(const VectorXd &pq, bool isPressureDemand, bool isLeakage, VectorXd &fDer)
// pq = [p, Qin1, Qin2, ..., Qout1, Qout2, ...]
// fDer = [df/dp, df/dQin1, df/dQin2, ..., df/dQout1, df/dQou2, ...]
{
  double out = 0.0;
  if(status) // if the node is active
  {
    if(isLeakage) // if the demands are depending the nodal pressure
    {
      if(pq(0)>leakageMinPressure)
      {
        consumption = demand;
        leakageAmount = leakageConstant*pow((pq(0)-leakageMinPressure),leakageExponent);
        out -= consumption;
        out -= leakageAmount;
        fDer(0) = -leakageConstant*leakageExponent*pow((pq(0)-leakageMinPressure),(leakageExponent-1));
      }
      else
      {
        out -= 0.0;
        consumption = 0.0;
        consumptionPercent = 0.0;
        fDer(0) = 0.0;
      }
    }
    if(isPressureDemand) // if the demands are depending the nodal pressure
    {
      if(pq(0)>=pdDesiredPressure)
      {
        out -= demand;
        consumption = demand;
        consumptionPercent = 100.0;
        fDer(0) = 0.0;
      }
      else if(pq(0)<pdDesiredPressure && pq(0)>pdMinPressure)
      {
        consumption = getConsumption(pq(0));
        consumptionPercent = consumption / demand;
        out -= consumption;
        fDer(0) = -consumption/(pdExponent*(pq(0)-pdMinPressure));
      }
      else
      {
        out -= 0.0;
        consumption = 0.0;
        consumptionPercent = 0.0;
        fDer(0) = 0.0;
      }
    }
    else // for constant demands
    {
      out -= demand;
      consumption = demand;
      consumptionPercent = 100.0;
      fDer(0) = 0.0;
    }

    for(int i=0; i<edgeIn.size(); i++)
    {
      out += pq(1 + i);
      fDer(1 + i) = 1.0;
    }

    for(int i=0; i<edgeOut.size(); i++)
    {
      out -= pq(1 + edgeIn.size() + i);
      fDer(1 + edgeIn.size() + i) = -1.0;
    }
  }
  else // if the node is NOT active
  {
    out = pq(0);
    consumption = 0.0;
    consumptionPercent = 0.0;
    fDer(0) = 1.0;
  }
  return out;
}
//--------------------------------------------------------------
double Node::functionParameterDerivative(bool isPressureDemand)
{
  double out;
  if(status) // if the node is active
  {
    if(isPressureDemand) // if the demands are depending the nodal pressure
    {
      if(head<pdMinPressure)
      {
        out = 0.0;
      }
      else if(head<pdDesiredPressure && head>pdMinPressure)
      {
        out = getConsumption(head) / demand;
      }
      else
      {
        out = -1.0;
      }
    }
    else // for constant demands
    {
      out = -1.0;
    }
  }
  else // if the node is NOT active
  {
    out = 0.0;
  }
  return out;
}

//--------------------------------------------------------------
double Node::getConsumption(double head)
{
  return demand*pow((head-pdMinPressure)/(pdDesiredPressure-pdMinPressure),1./pdExponent);
}

//--------------------------------------------------------------
void Node::setProperty(string prop, double value)
{
  if(prop == "demand")
    demand = value;
  else if(prop == "consumption")
  {
    consumption = value;
    if(demand == 0.0)
      consumptionPercent = 0.0;
    else
      consumptionPercent = 100.*consumption/demand;
  }
  else if(prop == "consumptionPercent")
  {
    consumptionPercent = value;
    consumption = consumptionPercent/100.*demand;
  }
  else if(prop == "userOutput")
    userOutput = value;
  else if(prop == "head")
    head = value;
  else if(prop == "density")
    density = value;
  else if(prop == "height" || prop == "geodeticHeight")
    geodeticHeight = value;
  else if(prop == "xPosition")
    xPosition = value;
  else if(prop == "yPosition")
    yPosition = value;
  else if(prop == "leakageExponent")
    leakageExponent = value;
  else if(prop == "leakageConstant")
    leakageConstant = value;
  else if(prop == "leakageMinPressure")
    leakageMinPressure = value;
  else
  {
    cout << endl << endl << "Node::getProperty() wrong argument:" << prop;
    cout << ", right values: demand|head|pressure|density|height|xPosition|yPosition" << endl << endl;
  }
}

//--------------------------------------------------------------
double Node::getProperty(string prop)
{
  double out = 0.0;

  if(prop == "demand")
    out = demand;
  else if(prop == "consumption")
  {
    if(status == 1)
      out = consumption;
    else
      out = 0.;
  }
  else if(prop == "leakage")
  {
    out = leakageAmount;
  }
  else if(prop == "consumptionPercent")
  {
    if(status == 1)
      out = consumptionPercent;
    else
      out = 0.;
  }
  else if(prop == "userOutput")
    out = userOutput;
  else if(prop == "pressure")
    out = head * density * 9.81;
  else if(prop == "head")
    out = head;
  else if(prop == "density")
    out = density;
  else if(prop == "height" || prop == "geodeticHeight")
    out = geodeticHeight;
  else if(prop == "xPosition")
    out = xPosition;
  else if(prop == "yPosition")
    out = yPosition;
  else if(prop == "segment")
    out = (double)segment;
  else if(prop == "status")
    out = (double)status;
  else if(prop == "leakageExponent")
    out = leakageExponent;
  else if(prop == "leakageConstant")
    out = leakageConstant;
  else if(prop == "leakageMinPressure")
    out = leakageMinPressure;
  else if(prop == "waterAge")
    out = waterAge.back();
  else if(prop == "chlorine")
    out = chlorineConcentration.back();
  else if(prop == "microbesWater" || prop == "Cf")
    out = microbesWater.back();
  else if(prop == "microbesWall" || prop == "Cb")
    out = microbesWall.back();
  else if(prop == "substratWater" || prop == "Sf")
    out = substratWater.back();
  else if(prop == "substratWall" || prop == "Sb")
    out = substratWall.back();
  else if(prop == "timeSteps" || prop == "TS")
    out = TS.back();
  else
  {
    cout << endl << endl << "Node::getProperty() wrong argument:" << prop;
    cout << ", right values: demand|head|pressure|density|height|xPosition|yPosition|userOutput" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
string Node::info(bool check_if_lonely)
{
  ostringstream strstrm;
  strstrm << "\n Node name       : " << name;
  strstrm << "\n open/closed     : ";
  strstrm << "\n height          : " << geodeticHeight << " m";
  strstrm << "\n head            : " << head << " m";
  strstrm << "\n pressure        : " << head*density * 9.81 << " Pa";
  strstrm << "\n desnity         : " << density << " kg/m3";
  strstrm << "\n demand          : " << demand << " l/s = " << demand * 3.6 << " m3/h";
  strstrm << "\n consumption     : " << consumption << " l/s = " << consumption * 3.6 << " m3/h";
  strstrm << "\n segment         : " << segment;
  strstrm << "\n incoming edges  : ";
  for(vector<int>::iterator it = edgeIn.begin(); it != edgeIn.end(); it++)
      strstrm << *it << " ";
  strstrm << "\n outgoing edges  : ";
  for(vector<int>::iterator it = edgeOut.begin(); it != edgeOut.end(); it++)
      strstrm << *it << " ";
  strstrm << endl;

  if (check_if_lonely && ((edgeIn.size() + edgeOut.size()) == 0))
  {
    strstrm << "\n!!! ERROR !!! Lonely node: " << name << " !!!\n";
    cout << strstrm.str();
    exit(-1);
  }

  return strstrm.str();
}

void Node::appendTimeSeries(string prop, double value)
{
  if (prop == "waterAge")
  {
    waterAge.push_back(value);
  }
  else if (prop == "chlorine")
  {
    chlorineConcentration.push_back(value);
  }
  else if (prop == "microbesWater" || prop == "Cf") 
  {
    microbesWater.push_back(value);
  }
    else if (prop == "microbesWall" || prop == "Cb")
  {
    microbesWall.push_back(value);
  }
    else if (prop == "substratWater" || prop == "Sf")
  {
    substratWater.push_back(value);
  }
    else if (prop == "substratWall" || prop == "Sb")
  {
    substratWall.push_back(value);
  }
    else if (prop == "timeSteps" || prop == "TS")
  {
    TS.push_back(value);
  }
}

vector<double> Node::getTimeSeries(string prop)
{
  if (prop == "waterAge")
  {
    return waterAge;
  }
  else  if (prop == "chlorine")
  {
    return chlorineConcentration;
  }
    else if (prop == "microbesWater" || prop == "Cf")
  {
    return microbesWater;
  }
    else if (prop == "microbesWall" || prop == "Cb")
  {
    return microbesWall;
  }
    else if (prop == "substratWater" || prop == "Sf")
  {
    return substratWater;
  }
    else if (prop == "substratWall" || prop == "Sb")
  {
    return substratWall;
  }
    else if (prop == "timeSteps" || prop == "TS")
  {
    return TS;
  }
}