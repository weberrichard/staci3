#include "Pipe.h"

using namespace std;

Pipe::Pipe(const string a_name, const string a_startNodeName, const string a_endNodeName, const double a_density, const double a_length, const double a_diameter, const double a_roughness, const double a_volumeFlowRate, bool a_isCheckValve, int a_frictionModel, double a_relativeViscosity) : Edge(a_name, a_diameter * a_diameter * M_PI / 4., a_volumeFlowRate, a_density)
{
  numberNode = 2;
  
  startNodeName = a_startNodeName;
  endNodeName = a_endNodeName;
  length = a_length;
  diameter = a_diameter;
  roughness = abs(a_roughness);
  relativeViscosity = a_relativeViscosity;

  frictionModel = a_frictionModel;
  setPipeConst();

  isCheckValve = a_isCheckValve;
  if(isCheckValve)
  {
    typeCode = 0;
    type = "PipeCV";
  }
  else
  {
    typeCode = 1;
    type = "Pipe";
  }

  status = 1; // i.e. open
}

//--------------------------------------------------------------
Pipe::~Pipe(){}

//--------------------------------------------------------------
string Pipe::info()
{
  ostringstream strstrm;
  strstrm << Edge::info();
  strstrm << "\n type                  : " << type;
  strstrm << "\n connection            : " << startNodeName << "(index:" << startNodeIndex << ") --> " << endNodeName << "(index:" << endNodeIndex << ")";
  strstrm << "\n length                : " << length << " [m]";
  strstrm << "\n diameter              : " << diameter << " [m]";
  strstrm << "\n roughness             : " << roughness << " [mm]";
  strstrm << "\n lambda                : " << lambda << " [-]";
  strstrm << "\n velocity              : " << getDoubleProperty("velocity") << " [m/s]";
  strstrm << "\n friction model        : ";
  if(frictionModel==0)
    strstrm << "H-W";
  else if(frictionModel==2)
    strstrm << "C-F";
  else
    strstrm << "!ERROR!";

  return strstrm.str();
}

//--------------------------------------------------------------
double Pipe::function(const VectorXd &ppq, VectorXd &fDer)// ppq = [Pstart, Pend, VolumeFlowRate]
{ 
   double out = 0.0;
   if(status == 1) // OPEN
   {
      fDer(0) = -1.0;
      fDer(1) = 1.0;
      if(frictionModel == 0) // H-W: HAZEN - WILLIAMS
      {
         double tmp = pipeConst * pow(abs(ppq(2)),0.85185);
         out = ppq(1) - ppq(0) + (endHeight-startHeight) +  ppq(2) * tmp;
         fDer(2) = 1.85185 * tmp;
      }
      else if(frictionModel == 1) // D-W: DARCY - WEISBACH
      {
         double lambda = getDarcy(ppq(2));
         out = ppq(1) - ppq(0) + (endHeight-startHeight) + lambda * pipeConst * ppq(2) * abs(ppq(2));
         fDer(2) = 2 * lambda * pipeConst * abs(ppq(2)); // d(lambda) / d(q) is neglected
      }
      else if(frictionModel == 2) // C-F: CONSTANT FRICTION COEFFICIENT
      {
         out = ppq(1) - ppq(0) + (endHeight-startHeight) + pipeConst * ppq(2) * abs(ppq(2));
         fDer(2) = 2 * pipeConst * abs(ppq(2));
      }

      if(fDer(2) < 2*pipeConst*2.832e-5)
         fDer(2) = 2*pipeConst*2.832e-5;
   }
   else // CLOSED, status is 0 or -1
   {
      out = ppq(1) - ppq(0) + (endHeight-startHeight) + 1e8 * ppq(2);
      fDer(0) = -1.0;
      fDer(1) =  1.0;
      fDer(2) =  1e8;
   }

  return out;
}

//--------------------------------------------------------------
void Pipe::setPipeConst()
{
   if(frictionModel == 0) // HAZEN
      pipeConst = 10.654 * pow(roughness,-1.85185) * pow(diameter,-4.87) * length;
   else if(frictionModel == 1) // DARCY
      pipeConst = length *8. / (diameter*diameter*diameter*diameter*diameter * gravity * M_PI*M_PI);
   else if(frictionModel == 2) // CONSTANT FRICTION COEFFICIENT
      pipeConst = roughness*length *8. / (diameter*diameter*diameter*diameter*diameter * gravity * M_PI*M_PI);
}

//--------------------------------------------------------------
double Pipe::getDarcy(double q)
{  
   double nu = relativeViscosity * 1e-6;
   double v = q / (diameter*diameter*M_PI/4.);
   double Re = v * diameter / nu;

   double out;
   if(Re < 2000.) // laminar: Hagen - Poiseuille
   {
      out = 64./Re;
   }
   else if(Re > 4000.) // turbulent: Colebrook - White approximation
   {  
      double den = log(roughness/(3.7*diameter) + 5.74/pow(Re,0.9))/log(10.);
      out = 0.25 / pow(den,2);
   }
   else // transitional: cubic interpolation like epanet2
   {
      // according epanet
      double y3 = -0.86859 * log(roughness/(3.7*diameter) + 5.74 / pow(4000.,0.9));
      double y2 = roughness / (3.7*diameter) + 5.74 / pow(Re,0.9);
      double fa = pow(y3,-2.);
      double fb = fa * (2-0.00514215/(y2*y3));
      double R = Re/2000.;
      double x1 = 7*fa-fb;
      double x2 = 0.128 - 17*fa + 2.5*fb;
      double x3 = -0.128 + 13*fa - 2*fb;
      double x4 = R*(0.032 - 3*fa + 0.5*fb);
      out = x1 + R*(x2 + R*(x3 + x4));
   }

   return out;
}

//--------------------------------------------------------------
void Pipe::initialization(int mode, double value)
{
   if(status == 1) // open
   {
      if(mode == 0)
         volumeFlowRate = 1.;
      else
         volumeFlowRate = value;
   }
   else // closed
   {
      volumeFlowRate = 0.0;
   }
}

//--------------------------------------------------------------
double Pipe::functionParameterDerivative(int parameter)
// 0: roughness coefficient, 1: diameter
{
  double out = 0.0;
  if(status == 1) // OPEN
  {
    if(parameter == 0) // roghness
    {
      if(frictionModel == 0) // H-W: HAZEN - WILLIAMS
      {
        out = pipeConst * (-1.85185) / roughness * volumeFlowRate * pow(abs(volumeFlowRate),0.85185);
      }
      else if(frictionModel == 2)// C-F: CONSTANT FRICTION COEFFICIENT
      {
        out = pipeConst / roughness * volumeFlowRate * abs(volumeFlowRate);
      }
    }
    else if(parameter == 1) // diameter
    {
      if(frictionModel == 0) // H-W: HAZEN - WILLIAMS
      {
        out = pipeConst * (-4.87) / diameter * volumeFlowRate * pow(abs(volumeFlowRate),0.85185);
      }
      else if(frictionModel == 2)// C-F: CONSTANT FRICTION COEFFICIENT
      {
        out = pipeConst * (-5) / diameter * volumeFlowRate * abs(volumeFlowRate);
      }
    }
    else
    {
      cout << endl << "!WARNING! Pipe::functionParameterDerivative(parameter), unkown input: parameter = " << parameter << endl << "Available parameters: 0: roughness, 1: diameter" << endl;
      cout << endl << "Name of pipe: " << name << endl;
    }
  }
  else // CLOSED
  {
    out = 0.0;
  }
  return out;
}

//--------------------------------------------------------------
double Pipe::getDoubleProperty(string prop)
{
  double out = 0.;
  if(prop == "roughness")
    out = roughness;
  else if(prop == "diameter")
    out = diameter;
  else if(prop == "length")
    out = length;
  else if(prop == "lambda")
    out = lambda;
  else if(prop == "pipeConst")
    out = pipeConst;
  else if(prop == "volume")
    out = referenceCrossSection*length;
  else if(prop == "segment")
    out = (double)segment;
  else if(prop == "status")
    out = (double)status;
  else if(prop == "userOutput")
    out = userOutput;
  else if(prop == "Rh")
    out = diameter / 2.;
  else if(prop == "massFlowRate" || prop == "mass_flow_rate")
    out = volumeFlowRate * density;
  else if(prop == "volumeFlowRate" || prop == "volume_flow_rate")
    out = volumeFlowRate;
  else if(prop == "volumeFlowRateLPS" || prop == "volume_flow_rate")
    out = volumeFlowRate*1000.;
  else if(prop == "volumeFlowRateAbs")
    out = abs(volumeFlowRate);
  else if(prop == "volumeFlowRateAbsLPS")
    out = abs(volumeFlowRate)*1000.;
  else if(prop == "velocity")
    out = volumeFlowRate / referenceCrossSection;
  else if(prop == "density")
    out = density;
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section" || prop == "Aref")
    out = referenceCrossSection;
  else if(prop == "startHeight")
    out = startHeight;
  else if(prop == "endHeight")
    out = endHeight;

  else
  {
    cout << endl << endl << "DOUBLE Pipe::getDoubleProperty() wrong argument:" << prop;
    cout << ", right values: diamater | roughness | length | massFlowRate | velocity | density | referenceCrossSection | user1 | user2 | startHeight | endHeight | volume" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
int Pipe::getIntProperty(string prop)
{
  int out = 0;
  if(prop == "frictionModel" || prop == "friction_model")
    out = frictionModel;
  else
  {
    cout << endl << endl << "INT Pipe::getIntProperty() wrong argument:" << prop;
    cout << ", right values: frictionModel" << endl << endl;
  }
  return out;
}

//--------------------------------------------------------------
void Pipe::setDoubleProperty(string prop, double value)
{
  if(prop == "roughness")
  {
    roughness = value;
    setPipeConst();
  }
  else if(prop == "diameter")
  {
    diameter = value;
    setPipeConst();
  }
  else if(prop == "length")
  {
    length = value;
    setPipeConst();
  }
  else if(prop == "volumeFlowRate")
    volumeFlowRate = value;
  else if(prop == "density")
    density = value;
  else if(prop == "referenceCrossSection" || prop == "reference_cross_section")
    referenceCrossSection = value;
  else if(prop == "startHeight")
    startHeight = value;
  else if(prop == "endHeight")
    endHeight = value;
  else if(prop == "userOutput")
    userOutput = value;
  else
  {  
    cout << endl << endl << "Pipe::setDoubleProperty( DOUBLE ) wrong argument:" << prop;
    cout << ", right values: diameter | roughness | length | volumeFlowRate | velocity | density | referenceCrossSection | userOutput | startHeight | endHeight" << endl << endl;
  }
}

//--------------------------------------------------------------
void Pipe::setIntProperty(string prop, int value)
{
  if(prop == "frictionModel" || prop == "friction_model")
    frictionModel = value;
  else
  {  
    cout << endl << endl << "Pipe::setIntProperty( DOUBLE ) wrong argument:" << prop;
    cout << ", right values: startNodeIndex | endNodeIndex | numberNode" << endl << endl;
  }
}
