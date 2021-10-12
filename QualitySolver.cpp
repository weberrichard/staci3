#include "QualitySolver.h"

//--------------------------------------------------------------------------//
QualitySolver::QualitySolver(string spr_filename) : SeriesHydraulics(spr_filename){}//
QualitySolver::~QualitySolver(){}                                           // 
//--------------------------------------------------------------------------//

int QualitySolver::Embedded_Prince_Dormand_solver(double t, double h, double tmax, double tolerance, double DX, string WhichSolver)
{
  int i;
  int last_interval = 0, edgeCounter = 0;
  int counter = 0., timestep = 0, attemptsTried = 0;

  bool repeat = true;

  vector<vector<double>> temp_x(2, vector<double>(SourceDimension)); 
  vector<double> x_previous(SourceDimension);
  double length = 0., maxError = 0.;
  double scale;
  vector<double> err;
  double xx, xx_min = 0.;

  double h_min_old = 0., h_min_new = 0.;
  vector<int> numberOfIntervals;
  vector < vector< vector<double> > > x;
  vector< vector<double> > h_actual;

  //--- Check for missdefined parameters ---//
  if (tmax < t || h <= 0.0) return -2;  
  if (tmax == t) return 0; 
  //----------------------------------------//

  //--------- Discretize the network (standard grid) ---------//
  numberOfIntervals = generateNumberOfIntervals(DX);
  //----------------------------------------------------------//

  //----------------- Initialization phase--------------------// 
  x = initializeWithStartingCondition(numberOfIntervals);
  
  //updateNodeVariables(x, WhichSolver, 0.);

  //----------------------------------------------------------//
  h_actual = updateInitializeTimestep(h_actual, h, numberOfIntervals);
  while ( t <= tmax ) 
  {
    repeat = true;
    while(repeat != false)
    { 
      edgeCounter = 0;
      h_min_old = getSmallestTimestep(h_actual);
      for (int i = 0; i < edges.size(); i++)
      {
        if (edges.at(i)->getEdgeStringProperty("type") == "Pipe")
        {
          for (int j = 0; j < SourceDimension; ++j)
          {
            if (edges.at(i)->getDoubleProperty("velocity") < 0) 
            {
              x[edgeCounter][0][j] = nodes.at(edges.at(i)->getEdgeIntProperty("endNodeIndex"))->getProperty(ParameterList[j]);
            }
            else
            {
              x[edgeCounter][0][j] = nodes.at(edges.at(i)->getEdgeIntProperty("startNodeIndex"))->getProperty(ParameterList[j]);
            }
          }
          counter = 0;
          h = h_min_old;
          h = min(h, tmax - t); 
          length = 0.;   
          length_max = edges.at(i)->getDoubleProperty("length"); 
          for (int stepCounter = 1; stepCounter < x[edgeCounter].size(); ++stepCounter)
          {
              scale = 1.0;
              length += DX;
              for (int j = 0; j < SourceDimension; ++j)
              {
                x_previous[j] = x[edgeCounter][stepCounter-1][j];
                temp_x[0][j] = x[edgeCounter][stepCounter][j];
              }
              for (int j = 0; j < attempts; j++)   
              { 
                maxError = 0.;
                xx = 0;
                xx_min = 0;
                err = Runge_Kutta(temp_x, t, h, x_previous, abs(edges.at(i)->getDoubleProperty("velocity")), DX, WhichSolver);
                for (int k = 0; k < SourceDimension; ++k) 
                { 
                  if (abs(err[k])>maxError)   
                  {
                    maxError = abs(err[k]);
                  }
                  xx = (temp_x[0][0] == 0.0) ? tolerance : fabs(temp_x[0][k]);
                  if(xx > xx_min)
                  { 
                    xx_min = xx;
                  }
                } 
                if (maxError == 0.0)
                {  
                    scale = maxScaleFactor; break;
                }
                scale = 0.8 * pow(( tolerance * xx_min /  maxError ),0.25);
                scale = min( max(scale,minScaleFactor), maxScaleFactor); 
                if ( maxError < ( tolerance * xx_min ) ) break; 
                h *= scale;
                if ( t + h > tmax )
                {
                  h = tmax - t; 
                }
                else if ( t + h + 0.5 * h > tmax )
                {
                  h = 0.5 * h;
                }
                attemptsTried = j;
              }
              if (attemptsTried >= attempts ) 
              { 
                h *= scale; return -1; 
              }
              for (int j = 0; j < SourceDimension; ++j)   
              {
                if(repeat == true) 
                {
                  x[edgeCounter][stepCounter][j] = temp_x[1][j];
                }
              }
              h_actual[edgeCounter][stepCounter] = h;
          } 
          edgeCounter += 1;
        }
      }
      h_min_new = getSmallestTimestep(h_actual);
      if (h_min_new == h_min_old)  
      { 
        repeat = false; 
        h = h_min_new; 
      }
    } 
    updateNodeVariables(x, WhichSolver, h);
    t += h;
    cout << "  Elapsed time: " << t << " velocity: " << edges.at(1)->getEdgeDoubleProperty("velocity") <<  endl;
    h *= scale;
    if ( last_interval ) break;
    if (  t + h > tmax ) { last_interval = 1; h = tmax - t; }
    else if ( t + h + 0.5 * h > tmax ) h = 0.5 * h;
  }
  return 0;
}

vector<double> QualitySolver::Runge_Kutta(vector< vector<double> >& x, double t0, double h, vector<double> x_previous, double flow, double DX, string WhichSolver)
{
    const double two_thirds = 2.0 / 3.0;
    const double one_seventwoninths = 1.0 / 729.0;
    const double one_twoninesevenzero = 1.0 / 2970.0;
    const double one_twofivetwozero = 1.0 / 2520.0;
    const double one_fiveninefourzero = 1.0 / 5940.0;
 
    double h5 = 0.2 * h;

    vector<double> k1, k2, k3, k4, k5, k6, out(SourceDimension);

    vector< vector<double> > RK_Inputs(6, vector<double>());

    for (int i = 0; i < SourceDimension; ++i)
    {
      RK_Inputs[0].push_back(x[0][i]);
    }

    k1 = source(t0, RK_Inputs[0], x_previous, flow, DX, WhichSolver);

    for (int i = 0; i < SourceDimension; ++i)
    {
      RK_Inputs[1].push_back(x[0][i] + h5 * k1[i]);
    }

    k2 = source(t0+h5, RK_Inputs[1], x_previous, flow, DX, WhichSolver);

    for (int i = 0; i < SourceDimension; ++i)
    {
      RK_Inputs[2].push_back(x[0][i] + h * ( 0.075 * k1[i] + 0.225 * k2[i]));
    }

    k3 = source(t0+0.3*h, RK_Inputs[2], x_previous, flow, DX, WhichSolver);

    for (int i = 0; i < SourceDimension; ++i)
    {
      RK_Inputs[3].push_back(x[0][i] + h * ( 0.3 * k1[i] - 0.9 * k2[i] + 1.2 * k3[i]));
    }

    k4 = source(t0+0.6*h, RK_Inputs[3], x_previous, flow, DX, WhichSolver);

    for (int i = 0; i < SourceDimension; ++i)
    { 
      RK_Inputs[4].push_back(x[0][i] + one_seventwoninths * h * ( 226.0 * k1[i] - 675.0 * k2[i] + 880.0 * k3[i] + 55.0 * k4[i]));
    }

    k5 = source(t0+two_thirds * h,  RK_Inputs[4], x_previous, flow, DX, WhichSolver);

    for (int i = 0; i < SourceDimension; ++i)
    {
      RK_Inputs[5].push_back(x[0][i] + one_twoninesevenzero * h * ( - 1991 * k1[i] + 7425.0 * k2[i] - 2660.0 * k3[i] - 10010.0 * k4[i] + 10206.0 * k5[i]));
    }

    k6 = source(t0+h, RK_Inputs[5], x_previous, flow, DX, WhichSolver);

    for (int i = 0; i < SourceDimension; ++i)
    {
      x[1][i] = x[0][i] +  one_fiveninefourzero * h * ( 341.0 * k1[i] + 3800.0 * k3[i] - 7975.0 * k4[i] + 9477.0 * k5[i] + 297.0 * k6[i]);
      out[i] = one_twofivetwozero * (77.0 * k1[i] - 400.0 * k3[i] + 1925.0 * k4[i] - 1701.0 * k5[i] + 99.0 * k6[i]);
    }
    return out;
}

vector<double> QualitySolver::source(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX, string WhichSolver)
{
  vector<double> out;
  if (WhichSolver == "waterAge")
  {
    out = wa->sourceTermWaterAge(t, x_actual, x_upwind, flow, DX);
  }
  else if (WhichSolver == "chlorine")
  {
    out = ch->sourceTermChlorine(t, x_actual, x_upwind, flow, DX);
  }
  else if (WhichSolver == "biofilm")  
  {
    out = bf->sourceTermBiofilm(t, x_actual, x_upwind, flow, DX);
  }
  return out;
}

vector<int> QualitySolver::generateNumberOfIntervals(double DX)
{ 
  vector<int> numberOfIntervalForAllPipe;
  for (int i = 0; i < edges.size(); ++i)
  {  
    if (edges.at(i)->getEdgeStringProperty("type") == "Pipe")
    {
      numberOfIntervalForAllPipe.push_back(round(edges.at(i)->getDoubleProperty("length")/DX)+1);
    }
  }
  return numberOfIntervalForAllPipe; 
}

vector < vector< vector<double> > > QualitySolver::initializeWithStartingCondition(vector<int> numberOfIntervals)
{
  vector<double> input(SourceDimension);
  vector< vector<double> > input2D;
  vector< vector< vector<double> > > out(numberOfIntervals.size(), vector< vector<double> >());
  for (int i = 0; i < numberOfIntervals.size(); ++i)
  {
    for (int j = 0; j < numberOfIntervals[i]; ++j)
    {
      for (int k = 0; k < SourceDimension; ++k)
      {
        input[k] = startingCondition[k];
      }
      input2D.push_back(input);
    }
    out[i] = input2D;
    input2D.erase(input2D.begin(),input2D.end());
  }
  for (int i = 0; i < nodes.size(); ++i)
  {
    nodes.at(i)->appendTimeSeries("timeSteps", 0);
    for (int j = 0; j < SourceDimension; ++j)
    {
      nodes.at(i)->appendTimeSeries(ParameterList[j],startingCondition[j]);
    }
  }
  return out;
}

vector< vector<double> > QualitySolver::updateInitializeTimestep(vector< vector<double> > out, double initTimestep, vector<int> numberOfIntervals)
{
  if (out.size() == 0)
  {
    for (int i = 0; i < numberOfIntervals.size(); ++i)
    {
      out.push_back(vector<double> ());
      for (int j = 0; j < numberOfIntervals[i]; ++j)
      {
        out[i].push_back(initTimestep);
      }
    }
  }
  else
  {
    for (int i = 0; i < out.size(); ++i)
    {
      for (int j = 0; j < out[i].size(); ++j)
      {
        out[i][j] = initTimestep;
      }
    }  
  }
  return out;
}

vector< vector<double> > QualitySolver::updateInitializeTimestep(vector< vector<double> > out, double initTimestep)
{
  for (int i = 0; i < out.size(); ++i) 
  {
    for (int j = 0; j < out[i].size(); ++j)
    {
      out[i][j] = initTimestep;
    }
  }  
  return out;
}

double QualitySolver::getSmallestTimestep(vector< vector<double> > input)
{
  double out = input[0][0]; 
  for (int i = 0; i < input.size(); ++i)
  {
    for (int j = 0; j < input[i].size(); ++j)
    {
      if (input[i][j] <= out)
      {
        out = input[i][j];
      }
    }
  }
  return out;  
}

void QualitySolver::updateNodeVariables(vector< vector< vector<double> > > x, string WhichSolver, double h)
{
  int counter = 0, counterOffset = 0, elementIndex = 0;
  double poolFlow, poolVolumeActual;
  bool pump = false, pool = false, press = false, valve = false, valveISO = false, begin = false, stopCount = false;
  vector<double> VolFlowRates, poolActual(SourceDimension), poolUpstreamNode(SourceDimension), temporal_value_2, nodalValue(SourceDimension);
  vector< vector<double> > temporal_value, NodeInputs;
  //---------------------------------------- Count the offset, for the calculated elements ------------------------------------//
  for (int i = 0; i < edges.size(); ++i)
  {
    if (edges.at(i)->getEdgeStringProperty("type") != "Pipe" && stopCount == false)
    {
      counterOffset += 1;
    }
    if (edges.at(i)->getEdgeStringProperty("type") == "Pipe")
    {
      stopCount = true;
    }
  }
  //----------------------------------------------------Identify the node's role-----------------------------------------------//
  for (int i = 0; i < nodes.size(); ++i)
  {
    for (int j = 0; j < SourceDimension; ++j)
    {
      nodalValue[j] = 0;
    }
    for (int j = 0; j < poolIndex.size(); ++j)
    {
      if (i == edges.at(poolIndex[j])->getEdgeIntProperty("startNodeIndex"))
      {
        pool = true;
        begin = true;
        elementIndex = j;
      }
      else if(i == edges.at(poolIndex[j])->getEdgeIntProperty("endNodeIndex"))
      {
        pool = true;
        elementIndex = j;
      }
    }
    for (int j = 0; j < presIndex.size(); ++j)
    {
      if (i == edges.at(presIndex[j])->getEdgeIntProperty("startNodeIndex"))
      {
        press = true;
        begin = true;
        elementIndex = j;
      }
      else if(i == edges.at(presIndex[j])->getEdgeIntProperty("endNodeIndex"))
      {
        press = true;
        elementIndex = j;
      }
    }
    for (int j = 0; j < pumpIndex.size(); ++j)
    {
      if (i == edges.at(pumpIndex[j])->getEdgeIntProperty("startNodeIndex"))
      {
        pump = true;
        begin = true;
        elementIndex = j;
      }
      else if(i == edges.at(pumpIndex[j])->getEdgeIntProperty("endNodeIndex"))
      {
        pump = true;
        elementIndex = j;
      }
    }
    for (int j = 0; j < valveIndex.size(); ++j)
    {
      if (i == edges.at(valveIndex[j])->getEdgeIntProperty("startNodeIndex"))
      {
        valve = true;
        begin = true;
        elementIndex = j;
      }
      else if(i == edges.at(valveIndex[j])->getEdgeIntProperty("endNodeIndex"))
      {
        valve = true;
        elementIndex = j;
      }
    }
    for (int j = 0; j < valveISOIndex.size(); ++j)
    {
      if (i == edges.at(valveISOIndex[j])->getEdgeIntProperty("startNodeIndex"))
      {
        begin = true;
        valveISO = true;
        elementIndex = j;
      }
      else if(i == edges.at(valveISOIndex[j])->getEdgeIntProperty("endNodeIndex"))
      {
        valveISO = true;
        elementIndex = j;
      }
    }
//---------------------------------------------------------------------------------------------------------//
    if (pool == false && press == false)
    {
      counter = 0;
      for (int j = 0; j < nodes.at(i)->edgeIn.size(); ++j)
      {
        if (edges.at(nodes.at(i)->edgeIn[j])->getEdgeStringProperty("type") == "Pipe" && (edges.at(nodes.at(i)->edgeIn[j])->getEdgeDoubleProperty("velocity") > 0))
        {
          NodeInputs.push_back(vector<double>());
          temporal_value = x[nodes.at(i)->edgeIn[j]-counterOffset];
          temporal_value_2 = temporal_value.back(); 
          for (int k = 0; k < SourceDimension; ++k)
          {
            NodeInputs[counter].push_back(temporal_value_2[k]);
            VolFlowRates.push_back(edges.at(nodes.at(i)->edgeIn[j])->getDoubleProperty("volumeFlowRate"));
            counter += 1;
          }
        }
      }
      for (int j = 0; j < nodes.at(i)->edgeOut.size(); ++j)
      {
        if (edges.at(nodes.at(i)->edgeOut[j])->getEdgeStringProperty("type") == "Pipe" && (edges.at(nodes.at(i)->edgeOut[j])->getEdgeDoubleProperty("velocity") < 0))
        {
          NodeInputs.push_back(vector<double>());
          temporal_value = x[nodes.at(i)->edgeIn[j]-counterOffset];
          temporal_value_2 = temporal_value.back();
          for (int k = 0; k < SourceDimension; ++k)
          {
            NodeInputs[counter].push_back(temporal_value_2[k]);
            VolFlowRates.push_back(edges.at(nodes.at(i)->edgeIn[j])->getDoubleProperty("volumeFlowRate"));
            counter += 1;
          }
        }
      }
      if(NodeInputs.size() == 0)
      {
        if (pump == true && begin == true)
        {
          pump = false;
          for (int j = 0; j < SourceDimension; ++j)
          {
            nodalValue[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("endNodeIndex"))->getProperty(ParameterList[j]);
          }
        }
        else if (pump == true && begin == false)
        {
          pump = false;
          for (int j = 0; j < SourceDimension; ++j)
          {
            nodalValue[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("startNodeIndex"))->getProperty(ParameterList[j]);
          }
        }
        else if (valve == true && begin == true)
        {
          valve = false;
          for (int j = 0; j < SourceDimension; ++j)
          {
            nodalValue[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("endNodeIndex"))->getProperty(ParameterList[j]);
          }
        }
        else if (valve == true && begin == false)
        {
          valve = false;
          for (int j = 0; j < SourceDimension; ++j)
          {
            nodalValue[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("startNodeIndex"))->getProperty(ParameterList[j]);
          }
        }
        else if (valveISO == true && begin == true)
        {
          valveISO = false;
          for (int j = 0; j < SourceDimension; ++j)
          {
            nodalValue[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("endNodeIndex"))->getProperty(ParameterList[j]);
          }
        }
        else if (valveISO == true && begin == false)
        {
          valveISO = false;
          for (int j = 0; j < SourceDimension; ++j)
          {
            nodalValue[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("startNodeIndex"))->getProperty(ParameterList[j]);
          }
        }
      }
      else
      {
        if (WhichSolver == "waterAge")
        {
          nodalValue = wa->nodeEquationWaterAge(NodeInputs, VolFlowRates);
        }
        else if (WhichSolver == "chlorine")
        {
          nodalValue = wa->nodeEquationWaterAge(NodeInputs, VolFlowRates);
        }
        else if(WhichSolver == "biofilm")
        {
          nodalValue = wa->nodeEquationWaterAge(NodeInputs, VolFlowRates);
        }
      }
    }
    else
    {
      if (pool == true)
      {
        pool = false;
        poolFlow = edges.at(elementIndex)->getDoubleProperty("volumeFlowRate");
        for (int j = 0; j < SourceDimension; ++j)
        {
          poolActual[j] = edges.at(elementIndex)->getDoubleProperty(ParameterList[j]);
          poolUpstreamNode[j] = nodes.at(edges.at(elementIndex)->getEdgeIntProperty("startNodeIndex"))->getProperty(ParameterList[j]);
        }
        poolVolumeActual = edges.at(elementIndex)->getDoubleProperty("waterLevel")*edges.at(elementIndex)->getDoubleProperty("referenceCrossSection");
        if (WhichSolver == "waterAge")
        {
          nodalValue = wa->poolEquationWaterAge(poolFlow, poolActual, poolUpstreamNode, poolVolumeActual, h);
        }
        else if (WhichSolver == "chlorine")
        {
          nodalValue = wa->poolEquationWaterAge(poolFlow, poolActual, poolUpstreamNode, poolVolumeActual, h);
        }
        else if(WhichSolver == "biofilm")
        {
          nodalValue = wa->poolEquationWaterAge(poolFlow, poolActual, poolUpstreamNode, poolVolumeActual, h);
        }        
        for (int j = 0; j < SourceDimension; ++j)
        {
          edges.at(elementIndex)->setDoubleProperty(ParameterList[j], nodalValue[j]);
        }
      }
      else if(press == true)
      {
        press = false;
        for (int j = 0; j < SourceDimension; ++j)
        {
          nodalValue[j] = edges.at(elementIndex)->getDoubleProperty(ParameterList[j]);
        }
      }
      else
      {
        cout << "CRITICAL FAILURE!!!" << endl;
      }
    }
    //----------Node update finally--------------//
    for (int j = 0; j < SourceDimension; ++j)
    {
      nodes.at(i)->appendTimeSeries(ParameterList[j], nodalValue[j]);
    }
    nodes.at(i)->appendTimeSeries("timeSteps", nodes.at(i)->getProperty("timeSteps")+h);
    //-------------------------------------------//
  }
}

void QualitySolver::setAdaptiveParameters(double e_min, double e_max, int modAttempts, double dt_inc, double timeMax)
{
  minScaleFactor = e_min;
  maxScaleFactor = e_max;
  attempts = modAttempts; 
  step_time = dt_inc;
  time_max = timeMax;
}

void QualitySolver::solveQuality(string WhichSolver)
{
  solverType = WhichSolver;
  solveSystem();
    
  if (WhichSolver == "waterAge")
  {
    startingCondition.push_back(0.);
    ParameterList = wa->listOfParameters;
    SourceDimension = wa->ModelDimension;
    cout << "------------ WaterAge solver initialized -------------" << endl;
    cout << "The return value of the solver: " << Embedded_Prince_Dormand_solver(time_start, step_time, time_max, tolerance, step_X, WhichSolver) << endl;
    cout << "------------------------------------------------------" << endl;
  } 
  else if (WhichSolver == "chlorine") 
  {
    SourceDimension = ch->ModelDimension;
    cout << "------------ chlorine solver initialized -------------" << endl;
    cout << Embedded_Prince_Dormand_solver(time_start, step_time, time_max, tolerance, step_X, WhichSolver);
    cout << "------------------------------------------------------" << endl;
  }
  else if (WhichSolver == "biofilm")
  {
    SourceDimension = bf->ModelDimension;
    startingCondition.push_back(0.);
    startingCondition.push_back(0.);  
    startingCondition.push_back(0.1);    
    startingCondition.push_back(0.1);
    ParameterList.push_back("Sf");
    ParameterList.push_back("Sb");
    ParameterList.push_back("Cf");
    ParameterList.push_back("Cb");
    cout << "------------ biofilm solver initialized --------------" << endl;
    cout << Embedded_Prince_Dormand_solver(time_start, step_time, time_max, tolerance, step_X, WhichSolver);
    cout << "------------------------------------------------------" << endl;
  }
  else
  {
    cout << "--------------- Error ------------------" << endl;
    cout << "-- please select an available solver  --" << endl;
    cout << "-- e.g. waterAge | chlorine | biofilm --" << endl;
    cout << "----------------------------------------" << endl;
  }
}