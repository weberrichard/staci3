#include "HydraulicSolver.h"

//--------------------------------------------------------------
HydraulicSolver::HydraulicSolver(string spr_filename) : Staci(spr_filename)
{
  maxIterationNumber = 100;
  maxPressureError = 1.524e-4*5.;
  maxMassFlowError = 2.832e-5*5.;
  relaxationFactor = 1.0;
  relaxationFactorIncrement = 1.1;
  relaxationFactorDecrement = 0.8;

  volumeFlowRateTolerance = 2.832e-5; // 1e-3 cfs
  headTolerance = 1.524e-4; // 5e-4 feet

  // setting the number of equations
  numberEquations = numberEdges + numberNodes;

  // Setting the heights of the edges
  for(int i=0; i<numberEdges; i++)
  {
    if(edges[i]->numberNode == 2)
    {
      int startIndex = edges[i]->startNodeIndex;
      edges[i]->startHeight = nodes[startIndex]->geodeticHeight;

      int endIndex = edges[i]->endNodeIndex;
      edges[i]->endHeight = nodes[endIndex]->geodeticHeight;
    }
    if(edges[i]->numberNode == 1)
    {
      int startIndex = edges[i]->startNodeIndex;
      edges[i]->startHeight = nodes[startIndex]->geodeticHeight;
    }
  }

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}

HydraulicSolver::~HydraulicSolver(){}

//--------------------------------------------------------------
bool HydraulicSolver::solveSystem()
{
  double e_mp, e_p, e_mp_r, e_p_r;
  int iter;
  bool isConv;
  vector<int> changedIndex;
  
  e_mp = 1e10, e_p = 1e10, e_mp_r = 1e10, e_p_r = 1e10;
  iter = 0;
  isConv = false;

  updateJacobian();

  computeError(f, e_mp, e_p, e_mp_r, e_p_r, isConv);

  if(printLevel > 2)
    iterInfo(iter, e_mp, e_p, changedIndex);

  while((iter<maxIterationNumber+1) && !isConv)
  {
    //checkJacobianMatrix();
    linearSolver(); // Solving Jac*dx = f

    changedIndex = edgeStatusUpdate();
    if(changedIndex.size() > 0)
      isConv = false;

    updateJacobian();
    computeError(f, e_mp, e_p, e_mp_r, e_p_r, isConv);

    if(printLevel > 2)
      iterInfo(iter+1, e_mp, e_p, changedIndex);
    if(printLevel > 3)
      iterInfoDetail();

    if(isConv)
    {
      changedIndex = edgeStatusUpdate(); // last check after convergence
      if(changedIndex.size() > 0) // edges status change found
      {
        updateJacobian();
        isConv = false;
      }
      else // everything converged and statuses not changed
      {
        break;
      }
    }
    
    if(iter>1)
      updateRelaxationFactor(e_mp, e_p, e_mp_r, e_p_r);

    iter++;
  }

  if(iter > maxIterationNumber)
  {
    iterInfo(iter, e_mp, e_p, changedIndex);
    iterInfoDetail();
    isConv == false;
  }
  else
  {
    for(int i=0; i<numberEdges; i++)
    {
      edges[i]->volumeFlowRate = x(i);
    }
    for(int i=0; i<numberNodes; i++)
    { 
      nodes[i]->head = x(numberEdges + i);
    }
    relaxationFactor = 1.0;
  }

  // Checking the pumps
  checkPumpOperatingPoint();
  
  return isConv;
}

//--------------------------------------------------------------
void HydraulicSolver::linearSolver()
{
  VectorXd dx;
  dx.resize(numberEquations);
  solver.factorize(jacobianMatrix); // performing LU decomposition
  //cout << endl << "EIGEN FACTORIZE: " << solver.info() << endl;
  //cout << endl << "JAC DET        : " << solver.absDeterminant() << endl;
  //cout << endl << solver.lastErrorMessage() << endl;
  dx = solver.solve(-f); // Solving Jac*dx = -f
  x += relaxationFactor*dx; // x_i+1 = x_i + relax*dx (~0.1<relax<~1.0)
}

//--------------------------------------------------------------
void HydraulicSolver::buildJacobian()
{
  // calculating the max rank of the jacobian matrix
  maxRank = 4;
  for(int i=0; i<numberNodes; i++)
  {
    if(nodes[i]->edgeIn.size() + nodes[i]->edgeOut.size() > maxRank)
    {
      maxRank = nodes[i]->edgeIn.size() + nodes[i]->edgeOut.size();
    }
  }

  int n = x.rows();
  jacobianMatrix.resize(n,n);
  jacobianMatrix.reserve(VectorXi::Constant(n,maxRank));

  for(int i=0; i<numberEdges; i++)
  {
    VectorXd funcDer = getEdgeFunction(i);

    int startIndex, endIndex;
    startIndex = edges[i]->startNodeIndex;

    jacobianMatrix.insert(i,numberEdges + startIndex) = funcDer(0);
    if(edges[i]->numberNode == 2)
    {
      endIndex = edges[i]->endNodeIndex;
      jacobianMatrix.insert(i,numberEdges + endIndex) = funcDer(1);
    }
    jacobianMatrix.insert(i,i) = funcDer(2);
  }

  for(int i=0; i<numberNodes; i++)
  {
    int nIn = nodes[i]->edgeIn.size(), nOut = nodes[i]->edgeOut.size();
    VectorXd funcDer = getNodeFunction(i);

    jacobianMatrix.insert(numberEdges + i, numberEdges + i) = funcDer(0);
    for(int j=0; j<nIn; j++)
    {
      jacobianMatrix.insert(numberEdges + i, nodes[i]->edgeIn[j]) = funcDer(1+j);
    }
    for(int j=0; j<nOut; j++)
    {
      jacobianMatrix.insert(numberEdges + i, nodes[i]->edgeOut[j]) = funcDer(1+nIn+j);
    }
  }
  solver.analyzePattern(jacobianMatrix);
}

//--------------------------------------------------------------
void HydraulicSolver::updateJacobian()
{
  for(int i=0; i<numberEdges; i++)
  {
    VectorXd funcDer = getEdgeFunction(i);

    int startIndex, endIndex;
    startIndex = edges[i]->startNodeIndex;

    jacobianMatrix.coeffRef(i,numberEdges + startIndex) = funcDer(0);
    if(edges[i]->numberNode == 2)
    {
      endIndex = edges[i]->endNodeIndex;
      jacobianMatrix.coeffRef(i,numberEdges + endIndex) = funcDer(1);
    }
    jacobianMatrix.coeffRef(i,i) = funcDer(2);
  }

  for(int i=0; i<numberNodes; i++)
  { 
    int nIn = nodes[i]->edgeIn.size(), nOut = nodes[i]->edgeOut.size();
    VectorXd funcDer = getNodeFunction(i);

    jacobianMatrix.coeffRef(numberEdges + i, numberEdges + i) = funcDer(0);
    for(int j=0; j<nIn; j++)
    {
      jacobianMatrix.coeffRef(numberEdges + i, nodes[i]->edgeIn[j]) = funcDer(1+j);
    }
    for(int j=0; j<nOut; j++)
    {
      jacobianMatrix.coeffRef(numberEdges + i, nodes[i]->edgeOut[j]) = funcDer(1+nIn+j);
    }
  }
}

//--------------------------------------------------------------
VectorXd HydraulicSolver::getEdgeFunction(int i)
{
  int startIndex, endIndex;
  startIndex = edges[i]->startNodeIndex;

  VectorXd ppq = VectorXd::Constant(3, 0.0); // contains [Pstart,Pend,Vf]
  ppq(0) = x(numberEdges + startIndex);
  if(edges[i]->numberNode == 2)
  {
    endIndex = edges[i]->endNodeIndex;
    ppq(1) = x(numberEdges + endIndex);
  }
  ppq(2) = x(i);

  VectorXd funcDer = VectorXd::Constant(3, 0.0);
  f(i) = edges[i]->function(ppq, funcDer);

  return funcDer;
}

//--------------------------------------------------------------
VectorXd HydraulicSolver::getNodeFunction(int i)
{
  int nIn = nodes[i]->edgeIn.size(), nOut = nodes[i]->edgeOut.size();
  VectorXd pq = VectorXd::Constant(1 + nIn + nOut, 0.0);

  pq(0) = x(numberEdges + i); // pq = [p, Qin1, Qin2, ..., Qout1, Qout2, ...]
  for(int j=0; j<nIn; j++)
  {
    pq(1+j) = x(nodes[i]->edgeIn[j]);
  }
  for(int j=0; j<nOut; j++)
  {
    pq(1+nIn+j) = x(nodes[i]->edgeOut[j]);
  }

  // calling the functions and function derivatives
  VectorXd funcDer = VectorXd::Constant(1 + nIn + nOut, 0.0);
  f(numberEdges + i) = nodes[i]->function(pq,isPressureDemand,isLeakage,funcDer);

  return funcDer;
}


//--------------------------------------------------------------
void HydraulicSolver::updateRelaxationFactor(double e_mp, double e_p, double & e_mp_r, double & e_p_r)
{
  double error, error_old;
  if (e_p > e_mp) {
    error = e_p;
    error_old = e_p_r;
  } else {
    error = e_mp;
    error_old = e_mp_r;
  }
  if (error <= error_old)
    relaxationFactor = relaxationFactor * relaxationFactorIncrement;
  else
    relaxationFactor = relaxationFactor * relaxationFactorDecrement;
  if (relaxationFactor < minRelaxationFactor)
    relaxationFactor = minRelaxationFactor;
  if (relaxationFactor > maxRelaxationFactor)
    relaxationFactor = maxRelaxationFactor;

  e_mp_r = e_mp;
  e_p_r = e_p;
}

//--------------------------------------------------------------
void HydraulicSolver::computeError(const VectorXd &f, double & e_mp, double & e_p, double & e_mp_r, double & e_p_r, bool & convergence)
{
  // Saving old errors
  e_mp_r = e_mp;
  e_p_r = e_p;

  // Determining new errors
  e_p = (f.head(numberEdges)).norm();
  e_mp = (f.tail(numberNodes)).norm();

  // Checking the convergence
  if ((e_p < maxPressureError) && (e_mp < maxMassFlowError)) 
    convergence = true;
}

//--------------------------------------------------------------
void HydraulicSolver::initialization()
{
  // Finding an inital base geodetic level of a pressure point or a pool
  double baseGeodeticHeight;
  if(presIndex.size() > 0)
  {
    int idx = edges[presIndex[0]]->startNodeIndex;
    baseGeodeticHeight = nodes[idx]->geodeticHeight;
  }
  else if(poolIndex.size() > 0)
  {
    int idx = edges[poolIndex[0]]->startNodeIndex;
    baseGeodeticHeight = nodes[idx]->geodeticHeight;
  }
  else
  {
    cout << endl << "!WARNING! No PressurePoint, nor Pool has been found in the network! \n EXITING! " << endl;
    exit(-1);
  }

  // Filling up the nodes
  for(int i = 0; i < numberNodes; i++)
  {
    double head = 50. + baseGeodeticHeight;
    nodes[i]->initialization(0., head);
  }

  // Filling up the edges
  for(int i=0; i<pipeIndex.size(); i++)
  {
    double v = 0.2;
    double D = edges[pipeIndex[i]]->getDoubleProperty("diameter");
    edges[pipeIndex[i]]->initialization(1, v*D*D*M_PI/4.);
  }
  for(int i=0; i<pumpIndex.size(); i++)
  {
    int idx = pumpIndex[i];
    int pumpType = edges[idx]->getIntProperty("pumpType");
    if(pumpType >= 0)
    {
      vector<double> qCurve = edges[idx]->getVectorProperty("qCurve");
      double min = *min_element(qCurve.begin(),qCurve.end());
      double max = *max_element(qCurve.begin(),qCurve.end());
      edges[idx]->initialization(1, (min+max)/2.);
    }
    else if(pumpType == -1)
    {
      edges[idx]->initialization(1, 0.1);
    }
  }
  for(int i=0; i<valveIndex.size(); i++)
  {
    if(edges[valveIndex[i]]->typeCode == 6) // FCV
    {
      edges[valveIndex[i]]->initialization(1, edges[valveIndex[i]]->setting*1.1);
    }
    else
    {
      edges[valveIndex[i]]->initialization(1, 1.);
    }
  }

  // Fill up the outflows (pressure points, pools)
  double sumCons = 0.;
  for(int i=0; i<numberNodes; i++)
    sumCons += nodes[i]->demand;
  int np = presIndex.size() + poolIndex.size();
  for(int i=0; i<presIndex.size(); i++)
  {
    edges[presIndex[i]]->initialization(1, sumCons/np);
  }
  for(int i=0; i<poolIndex.size(); i++)
  {
    edges[poolIndex[i]]->initialization(1, sumCons/np);
  }
}

//--------------------------------------------------------------
void HydraulicSolver::iterInfo(int iter, double e_mp, double e_p, vector<int> idx)
{
  printf(" iter #%3i/%4i | e_mp=%6.3e | e_p=%6.3e | relax=%5.2f | status changes:%3zu", iter, maxIterationNumber, e_mp, e_p, relaxationFactor, idx.size());
  cout << endl;
}

//--------------------------------------------------------------
void HydraulicSolver::iterInfoDetail()
{   
  printf("\n==========================================================");
  printf("\n   WORSTLY CONVERVGED ELEMENTS    ");
  printf("\n=====================+====================================");
  vector<int> idx;
  vector<double> fs = findMaxValues(f.head(numberEdges), idx, 10, true);
  printf("\n      Edge name      |      f     |  vf [m3/s]  | status |");
  printf("\n=====================+====================================");
  for(int i=0; i<fs.size(); i++)
  {
    printf("\n%20s | %10.3e | %10.3e | %6i |",edges[idx[i]]->name.c_str(), fs[i], x(idx[i]), edges[idx[i]]->status);
  }
  printf("\n=====================+====================================");

  fs = findMaxValues(f.tail(numberNodes), idx, 10, true);
  printf("\n      Node name      |     f    |  p  [m]  ");
  printf("\n=====================+====================================");
  for(int i=0; i<fs.size(); i++)
  {
    printf("\n%20s | %8.3e | %8.3f ",nodes[idx[i]]->name.c_str(), fs[i], x(numberEdges + idx[i]));
  }
  printf("\n=====================+====================================");
  cout << endl;
}

//--------------------------------------------------------------
void HydraulicSolver::listResult()
{
  cout << scientific << setprecision(3) << showpos;
  cout << endl << endl << "RESULTS:";
  for(int i = 0; i < numberEdges; i++){
    printf("\n %-12s:",edges[i]->name.c_str());
    cout << "  mp = " << edges[i]->volumeFlowRate << " kg/s" << "   Q = " << edges[i]->volumeFlowRate << " l/s" << "   v = " << edges[i]->getDoubleProperty("velocity") << " m/s";
  }
  cout << endl << "\t" << "-----------------------------------------------------------------";
  for(int i = 0; i < numberNodes; i++){
    printf("\n %-12s:", nodes[i]->name.c_str());
    cout << "  p = " << nodes[i]->head * 1000. * 9.81 / 1.e5 << " bar" << "     H = " << nodes[i]->head << " m" << "      H+height=" << nodes[i]->head + nodes[i]->geodeticHeight << " m";
  }
  cout << endl << endl;
}

//--------------------------------------------------------------
void HydraulicSolver::checkPumpOperatingPoint()
{
  for(int i=0; i<pumpIndex.size(); i++)
  {
    edges[i]->checkPump();
  }
}

//--------------------------------------------------------------
vector<int> HydraulicSolver::edgeStatusUpdate()
{
  vector<int> changedIndex;
  // checking the pipes with check valves
  for(int i=0; i<pipeCVIndex.size(); i++)
  {
    int idx = pipeCVIndex[i]; // index of the pipe with cv in the edges list
    int ii = edges[idx]->startNodeIndex;
    int ij = edges[idx]->endNodeIndex;
    double pi = x(numberEdges + ii);
    double pj = x(numberEdges + ij);
    double hi = nodes[ii]->geodeticHeight;
    double hj = nodes[ij]->geodeticHeight;
    double vf = x(idx);
    int status = edges[idx]->status;

    if(abs(pi+hi - pj-hj) > headTolerance)
    {
      if(pi+hi - pj-hj < -headTolerance)
      {
        edges[idx]->status = 0;
      }
      else if(vf < -volumeFlowRateTolerance)
      {
        edges[idx]->status = 0;
      }
      else
      {
        edges[idx]->status = 1;
      }
    }
    else
    {
      if(vf < -volumeFlowRateTolerance)
      {
        edges[idx]->status = 0;
      }
    }

    if(status != edges[idx]->status)
    {
      changedIndex.push_back(idx);
      if(edges[idx]->status == 0)
      {
        x(idx) = 0.0; // setting new "initial" condition
        closeIsolatedNode(ii);
        closeIsolatedNode(ij); 
      }
      else
      {
        nodes[edges[idx]->startNodeIndex]->status = 1;
        nodes[edges[idx]->endNodeIndex]->status = 1;
      }
    }
  }

  for(int i=0; i<valveIndex.size(); i++)
  {
    int idx = valveIndex[i]; // index of the valves in the edges list
    int ii = edges[idx]->startNodeIndex;
    int ij = edges[idx]->endNodeIndex;
    double pi = x(numberEdges + ii);
    double pj = x(numberEdges + ij);
    double hi = nodes[ii]->geodeticHeight;
    double hj = nodes[ij]->geodeticHeight;
    double vf = x(idx);
    int old_status = edges[idx]->status;
    double setting = edges[idx]->setting;
    if(edges[idx]->typeCode == 3) // PRV: Pressure Reduction Valve
    {
      if(old_status == 2) // ACTIVE
      {
        if(vf < -volumeFlowRateTolerance)
        {
          edges[idx]->status = 0;
        } 
        else if(pi + hi < setting + hj - headTolerance)// todo: minor loss
        {
          edges[idx]->status = 1;
        }
        else
        {
          edges[idx]->status = 2;
        }
      }
      else if(old_status == 1) // OPEN
      {
        if(vf < -volumeFlowRateTolerance)
        {
          edges[idx]->status = 0;
        } 
        else if(pj + hj > setting + hj + headTolerance)// todo: minor loss
        {
          edges[idx]->status = 2;
        }
        else
        {
          edges[idx]->status = 1;
        }
      }
      else if(old_status == 0) // CLOSED
      {
        if(pi + hi > pj + hj + headTolerance && pi + hi < setting + hj - headTolerance)
        {
          edges[idx]->status = 1;
        }
        else if(pi + hi > pj + hj + headTolerance && pj + hj < setting + hj - headTolerance)
        {
          edges[idx]->status = 2;
        }
        else
        {
          edges[idx]->status = 0;
        }
      }

      if(old_status != edges[idx]->status)
      {
        changedIndex.push_back(idx);
        if(edges[idx]->status == 0)
        {
          x(idx) = 0.0;
          closeIsolatedNode(ii);
          closeIsolatedNode(ij);
        }
        else
        {
          nodes[edges[idx]->startNodeIndex]->status = 1;
          nodes[edges[idx]->endNodeIndex]->status = 1;
          if(edges[idx]->status == 2)
          {
            x(numberEdges + ij) = setting;
          }
        }
      }
    }
    else if(edges[idx]->typeCode == 4) // PSV: Pressure Sustaining Valve
    {
      if(old_status == 2) // ACTIVE
      {
        if(vf < -volumeFlowRateTolerance)
        {
          edges[idx]->status = 0;
        } 
        else if(pj + hj > setting + hi + headTolerance)// todo: minor loss
        {
          edges[idx]->status = 1;
        }
        else
        {
          edges[idx]->status = 2;
        }
      }
      else if(old_status == 1) // OPEN
      {
        if(vf < -volumeFlowRateTolerance)
        {
          edges[idx]->status = 0;
        } 
        else if(pi + hi < setting + hi - headTolerance)// todo: minor loss
        {
          edges[idx]->status = 2;
        }
        else
        {
          edges[idx]->status = 1;
        }
      }
      else if(old_status == 0) // CLOSED
      {
        if(pj + hj > setting + hi + headTolerance && pi + hi > pj + hj + headTolerance)
        {
          edges[idx]->status = 1;
        }
        else if(pi + hi > setting + hi + headTolerance && pi + hi > pj + hj + headTolerance)
        {
          edges[idx]->status = 2;
        }
        else
        {
          edges[idx]->status = 0;
        }
      }

      if(edges[idx]->status != old_status)
      {
        changedIndex.push_back(idx);

        // setting new "initial" condition
        if(edges[idx]->status == 0)
          x(idx) = 0.0; 
        else if(edges[idx]->status == 2)
          x(numberEdges + ii) = setting;
      }
    }
    else if(edges[idx]->typeCode == 6) // FCV: Flow Control Valve
    {
      if(pi + hi - pj - hj < -headTolerance)
      {
        edges[idx]->status = 1; // OPEN
      }
      else if(vf < -volumeFlowRateTolerance)
      {
        edges[idx]->status = 1; // OPEN
      }
      else if(old_status == 1 && vf >= setting)
      {
        edges[idx]->status = 2; // ACTIVE
      }

      if(edges[idx]->status != old_status)
      {
        changedIndex.push_back(idx);
        
        // setting new "initial" condition
        if(edges[idx]->status == 0)
        {
          x(idx) = 0.0; 
          closeIsolatedNode(ii);
          closeIsolatedNode(ij);
        }
        else
        {
          nodes[edges[idx]->startNodeIndex]->status = 1;
          nodes[edges[idx]->endNodeIndex]->status = 1;
          if(edges[idx]->status == 2)
          {
            x(numberEdges + ij) = setting;
          }
        }
      }
    }
  }

  if(printLevel > 1)
  {
    if(changedIndex.size() > 0)
    {
      printf(" ************************************* CHANGE BY HYDRAULICS *************************************\n");
      printf(" |       Edge ID        |   type   | status |  setting   | vf [m3/s]  |   p1 [m]   |    p0 [m]  |\n");
      for(int i=0; i<changedIndex.size(); i++)
      {
        int idx = changedIndex[i];
        printf(" | %-20s | %8s |  %3i   | %10.3e | %10.3e | %10.3e | %10.3e | \n", edges[idx]->name.c_str(), edges[idx]->type.c_str(), edges[idx]->status, edges[idx]->setting, x(idx), x(numberEdges + edges[idx]->endNodeIndex), x(numberEdges + edges[idx]->startNodeIndex));
      }
      cout << endl;
    }
  }

  return changedIndex;
}

//--------------------------------------------------------------
void HydraulicSolver::closeIsolatedNode(int idx)
{
  int n_in = nodes[idx]->edgeIn.size();
  int n_out = nodes[idx]->edgeOut.size();
  bool node_active = false;

  for(int i=0; i<n_in; i++)
  {
    if(edges[nodes[idx]->edgeIn[i]]->status >= 1)
      node_active += true;
  }
  for(int i=0; i<n_out; i++)
  {
    if(edges[nodes[idx]->edgeOut[i]]->status >= 1)
      node_active += true;
  }

  if(node_active)
    nodes[idx]->status = 1;
  else
    nodes[idx]->status = 0;
}

//--------------------------------------------------------------
void HydraulicSolver::checkJacobianMatrix()
{
  bool ok = true;
  for(int i=0; i<jacobianMatrix.cols(); i++)
  {
    if(jacobianMatrix.col(i).norm() == 0)
    {
      ok = false;
      cout << "!WARNING! There are only zeros in col " << i;
      if(i<numberEdges)
        cout << "  edge name: " << edges[i]->name << endl;
      else
        cout << "  node name: " << nodes[i-numberEdges]->name << endl;
    }
  }

  if(ok)
    cout << endl << " There is no column in the jacobianMatrix with zero norm. " << endl;
  else
    exit(-1);
}


//--------------------------------------------------------------
void HydraulicSolver::addNewPipe(string name, string nodeFrom, string nodeTo, double density, double length, double diameter, double roughness, double volumeFlowRate, bool isCheckVavle, int frictionModel, double relativeViscosity)
{
  edges.push_back(new Pipe(name, nodeFrom, nodeTo, density, length, diameter, roughness, volumeFlowRate, isCheckVavle, frictionModel, relativeViscosity));

  numberNodes = nodes.size();
  numberEdges = edges.size();
  numberEquations = numberEdges + numberNodes;
  
  buildSystem();
  buildIndexing();

  int n=edges.size();
  edges[n-1]->startHeight = nodes[edges[n-1]->startNodeIndex]->geodeticHeight;
  edges[n-1]->endHeight = nodes[edges[n-1]->endNodeIndex]->geodeticHeight;

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}

//--------------------------------------------------------------
void HydraulicSolver::addNewISOValves(vector<string> valveName, vector<string> pipeName, vector<bool> isStart, double density, vector<double> referenceCrossSection, double volumeFlowRate)
{ 
  vector<int> pipeISO;
  for(unsigned int i=0; i<pipeName.size(); i++)
  {
    int pi = edgeIDtoIndex(pipeName[i]); // pipe index
    pipeISO.push_back(pi);
  }

  addNewISOValves(valveName, pipeISO, isStart, density, referenceCrossSection, volumeFlowRate);
}

//--------------------------------------------------------------
void HydraulicSolver::addNewISOValves(vector<string> valveName, vector<int> pipeISO, vector<bool> isStart, double density, vector<double> referenceCrossSection, double volumeFlowRate)
{

  for(unsigned int i=0; i<valveName.size(); i++)
  {
    int pi = pipeISO[i]; // pipe index
    int ns = nodeIDtoIndex(edges[pi]->startNodeName); // start node index
    int ne = nodeIDtoIndex(edges[pi]->startNodeName); // end node index
    double x = (nodes[ns]->xPosition + nodes[ne]->xPosition)*.5; // coordinates of additional node
    double y = (nodes[ns]->yPosition + nodes[ne]->yPosition)*.5;
    string newNode = valveName[i] + "_N";
    if(isStart[i])
    {
      string sn = edges[pi]->startNodeName;
      edges[pi]->startNodeName = newNode;
      nodes.push_back(new Node(newNode, x, y, edges[pi]->startHeight, 0.0, 0.0, density));
      edges.push_back(new ValveISO(valveName[i], sn, newNode, density, referenceCrossSection[i], volumeFlowRate));
      int n=edges.size();
      edges[n-1]->startHeight = edges[pi]->startHeight;
      edges[n-1]->endHeight = edges[pi]->startHeight;
      numberEdges++;
      numberNodes++;
    }
    else
    {
      string en = edges[pi]->endNodeName;
      edges[pi]->endNodeName = newNode;
      nodes.push_back(new Node(newNode, x, y, edges[pi]->endHeight, 0.0, 0.0, density));
      edges.push_back(new ValveISO(valveName[i], en, newNode, density, referenceCrossSection[i], volumeFlowRate));
      int n=edges.size();
      edges[n-1]->startHeight = edges[pi]->endHeight;
      edges[n-1]->endHeight = edges[pi]->endHeight;
      numberEdges++;
      numberNodes++;
    }
  }

  numberEquations = numberEdges + numberNodes;

  buildSystem();
  buildIndexing();

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}
//--------------------------------------------------------------
void HydraulicSolver::addNewPRVValves(vector<string> valveName, vector<string> pipeName, vector<bool> isStart, double density, vector<double> referenceCrossSection, vector<double> PRVSetting, vector<double> Losses, double volumeFlowRate)
{ 
  vector<int> pipeISO;
  for(unsigned int i=0; i<pipeName.size(); i++)
  {
    int pi = edgeIDtoIndex(pipeName[i]); // pipe index
    pipeISO.push_back(pi);
  }

  addNewPRVValves(valveName, pipeISO, isStart, density, referenceCrossSection, PRVSetting, Losses, volumeFlowRate);
}

//--------------------------------------------------------------
void HydraulicSolver::addNewPRVValves(vector<string> valveName, vector<int> pipeISO, vector<bool> isStart, double density, vector<double> referenceCrossSection, vector<double> PRVSetting, vector<double> Losses, double volumeFlowRate)
{

  for(unsigned int i=0; i<valveName.size(); i++)
  {
    int pi = pipeISO[i]; // pipe index
    int ns = nodeIDtoIndex(edges[pi]->startNodeName); // start node index
    int ne = nodeIDtoIndex(edges[pi]->startNodeName); // end node index
    double x = (nodes[ns]->xPosition + nodes[ne]->xPosition)*.5; // coordinates of additional node
    double y = (nodes[ns]->yPosition + nodes[ne]->yPosition)*.5;
    string newNode = valveName[i] + "_N";
    if(isStart[i])
    {
      string sn = edges[pi]->startNodeName;
      edges[pi]->startNodeName = newNode;
      nodes.push_back(new Node(newNode, x, y, edges[pi]->startHeight, 0.0, 0.0, density));
      edges.push_back(new ValvePRV(valveName[i], sn, newNode, density, referenceCrossSection[i], PRVSetting[i], Losses[i], volumeFlowRate));
      int n=edges.size();
      edges[n-1]->startHeight = edges[pi]->startHeight;
      edges[n-1]->endHeight = edges[pi]->startHeight;
      numberEdges++;
      numberNodes++;
    }
    else
    {
      string en = edges[pi]->endNodeName;
      edges[pi]->endNodeName = newNode;
      nodes.push_back(new Node(newNode, x, y, edges[pi]->endHeight, 0.0, 0.0, density));
      edges.push_back(new ValvePRV(valveName[i], en, newNode, density, referenceCrossSection[i], PRVSetting[i], Losses[i], volumeFlowRate));
      int n=edges.size();
      edges[n-1]->startHeight = edges[pi]->endHeight;
      edges[n-1]->endHeight = edges[pi]->endHeight;
      numberEdges++;
      numberNodes++;
    }
  }

  numberEquations = numberEdges + numberNodes;

  buildSystem();
  buildIndexing();

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}
//--------------------------------------------------------------
void HydraulicSolver::addNewPressurePoint(string name, double referenceCrossSection, string startNodeName, double density, double head, double totalHead, double volumeFlowRate)
{
  edges.push_back(new PressurePoint(name, referenceCrossSection, startNodeName, density, head, totalHead, volumeFlowRate));

  numberNodes = nodes.size();
  numberEdges = edges.size();
  numberEquations = numberEdges + numberNodes;
  
  buildSystem();
  buildIndexing();
  int n=edges.size();
  edges[n-1]->startHeight = nodes[edges[n-1]->startNodeIndex]->geodeticHeight;
  //edges[n-1]->endHeight = nodes[edges[n-1]->endNodeIndex]->geodeticHeight;
  cout <<  " megvan44 " << endl;

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}

//--------------------------------------------------------------
void HydraulicSolver::deleteISOValves(vector<string> valveName)
{
  vector<int> valveIndex;
  for(unsigned int i=0; i<valveName.size(); i++)
  { 
    int idx = edgeIDtoIndex(valveName[i]);
    valveIndex.push_back(idx);
  }

  deleteISOValves(valveIndex);
}

//--------------------------------------------------------------
void HydraulicSolver::deleteISOValves(vector<int> valveIndex)
{
  vector<string> nodesToDelete;
  vector<string> changeToNodes;
  vector<double> nodeHeight;
  for(unsigned int i=0; i<valveIndex.size(); i++)
  { 
    int idx = valveIndex[i];
    nodesToDelete.push_back(edges[idx]->endNodeName);
    changeToNodes.push_back(edges[idx]->startNodeName);
    idx = nodeIDtoIndex(nodesToDelete[i]);
    nodeHeight.push_back(nodes[idx]->geodeticHeight);
    if(nodes[idx]->demand != 0.)
    {
      cout << endl << " !WARNING! Node " << nodes[idx]->name << " is to be deleted, while having " << nodes[idx]->demand << " demand" << endl;
    }
  }

  // deleting edges
  for(unsigned int i=0; i<edges.size(); i++)
  {
    for(unsigned int j=0; j<valveIndex.size(); j++)
    {
      if(i == valveIndex[j])
      {
        edges.erase(edges.begin()+i);
        i--;
      }
    }
  }

  // deleting unnecessary nodes
  for(unsigned int i=0; i<nodes.size(); i++)
  {
    for(unsigned int j=0; j<nodesToDelete.size(); j++)
    {
      if(nodes[i]->name == nodesToDelete[j])
      {
        nodes.erase(nodes.begin()+i);
        i--;
      }
    }
  }

  // changing deleted nodes
  for(unsigned int i=0; i<edges.size(); i++)
  {
    for(unsigned int j=0; j<nodesToDelete.size(); j++)
    {
      if(edges[i]->startNodeName == nodesToDelete[j])
      {
        edges[i]->startNodeName = changeToNodes[j];
        edges[i]->startHeight = nodeHeight[j];
      }
      if(edges[i]->endNodeName == nodesToDelete[j])
      {
        edges[i]->endNodeName = changeToNodes[j];
        edges[i]->endHeight = nodeHeight[j];
      }
    }
  }

  numberNodes = nodes.size();
  numberEdges = edges.size();
  numberEquations = numberEdges + numberNodes;

  buildSystem();

  buildIndexing();

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}
//--------------------------------------------------------------
void HydraulicSolver::deletePRVValves(vector<string> valveName)
{
  vector<int> valveIndex;
  for(unsigned int i=0; i<valveName.size(); i++)
  { 
    int idx = edgeIDtoIndex(valveName[i]);
    valveIndex.push_back(idx);
  }

  deleteISOValves(valveIndex);
}

//--------------------------------------------------------------
void HydraulicSolver::deletePRVValves(vector<int> valveIndex)
{
  vector<string> nodesToDelete;
  vector<string> changeToNodes;
  vector<double> nodeHeight;
  for(unsigned int i=0; i<valveIndex.size(); i++)
  { 
    int idx = valveIndex[i];
    nodesToDelete.push_back(edges[idx]->endNodeName);
    changeToNodes.push_back(edges[idx]->startNodeName);
    idx = nodeIDtoIndex(nodesToDelete[i]);
    nodeHeight.push_back(nodes[idx]->geodeticHeight);
    if(nodes[idx]->demand != 0.)
    {
      cout << endl << " !WARNING! Node " << nodes[idx]->name << " is to be deleted, while having " << nodes[idx]->demand << " demand" << endl;
    }
  }

  // deleting edges
  for(unsigned int i=0; i<edges.size(); i++)
  {
    for(unsigned int j=0; j<valveIndex.size(); j++)
    {
      if(i == valveIndex[j])
      {
        edges.erase(edges.begin()+i);
        i--;
      }
    }
  }

  // deleting unnecessary nodes
  for(unsigned int i=0; i<nodes.size(); i++)
  {
    for(unsigned int j=0; j<nodesToDelete.size(); j++)
    {
      if(nodes[i]->name == nodesToDelete[j])
      {
        nodes.erase(nodes.begin()+i);
        i--;
      }
    }
  }

  // changing deleted nodes
  for(unsigned int i=0; i<edges.size(); i++)
  {
    for(unsigned int j=0; j<nodesToDelete.size(); j++)
    {
      if(edges[i]->startNodeName == nodesToDelete[j])
      {
        edges[i]->startNodeName = changeToNodes[j];
        edges[i]->startHeight = nodeHeight[j];
      }
      if(edges[i]->endNodeName == nodesToDelete[j])
      {
        edges[i]->endNodeName = changeToNodes[j];
        edges[i]->endHeight = nodeHeight[j];
      }
    }
  }

  numberNodes = nodes.size();
  numberEdges = edges.size();
  numberEquations = numberEdges + numberNodes;

  buildSystem();

  buildIndexing();

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}