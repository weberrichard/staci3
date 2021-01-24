#include "Sensitivity.h"

Sensitivity::Sensitivity(string fileName) : HydraulicSolver(fileName){}
Sensitivity::~Sensitivity(){}

//-----------------------------------------------------------------
bool Sensitivity::calculateSensitivity(string parameter)
{
  bool convergence = solveSystem();

  if(convergence)
  {
    if(parameter == "diameter" || parameter == "roughness")
    {
      int parNumber;
      if(parameter == "roughness")
        parNumber = 0;
      else if(parameter == "diameter")
        parNumber = 1;

      massFlowRateSensitivity.resize(numberEdges,numberEdges);
      pressureSensitivity.resize(numberNodes,numberEdges);
      for(int i=0; i<numberEdges; i++)
      {
        SparseVector<double> funcParDer;
        funcParDer.resize(numberNodes + numberEdges);
        funcParDer.coeffRef(i) = edges[i]->functionParameterDerivative(parNumber);
        VectorXd resultDerivative = solver.solve(-funcParDer);
        massFlowRateSensitivity.col(i) = resultDerivative.head(numberEdges);
        pressureSensitivity.col(i) = resultDerivative.tail(numberNodes);
      }
      for (int i = 0; i < numberEdges; ++i)
      {
        for (int j = 0; j < numberEdges; ++j)
        {
           rowSum += massFlowRateSensitivity(i,j);
        }
        if(abs(rowSum) > abs(rowSumLargest))
        {
          rowSumLargest = rowSum;
          setLargestSensitivity(rowSumLargest);
        }
        cout << "rowSum" << rowSum << endl;
        edges[i]->setMFRSens(rowSum);
        rowSum = 0.;
      }
    }
    else if(strcmp(parameter.c_str(), "demand") == 0)
    {
      massFlowRateSensitivity.resize(numberEdges,numberNodes);
      pressureSensitivity.resize(numberNodes,numberNodes);
      for(int i=0; i<numberNodes; i++)
      {
        SparseVector<double> funcParDer;
        funcParDer.resize(numberNodes + numberEdges);
        funcParDer.coeffRef(numberEdges + i) = nodes[i]->functionParameterDerivative(isPressureDemand);
        VectorXd resultDerivative = solver.solve(-funcParDer);
        massFlowRateSensitivity.col(i) = resultDerivative.head(numberEdges);
        pressureSensitivity.col(i) = resultDerivative.tail(numberNodes);
      }
      for (int i = 0; i < numberNodes; ++i)
      {
        for (int j = 0; j < numberNodes; ++j)
        {
           rowSum += pressureSensitivity(i,j);
        }
        if(abs(rowSum) > abs(rowSumLargest))
        {
          rowSumLargest = rowSum;
          setLargestSensitivity(rowSumLargest);
        }
        nodes[i]->setPsens(rowSum);
        rowSum = 0.;
      }
    }
    else
    {
      cout << "\n\n !!!WARNING!!! calculateSensitivity(string parameter) -> unknown parameter: " << parameter;
      cout << "\n Available values: diameter | friction_coeff | demand";
      cout << "\n Skipping sensitivity calculations, then continouing ..." << endl << endl;
    }

  } // end of if(convergence)
  else
    cout << endl << "[*] Sensitivity (" << parameter << "): hydraulic solver has NOT convergenved :(" << endl;

  return convergence;
}

void Sensitivity::fillEdgeSensitivityForPlot(double OverrideMax)
{
  int numberEdges = edges.size();
  double rowSum, rowSumLargest = OverrideMax;
  cout << "max sens: " << rowSumLargest << endl;
  for (int i = 0; i < numberEdges; ++i)
  {
    rowSum = edges[i]->getMFRSens();
    edges[i]->setMFRSens(rowSum/rowSumLargest);
  }
}

void Sensitivity::fillPipeSensitivityForPlot()
{
  int numberEdges = edges.size();
  double rowSum, rowSumLargest = getLargestSensitivity();
  cout << "max sens: " << rowSumLargest << endl;
  for (int i = 0; i < numberEdges; ++i)
  {
    if(edges[i]->type == "Pipe")
    {
      rowSum = edges[i]->getMFRSens();
      edges[i]->setMFRSens(abs(rowSum/rowSumLargest));
      edges[i]->setDoubleProperty("MFRSensitivity",abs(rowSum/rowSumLargest));
      cout <<"rowSum2: " << edges[i]->getMFRSens() << " : " << edges[i]->getDoubleProperty("MFRSensitivity") << endl;
    }
  }
}

void Sensitivity::fillEdgeSensitivityForPlot()
{
  int numberEdges = edges.size();
  double rowSum, rowSumLargest = getLargestSensitivity();
  cout << "max sens: " << rowSumLargest << endl;
  for (int i = 0; i < numberEdges; ++i)
  {
    rowSum = edges[i]->getMFRSens();
    edges[i]->setMFRSens(abs(rowSum/rowSumLargest));
    cout <<"rowSum2: " <<  edges[i]->getMFRSens() << endl;
  }
}

void Sensitivity::fillNodalSensitivityForPlot()
{
  int numberNodes = nodes.size();
  double rowSum, rowSumLargest = getLargestSensitivity();
  cout << "max sens: " << rowSumLargest << endl;
  for (int i = 0; i < numberNodes; ++i)
  {
    rowSum = nodes[i]->getPsens();
    nodes[i]->setPsens(rowSum/rowSumLargest);
  }
}

void Sensitivity::fillNodalSensitivityForPlotNotNormalized()
{
  int numberNodes = nodes.size();
  double rowSum;
  //cout << "max sens: " << rowSumLargest << endl;
  for (int i = 0; i < numberNodes; ++i)
  {
    rowSum = nodes[i]->getPsens();
    nodes[i]->setPsens(rowSum);
  }
}

double Sensitivity::CalculateAverageSensitivity()
{
  int numberNodes = nodes.size();
  double rowsum = 0.;
  double nodal = 0.;
  double AverageLocal;
  for (int i = 0; i < numberNodes; ++i)
  {
    rowsum = nodes[i]->getPsens();
    //cout << "node " << i << ". loc. sens: " << rowsum << endl; 
    nodal += rowsum;
  }
  AverageLocal = nodal/numberNodes;
  //cout << " Av. Local: " << AverageLocal << endl;
  return AverageLocal;
}

void Sensitivity::fillNodalSensitivityForPlot(double OverrideMax)
{
  int numberNodes = nodes.size();
  double rowSum, rowSumLargest = OverrideMax;
  for (int i = 0; i < numberNodes; ++i)
  {
    rowSum = nodes[i]->getPsens();
    nodes[i]->setPsens(rowSum/rowSumLargest);
  }
}