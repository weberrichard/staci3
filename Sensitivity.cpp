#include "Sensitivity.h"

Sensitivity::Sensitivity(string fileName) : HydraulicSolver(fileName){}
Sensitivity::~Sensitivity(){}

//-----------------------------------------------------------------
bool Sensitivity::calculateSensitivity(string parameter)
{
  bool convergence = solveSystem();
  double rowSum = 0., rowSumLargest = 0.;
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
