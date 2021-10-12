#include "WaterAge.h"

//--------------------------------------------------------------------------//
WaterAge::WaterAge(){}//
WaterAge::~WaterAge(){}                                           // 
//--------------------------------------------------------------------------//

vector<double> WaterAge::sourceTermWaterAge(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX)
{
	vector<double> out(ModelDimension);
	for (int i = 0; i < ModelDimension; ++i)
	{ 
		out[i] = -1*flow*(x_actual[i]-x_upwind[i])/DX + 1; 
	}
    return out;
}

vector<double> WaterAge::nodeEquationWaterAge(vector< vector<double> > nodeInputs, vector<double> VolFlowRates)
{
	double SummarizedVolume = 0.;
	vector<double> nodalValue(ModelDimension);
	for (int j = 0; j < ModelDimension; ++j)
	{
		nodalValue[j] = 0.;
	}
	for (int i = 0; i < nodeInputs.size(); ++i)
	{
		for (int j = 0; j < nodeInputs[i].size(); ++j)
		{
			nodalValue[j] += nodeInputs[i][j]*VolFlowRates[i];
		}
		SummarizedVolume += VolFlowRates[i];
	}
	for (int j = 0; j < ModelDimension; ++j)
	{
		nodalValue[j] /= SummarizedVolume;
	}
	return nodalValue;
}

vector<double> WaterAge::poolEquationWaterAge(double poolFlow, vector<double> poolActual, vector<double> poolUpstreamNode, double poolVolumeActual, double h)
{
	vector<double> nodalValue(ModelDimension);
    if(poolFlow < 0)
    {
      for (int i = 0; i < ModelDimension; ++i)
      {
        nodalValue[i] = poolActual[i]+h;
      }
    }
    else
    {
    	for (int i = 0; i < ModelDimension; ++i)
		{
			nodalValue[i] = ((poolUpstreamNode[i]*poolFlow*h)+(poolVolumeActual*poolActual[i]))/(poolVolumeActual+poolFlow*h);
		}
    }
    return nodalValue;
}