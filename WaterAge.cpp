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

vector<double> WaterAge::nodeEquationWaterAge(vector<double> &returnNodalValue, vector< vector<double> > nodeInputs, vector<double> VolFlowRates)
{
	double SummarizedVolume = 0.;
	for (int j = 0; j < ModelDimension; ++j)
	{
		returnNodalValue[j] = 0.;
	}
	for (int i = 0; i < nodeInputs.size(); ++i)
	{
		for (int j = 0; j < nodeInputs[i].size(); ++j)
		{
			returnNodalValue[j] += nodeInputs[i][j]*VolFlowRates[i];
		}
		SummarizedVolume += VolFlowRates[i];
	}
	for (int j = 0; j < ModelDimension; ++j)
	{
		returnNodalValue[j] = returnNodalValue[j]/SummarizedVolume;
	}
	return returnNodalValue;
}

vector<double> WaterAge::poolEquationWaterAge(double poolFlow, vector<double> poolActual, vector<double> poolUpstreamNode, double poolVolumeActual, double h)
{
	vector<double> returnNodalValue(ModelDimension);
    if(poolFlow < 0)
    {
		for (int i = 0; i < ModelDimension; ++i)
		{
		returnNodalValue[i] = poolActual[i]+h;
		}
    }
    else
    {
    	for (int i = 0; i < ModelDimension; ++i)
		{
			returnNodalValue[i] = ((poolUpstreamNode[i]*poolFlow*h)+(poolVolumeActual*poolActual[i]))/(poolVolumeActual+poolFlow*h);
		}
    }
    return returnNodalValue;
}