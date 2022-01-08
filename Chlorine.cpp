#include "Chlorine.h"

//--------------------------------------------------------------------------//
chlorine::chlorine(){}//
chlorine::~chlorine(){}                                           // 
//--------------------------------------------------------------------------//

vector<double> chlorine::sourceTermChlorine(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX)
{
	vector<double> out;
	for (int i = 0; i < ModelDimension; ++i)
	{
		//cout << "x_act: " << x_actual[i] << " x_upwind: " << x_upwind[i] << endl;
		out.push_back(-1*flow*(x_actual[i]-x_upwind[i])/DX - kb*pow(x_actual[i],n));
		//cout << out[i] << endl; 
	}
    return out;
}

void chlorine::setParameters(std::string which, double value)
{
	if (which == "bulk_coefficient" || which == "kb")
	{
		kb = value; 
	}  
	else if (which == "decay_order" || "n")
	{
		n = value;
	}
	else
	{
		cout << "------------------------- ERROR -------------------------" << endl;
		cout << "- invalid parameter, the available model parameters are -" << endl;
		cout << "------- bulk_coefficient (kb) and decay_order (n) -------" << endl;
	}
}