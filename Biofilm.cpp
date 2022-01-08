#include "Biofilm.h"

//--------------------------------------------------------------------------//
biofilm::biofilm(){}//
biofilm::~biofilm(){}                                           // 
//--------------------------------------------------------------------------//

vector<double> biofilm::sourceTermBiofilm(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX)
{
	vector<double> out(ModelDimension);
	out[0] = (-1*flow*(x_actual[0]-x_upwind[0])/DX - ((12*Diff/Dcs/W)*(x_actual[0]-x_actual[1])) - (mu_f*x_actual[0]/Ycs/(Kf+x_actual[0]))*x_actual[2]);//kb*pow(x_actual[i],n));
	out[1] = (12*Diff/Dcs/W)*(x_actual[0]-x_actual[1]) - (mu_b*x_actual[1]/Ycs/(Kb+x_actual[1])*x_actual[3]);
	out[2] = -1*flow*(x_actual[2]-x_upwind[2])/DX + mu_f*(1-(x_actual[2]/(C_fm*x_actual[0]+0.001)))*x_actual[2];
	out[3] = mu_b*(1-(x_actual[3]/(C_bm*x_actual[1]+0.001)))*x_actual[3];  
    return out;
}

void biofilm::setParameters(std::string which, double value) 
{
	if (which == "Diff")
	{
		Diff = value;   
	}  
	else if (which == "Dcs")    
	{
		Dcs = value;
	}
	else if (which == "W") 
	{
		W = value;
	}
	else if (which == "mu_f")
	{
		mu_f = value; 
	}
	else if (which == "mu_b")
	{
		mu_b = value;
	}
	else if (which == "Ycs")
	{
		Ycs = value;
	}
	else if (which == "Yf")
	{
		Yf = value;
	}
	else if (which == "Yb")
	{
		Yb = value;
	}
	else if (which == "Kf")
	{
		Kf = value;
	}
	else if (which == "Kb")
	{
		Kb = value;
	}
	else if (which == "C_fm")
	{
		C_fm = value;
	}
	else if (which == "C_bm")  
	{
		C_bm = value;
	}
	else
	{
		cout << "------------------------- ERROR -------------------------" << endl;
		cout << "- invalid parameter, the available model parameters are -" << endl;
		cout << "------- Diff, Dcs, W, mu_f, Ycs, Kf, C_fm, C_bm ---------" << endl;
	}
}

/*vector<double> biofilm::calculateNodalMixture(vector< vector<double> > nodalInputs, vector< vector<double> > inFlows)
{
	vector< vector<double> > nodalOutput;
	double denominator = 0.;
	for (int i = 0; i < nodalInputs.size(); ++i)
	{
		nodalOutput.push_back(0.);
		for (int j = 0; j < nodalInputs[i].size(); ++j)
		{
			for (int k = 0; k < nodalInputs[i][j]; ++k)
			{
				for (int i = 0; i < count; ++i)
				{
					
				}
				nodalOutput[i][k] += nodalInputs[i][j][k]*inFlows[i][j];
				denominator += inFlows;
			}
		}
		nodalOutput[i] /= denominator;
		denominator = 0.;
	}
	return nodalOutput;
}

vector<double> biofilm:calculatePoolParameters*/