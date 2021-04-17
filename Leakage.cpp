#include "Leakage.h"

//--------------------------------------------------------------
Leakage::Leakage(string spr_filename) : HydraulicSolver(spr_filename){}
Leakage::~Leakage(){}

//--------------------------------------------------------------
bool Leakage::calculateLeakage()
{
	// calculating how much water is lost through NODAL leakage
	isLeakage = true;
	bool succesfull = solveSystem();
	double summarizedDemand = 0.0;
	for(int i=0; i<numberNodes; i++)
	{
		summarizedLeakage += nodes[i]->getProperty("leakage");
		cout << i << ". node leakage:" << nodes[i]->getProperty("leakage") << endl;
		summarizedDemand += nodes[i]->getProperty("demand");
		cout << i << ". node demand:" << nodes[i]->getProperty("demand") << endl;
	}
	cout << "summarized leakage: " << summarizedLeakage << endl;
	cout << "summarized demand: " << summarizedDemand << endl;
	return succesfull;
}

double Leakage::calculateUnservedDemands()
{
	double sumd=0.; // in m3/s
	for(int i=0; i<numberNodes; i++)
	{
		sumd += nodes[i]->demand;
	}
	isPressureDemand = true;
	bool successfull = solveSystem();

	double b = 0., bi = 0.;
	for(int j=0; j<numberNodes; j++)
	{
		b += nodes[j]->getProperty("consumption");
	}
	bi = (sumd - b)/sumd;
	if (successfull == false)
	{
		cout << "Convergence error (PD SOLVER)" << endl << "Overwrite with 10, the original value is: " << bi << endl;
		bi = 10.;
	}
	return bi;
}

void Leakage::setLeakageConstant(double C)
{
	for(int i=0; i<numberNodes; i++)
	{
		cout << "Leakage constant modified from: " << nodes[i]->getProperty("leakageConstant") << " to: ";
		nodes[i]->setProperty("leakageConstant",C);
		cout << nodes[i]->getProperty("leakageConstant")  << " , " << C << endl;
	}
}

/*void Leakage::setLeakageConstant(vector<int> ID, vector<double> C)
{
	for (int i = 0; i < ID.size(); ++i)
	{
		for(int j=0; j<numberNodes; j++)
		{
			if (ID[i] == j)
			{
				cout << "Leakage constant modified from: " << nodes[j]->getProperty("leakageConstant") << " to: ";
				nodes[j]->setProperty("leakageConstant",C);
				cout << nodes[j]->getProperty("leakageConstant") << endl;
			}
		}
	}

}*/