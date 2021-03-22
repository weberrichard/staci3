#include "Leakage.h"

//--------------------------------------------------------------
Leakage::Leakage(string spr_filename) : HydraulicSolver(spr_filename){}
Leakage::~Leakage(){}

//--------------------------------------------------------------
void Leakage::calculateLeakage()
{
	// calculating how much water is lost through NODAL leakage
	isLeakage = true;
	solveSystem();
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
}

double Leakage::calculateUnservedDemands()
{
	double sumd=0.; // in m3/s
	for(int i=0; i<numberNodes; i++)
	{
		sumd += nodes[i]->demand;
	}
	isPressureDemand = true;
	solveSystem();

	double b = 0., bi = 0.;
	for(int j=0; j<numberNodes; j++)
	{
		b += nodes[j]->getProperty("consumption");
	}
	bi = (sumd - b)/sumd;
	return bi;
}