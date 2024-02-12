#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../Vulnerability.h"

using namespace std;
using namespace Eigen;

double standardizeDiameter(const vector<double>& st_diam, double diam);
void standardizeDiameter(const vector<double>& st_diam, vector<double> &diam);

int main(int argc, char* argv[])
{
	// Name of containing folder of staci file
	string caseFolder = "../../Networks/Sopron/";

	// standard diameter values
	vector<double> st_diam{50.,63.,75.,100.,125.,150.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.}; 

	string caseName, runType;
	if(argc < 3)
	{
		cout << "Wrong number of inputs " << argc << ", correct: 2/3 (caseName, runType[pn,sn,of,of2,of_p1,of_p2], standardize(sn)/pref(of_p1,of_p2))" << endl;
		exit(-1);
	}
	else
	{
		caseName = argv[1];
		runType = argv[2];
	}
	
	Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");
	if(runType == "pn") // pipe number
	{
		int pn = wds->pipeIndex.size();
		ofstream wFile;
		wFile.open("pn.txt");
		wFile << pn << endl;
		wFile.close();
	}
	else if(runType == "sn") // just saving the network with diameter.txt values
	{
		vector<double> diameter = readVectorDouble("diameter.txt");

		bool do_st;
		istringstream("1") >> do_st;
		if(do_st)
		{
			standardizeDiameter(st_diam, diameter);
		}

		for(int i=0; i<wds->pipeIndex.size(); i++)
		{  
			wds->edges[wds->pipeIndex[i]]->setDoubleProperty("diameter",diameter[i]*1.e-3);
		}

		wds->saveSystem(caseFolder + caseName + "_opt.inp");
	}
	else if(runType == "of") // classic optimisation
	{
		double cost=0., cost_orig=0.;
		double Arp=13.,Brp=29.,Crp=1200.;
		vector<double> diameter = readVectorDouble("diameter.txt");
		for(int i=0; i<wds->pipeIndex.size(); i++)
		{  
			double length=wds->edges[wds->pipeIndex[i]]->getDoubleProperty("length");

			double diam_orig =  wds->edges[wds->pipeIndex[i]]->getDoubleProperty("diameter");
			cost_orig += (Arp + Brp*diam_orig + Crp*diam_orig*diam_orig)*length;

			wds->edges[wds->pipeIndex[i]]->setDoubleProperty("diameter",diameter[i]*1.e-3);
			cost += (Arp + Brp*diameter[i]*1.e-3 + Crp*diameter[i]*diameter[i]*1.e-6)*length;
		}
		double cost_rel = cost/cost_orig;

		wds->isPressureDemand = true;
		wds->solveSystem();
		double rs=0., sumd=0.;
		for(int i=0; i<wds->nodes.size(); i++)
		{
			double d = wds->nodes[i]->getProperty("demand");
			if(d!=0)
			{
				double c = wds->nodes[i]->getProperty("consumption");
				rs += (d-c);
				sumd += d;
			}
		}
		rs /= sumd;

		double a = 1000.;
		double b = 1.;
		double of = a*rs + b*cost_rel;

		ofstream wFile;
		wFile.open("of.txt");
		wFile << of << endl;
		wFile.close();
	}
	else if(runType == "of2") // backup optimisation
	{
		double cost=0., cost_orig=0.;
		double Arp=13.,Brp=29.,Crp=1200.;
		vector<double> diameter = readVectorDouble("diameter2.txt");
		for(int i=0; i<wds->pipeIndex.size(); i++)
		{  
			double length=wds->edges[wds->pipeIndex[i]]->getDoubleProperty("length");

			double diam_orig =  wds->edges[wds->pipeIndex[i]]->getDoubleProperty("diameter");
			cost_orig += (Arp + Brp*diam_orig + Crp*diam_orig*diam_orig)*length;

			wds->edges[wds->pipeIndex[i]]->setDoubleProperty("diameter",diameter[i]*1.e-3);
			cost += (Arp + Brp*diameter[i]*1.e-3 + Crp*diameter[i]*diameter[i]*1.e-6)*length;
		}
		double cost_rel = cost/cost_orig;

		// cout << "cost: " << cost << endl;
		// cout << "cost_orig: " << cost_orig << endl;

		wds->calculateVulnerability();
		double rs = wds->backupRelativeShortfall;

		// cout << "rs: " << rs << endl;

		double a = 1000.;
		double b = 1.;
		double of2 = a*rs + b*cost_rel;

		//cout << "of2: " << of2 << endl;

		ofstream wFile;
		wFile.open("of2.txt");
		wFile << of2 << endl;
		wFile.close();
	}
	else if(runType == "of_p1" || runType == "of_p2") // objective function for pressure with standard diameters 
	{
		// loading diameters, standardising, calculating cost
		double cost=0., cost_orig=0.;
		double Arp=13.,Brp=29.,Crp=1200.;
		vector<double> diameter;
		if(runType == "of_p1")
		{
			diameter = readVectorDouble("diameter.txt");
		}
		else
		{
			diameter = readVectorDouble("diameter2.txt");
		}
		
		// standardised diameter values
		standardizeDiameter(st_diam,diameter);

		for(int i=0; i<wds->pipeIndex.size(); i++)
		{  
			double length=wds->edges[wds->pipeIndex[i]]->getDoubleProperty("length");

			double diam_orig =  wds->edges[wds->pipeIndex[i]]->getDoubleProperty("diameter");
			cost_orig += (Arp + Brp*diam_orig + Crp*diam_orig*diam_orig)*length;

			wds->edges[wds->pipeIndex[i]]->setDoubleProperty("diameter",diameter[i]*1.e-3);
			cost += (Arp + Brp*diameter[i]*1.e-3 + Crp*diameter[i]*diameter[i]*1.e-6)*length;
		}
		double cost_rel = cost/cost_orig;

		if(runType == "of_p1" || runType == "of_p2") // for the single network
		{
			wds->isPressureDemand = false;
			wds->solveSystem();

			double pres_rel=0.;
			double pref = stod(argv[3],0);
			for(int i=0; i<wds->nodes.size(); i++)
			{
				pres_rel += abs(wds->nodes[i]->head - pref);
			}
			pres_rel /= pref;
			pres_rel /= wds->numberNodes;
			double of_p1 = 1./pres_rel;

			wds->presRef = stod(argv[3],0);
			wds->calculateVulnerability();
			double of_p2 = 1./wds->backupPressure;

			// writing the objective function value
			if(runType == "of_p1")
			{
				ofstream wFile;
				wFile.open("of_p1.txt");
				wFile << cost_rel << endl;
				wFile << of_p1 << endl;
				wFile.close();
			}
			else
			{
				ofstream wFile;
				wFile.open("of_p2.txt");
				wFile << cost_rel << endl;
				wFile << of_p2 << endl;
				wFile.close();
			}

			// writing every other value to log file
			FILE *wfile;
			wfile = fopen((caseName+"_log.txt").c_str(),"a");

			fprintf(wfile, "%8.5e,%8.5e,%8.5e,%8.5e,",cost_rel,of_p1,of_p2,pref);
			for(int i=0; i<wds->pipeIndex.size(); i++)
			{
			 	fprintf(wfile, "%8.5e,",wds->edges[wds->pipeIndex[i]]->getDoubleProperty("diameter"));
			}
			fprintf(wfile,"\n");
  			fclose(wfile);
		}
	}
	
	return 0;
}

double standardizeDiameter(const vector<double>& st_diam, double diam)
{
	double new_diam=diam;
	for(int i=1; i<st_diam.size(); i++)
	{
		if (st_diam[i-1] < diam && st_diam[i] >= diam)
		{
			new_diam = st_diam[i];
		}
	}
	return new_diam;
}

void standardizeDiameter(const vector<double>& st_diam, vector<double> &diam)
{
	for(int i=0; i<diam.size(); i++)
	{
		diam[i] = standardizeDiameter(st_diam, diam[i]);
	}
}