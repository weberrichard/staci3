/*===================================================================*\
                                  BIWS
                            ---------------

  Battle of intermittent water supply challnge
  
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef BIWS_H
#define BIWS_H

#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <iostream>
#include <vector>
#include "SeriesHydraulics.h"

using namespace std;

class BIWS
{
public:
	BIWS(string fileName);
	virtual ~BIWS();

	// main STACI case
	vector<SeriesHydraulics*> wds;

	// calculating fitness functions: ff and I
	void evaluate();
	void readCostTable(string fname, string whichCost);
	// modifying stuff
	double leakageRepair(int year, string leakageID_local, double leakageEdgeLength_local, bool Activate);
	double replacePipe(int year, string pipeID, double newDiameter, bool Activate);
	double increaseTank(int year, string tankID, double newVolume, bool Activate);
	double installValve(int year, string pipeID, bool isStart, string type, double setting, bool Activate);
	double replacePump(int year, string pumpID, double newQ, double newH, bool Activate);
	double installFrequencyInverter(int year, string pumpID, double RevRate, bool Activate);

	// fitnes functions
	vector<double> I;
	vector<vector<double> > ff;

	bool PipeReplacementCost_Read = false;//->Menjen privatebe
	bool ValvePlacementCost_Read = false;

	vector<vector<string> > PipeReplacementCost;
	vector<vector<string> > ValvePlacementCost;

	// printing fitness functions to console
	void printFitnessFunction();
	void undo_LeakageRepair(int year, string leakageID_local, double leakageEdgeLength_local, double leakageCoefficient_local, bool Activate);
	void undo_PipeReplace(int year, string leakageID_local);
	vector<vector<string> > leakageEdgeID, failsafe_leakageEdgeID;
	vector<vector<double> > leakageEdgeLength, failsafe_leakageEdgeLength;
	vector<vector<double> > leakageCoefficient, failsafe_leakageCoefficient;

private:
	// constants of biws
	// number of nodes with demand
	double nDem;

	// number of years
	int nYear = 6;

	// original and new nominal efficiency of pumps
	double oldEta = 0.65;
	double newEta = 0.8;

	// roughness of new pipeline
	double newHW = 120.;

	// pressure-dependent demands
	double pRef = 20.;
	double pdExponent = 2.0;
	double pdDesiredPressure = 10.;
	double pdMinPressure = 0.0;
	double roughness_local = 0.;
	// leakge constants
	double leakageExponent = 1.0;
	double leakageMinPressure = 0.0;
	vector<vector<double> > k0; // for saving initial leakage const
};

#endif
