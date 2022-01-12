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

#include "SeriesHydraulics.h"

using namespace std;

class BIWS
{
public:
	BIWS(string fileName);
	~BIWS();

	// main STACI case
	SeriesHydraulics *wds;

	// calculating fitness functions: ff and I
	void evaluate();

	// modifying stuff
	double leakageRepair(int year, string leakageID);
	double replacePipe(int year, string pipeID, double newDiameter);
	double increaseTank(int year, string tankID, double newVolume);
	double installValve(int year, string pipeID, bool isStart, string type);
	double replacePump(int year, string pumpID, double newQ, double newH);
	double installFrequencyInverter(int year, string pumpID);

	// fitnes functions
	vector<double> I;
	vector<vector<double> > ff;

	// printing fitness functions to console
	void printFitnessFunction();

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
	
	// leakge constants
	double leakageExponent = 1.0;
	double leakageMinPressure = 0.0;
	vector<string> leakageEdgeID;
	vector<double> leakageEdgeLength;
	vector<double> leakageCoefficient;
	vector<double> k0; // for saving initial leakage const
};

#endif
