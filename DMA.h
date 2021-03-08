/*===================================================================*\
                                  DMA
                            ---------------
	
	This class is capable of calculating the DMA zones.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef DMA_H
#define DMA_H

#include "HydraulicSolver.h"

class DMA : public HydraulicSolver
{
public:
	DMA(string spr_filename);
	~DMA();

	void determineDMAZones();
	int numberDMAZones; // number of DMA zones

	void addNewFlowMeter(vector<string> flowMeterName, vector<string> pipeName, vector<bool> isStart, double density, vector<double> referenceCrossSection, double volumeFlowRate);
   void addNewFlowMeter(vector<string> flowMeterName, vector<int> pipeIndex, vector<bool> isStart, double density, vector<double> referenceCrossSection, double volumeFlowRate);

private:
	vector<int> edgeVector; // edge vector of the original network without flow meters
	vector<vector<int> > zoneVector; // contains the original node indicies of pipes
	vector<int> zoneEdgeVector; // edge vector of the DMA zones
	void updateEdgeVector();
};

#endif