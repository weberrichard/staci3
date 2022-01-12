#include "BIWS.h"

//--------------------------------------------------------------
BIWS::BIWS(string fileName)
{
	// setting fitness functions sizes
	I.resize(9,0.);
	ff.resize(9, vector<double>(6,0.));

	// loading leakage data
	leakageEdgeID = readVectorString("Leakages_edge_id.csv");
	leakageEdgeLength = readVectorDouble("Leakages_edge_length.csv");
	leakageCoefficient = readVectorDouble("Leakages_coefficient.csv");

	// loading the INP file
	wds = new SeriesHydraulics(fileName);

	// hydraulic model setting
	wds->isPressureDemand = true;
	wds->isLeakage = true;
	for(int i=0; i<wds->nodes.size(); i++)
	{
		wds->nodes[i]->pdExponent = pdExponent;
		wds->nodes[i]->pdDesiredPressure = pdDesiredPressure;
		wds->nodes[i]->pdMinPressure = pdMinPressure;
		
		wds->nodes[i]->leakageExponent = leakageExponent;
		wds->nodes[i]->leakageMinPressure = leakageMinPressure;
		wds->nodes[i]->leakageConstant = 0.0; // zero for now
	}

	// setting leakage model coefficients
	for(int j=0; j<leakageEdgeID.size(); j++)
	{
		int edge_idx = wds->edgeIDtoIndex(leakageEdgeID[j]);
		int sn_idx = wds->edges[edge_idx]->startNodeIndex;
		int en_idx = wds->edges[edge_idx]->endNodeIndex;
		double l = wds->edges[edge_idx]->getDoubleProperty("length");

		wds->nodes[sn_idx]->leakageConstant += (l-leakageEdgeLength[j])/l*leakageCoefficient[j];
		wds->nodes[en_idx]->leakageConstant += leakageEdgeLength[j]/l*leakageCoefficient[j];
	}

	// saving initial leakage coefficients
	k0.resize(wds->nodes.size(),0.);
	for(int j=0; j<wds->nodes.size(); j++)
	{
		k0[j] = wds->nodes[j]->leakageConstant;
	}

	// counting the number of nodes with demands
	nDem=0.;
	for(int j=0; j<wds->nodes.size(); j++)
	{
		if(wds->nodes[j]->demand > 0.)
		{
			nDem += 1.;
		}
	}

	// setting original nominal efficiency
	for(int j=0; j<wds->pumpIndex.size(); j++)
	{
		wds->edges[wds->pumpIndex[j]]->setDoubleProperty("efficiencyNominal",oldEta);
	}
}

//--------------------------------------------------------------
BIWS::~BIWS(){}

//--------------------------------------------------------------
void BIWS::evaluate()
{
	// for partial results of I and ff
	vector<double> Vl(nYear,0.), Vs(nYear,0.); // I3
	vector<double> Vc(nYear,0.), Vd(nYear,0.); // I4
	vector<vector<double> > SR(nYear); // I9

	// start of main cycle
	for(int i=0; i<nYear; i++)
	{
		// leakage growth over years
		if(i > 0)
		{
			for(int j=0; j<wds->nodes.size(); j++)
			{
				wds->nodes[j]->leakageConstant = k0[j]*exp(0.25*(double)i*52./260.);
			}
		}

		// hydraulic solver
		wds->seriesSolve();

		// fitness function
		// I1 - service effectiveness
		// I2 - service continouity
		// I6 - full supply continouity
		double nij=0.;
		for(int j=0; j<wds->nodes.size(); j++)
		{
			double wij=1.;
			double dij=1.;
			for(int k=0; k<wds->vectorTime.size(); k++)
			{
				if((int)wds->vectorTime[k]%3600<1e-2)
				{
					if(wds->nodes[j]->vectorHead[k] > 0.0 && wds->nodes[j]->demand > 0.)
					{
						nij += 52.; // as there are 52 weeks
					}
					if(wds->nodes[j]->vectorHead[k] < 0. && wds->nodes[j]->demand > 0.)
					{
						wij = 0.;
					}
					if(wds->nodes[j]->vectorHead[k] < 10. && wds->nodes[j]->demand > 0.)
					{
						dij = 0.;
					}
				}
			}
			ff[1][i] += wij/nDem;
			ff[5][i] += dij/nDem;
		}
		ff[0][i] += nij/(nDem*24.*364.);
		I[0] += ff[0][i]/nYear;
		I[1] += ff[1][i]/nYear; // todo 6 or 5 ?
		I[5] += ff[5][i]/nYear;

		// I3 - leak vs. input water
		// I4 - consumption vs. demand
		// I5 - pressure
		// I9 - supply equity, calculating and saving SR values
		for(int j=0; j<wds->presIndex.size(); j++)
		{
			for(int k=0; k<wds->vectorTime.size(); k++)
			{
				if((int)wds->vectorTime[k]%3600 < 1e-2)
				{
					if(wds->edges[wds->presIndex[j]]->vectorVolumeFlowRate[k] < 0.)
					{
						Vs[i] += wds->edges[wds->presIndex[j]]->vectorVolumeFlowRate[k];
					}
				}
			}
		}
		for(int j=0; j<wds->poolIndex.size(); j++)
		{
			for(int k=0; k<wds->vectorTime.size(); k++)
			{
				if((int)wds->vectorTime[k]%3600 < 1e-2)
				{
					if(wds->edges[wds->poolIndex[j]]->vectorVolumeFlowRate[k] < 0.)
					{
						Vs[i] += wds->edges[wds->poolIndex[j]]->vectorVolumeFlowRate[k];
					}
				}
			}
		}
		Vs[i] *= -1.;

		for(int j=0; j<wds->nodes.size(); j++)
		{
			double vd=0., vc=0.; // for I9
			for(int k=0; k<wds->vectorTime.size(); k++)
			{
				if((int)wds->vectorTime[k]%3600 < 1e-2)
				{
					Vl[i] += wds->nodes[j]->vectorLeakage[k];
					Vc[i] += wds->nodes[j]->vectorConsumption[k];
					Vd[i] += wds->nodes[j]->vectorDemand[k];
					if(wds->nodes[j]->demand > 0.)
					{
						ff[4][i] += min(wds->nodes[j]->vectorHead[k],pRef);
						vc += wds->nodes[j]->vectorConsumption[k];
						vd += wds->nodes[j]->vectorDemand[k];
					}
				}
			}
			if(vd > 0.)
				SR[i].push_back(vc/vd);
		}
		ff[2][i] = Vl[i]/Vs[i];
		ff[3][i] = Vc[i]/Vd[i];
		ff[4][i] /= pRef*168.*nDem;
		I[4] += ff[4][i]/nYear;

		// I7 - negative pressures
		vector<double> Lmax(wds->edges.size(),.0);
		for(int j=0; j<wds->nodes.size(); j++)
		{
			for(int k=0; k<wds->vectorTime.size(); k++)
			{
				if((int)wds->vectorTime[k]%3600 < 1e-2)
				{
					if(wds->nodes[j]->vectorHead[k] < 0.)
					{
						for(int l=0; l<wds->nodes[j]->edgeIn.size(); l++)
						{
							int edge_idx = wds->nodes[j]->edgeIn[l];
							int tc = wds->edges[edge_idx]->typeCode;
							if(tc == 0 || tc == 1)
							{
								int si = wds->edges[edge_idx]->startNodeIndex;
								double length = wds->edges[edge_idx]->getDoubleProperty("length");
								if(wds->nodes[si]->vectorHead[k] < 0.)
								{
									if(Lmax[edge_idx]<length)
										Lmax[edge_idx] = length;
								}
								else
								{
									double L = -wds->nodes[j]->vectorHead[k]*length/(wds->nodes[si]->vectorHead[k]-wds->nodes[j]->vectorHead[k]);
									if(Lmax[edge_idx]<L)
										Lmax[edge_idx] = L;
								}
							}
						}
						for(int l=0; l<wds->nodes[j]->edgeOut.size(); l++)
						{
							int edge_idx = wds->nodes[j]->edgeOut[l];
							int tc = wds->edges[edge_idx]->typeCode;
							if(tc == 0 || tc == 1)
							{
								int ei = wds->edges[edge_idx]->endNodeIndex;
								double length = wds->edges[edge_idx]->getDoubleProperty("length");
								if(wds->nodes[ei]->vectorHead[k] < 0.)
								{
									if(Lmax[edge_idx]<length)
										Lmax[edge_idx] = length;
								}
								else
								{
									double L = -wds->nodes[j]->vectorHead[k]*length/(wds->nodes[ei]->vectorHead[k]-wds->nodes[j]->vectorHead[k]);
									if(Lmax[edge_idx]<L)
										Lmax[edge_idx] = L;
								}
							}
						}
					}
				}
			}
		}
		double lmax=0.;
		for(int j=0; j<Lmax.size(); j++)
		{
			lmax += Lmax[j];
		}
		ff[6][i] = lmax;
		I[6] += lmax/6.;

		// I8 - energy consumption of the pumps
		for(int j=0; j<wds->pumpIndex.size(); j++)
		{
			int si = wds->edges[wds->pumpIndex[j]]->startNodeIndex;
			int ei = wds->edges[wds->pumpIndex[j]]->endNodeIndex;
			vector<double> qCurve = wds->edges[wds->pumpIndex[j]]->getVectorProperty("qCurve");
			double eta0 = wds->edges[wds->pumpIndex[j]]->getDoubleProperty("efficiencyNominal");
			for(int k=0; k<wds->vectorTime.size(); k++)
			{
				if((int)wds->vectorTime[k]%3600 < 1e-2)
				{		
					double dp = wds->nodes[ei]->vectorHead[k]-wds->nodes[si]->vectorHead[k];
					double q = wds->edges[wds->pumpIndex[j]]->vectorVolumeFlowRate[k];
					double eta = eta0*(2.*q/qCurve[0]-pow(q/qCurve[0],2));
					ff[7][i] += dp*q*9.81/eta; // kWh is needed, so 1000 [kg/m3] / 1000 [W/kW]
				}
			}
		}
		I[7] += ff[7][i];
	}

	// I3 calculation
	double Vss=0., Vls=0.;
	for(int i=0; i<Vs.size(); i++)
	{
		Vss += Vs[i];
		Vls += Vl[i];
	}
	I[2] = Vls/Vss;

	// I4 calculation
	double Vcs=0., Vds=0.;
	for(int i=0; i<Vc.size(); i++)
	{
		Vcs += Vc[i];
		Vds += Vd[i];
	} 
	I[3] = Vcs/Vds;

	// I9 calculation
	double ASR = 0.;
	vector<double> vASR(6,0.);
	for(int i=0; i<SR.size(); i++)
	{
		for(int j=0; j<SR[i].size(); j++)
		{
			ASR += SR[i][j];
			vASR[i] += SR[i][j];
		}
		vASR[i] /= nDem;
	}
	ASR /= 6.*nDem;
	double ADEV = 0.;
	vector<double> vADEV(6,0.);
	for(int i=0; i<SR.size(); i++)
	{
		for(int j=0; j<SR[i].size(); j++)
		{
			ADEV += abs(SR[i][j]-ASR);
			vADEV[i] += abs(SR[i][j]-vASR[i]);
		}
		vADEV[i] /= nDem;
		ff[8][i] = 1-vADEV[i]/vASR[i];
	}
	ADEV /= 6.*nDem;
	I[8] = 1-ADEV/ASR;

}

//--------------------------------------------------------------
double BIWS::leakageRepair(int year, string leakageID)
{
	double cost=0.;
	return cost;
}

//--------------------------------------------------------------
double BIWS::replacePipe(int year, string pipeID, double newDiameter)
{
	// use newHW variable for Hazen coef
	double cost=0.;
	return cost;
}

//--------------------------------------------------------------
double BIWS::increaseTank(int year, string tankID, double newVolume)
{

	double cost=0.;
	return cost;
}

//--------------------------------------------------------------
double BIWS::installValve(int year, string pipeID, bool isStart, string type)
{
	double cost=0.;
	return cost;
}

//--------------------------------------------------------------
double BIWS::replacePump(int year, string pumpID, double newQ, double newH)
{
	double cost=0.;
	return cost;
}

//--------------------------------------------------------------
double BIWS::installFrequencyInverter(int year, string pumpID)
{
	double cost=0.;
	return cost;
}

//--------------------------------------------------------------
void BIWS::printFitnessFunction()
{
	// printing fitness function to console
	printf("------------------------------------------------------------------------------------------\n");
	printf(" year |     0     |     1     |     2     |     3     |     4     |     5     ||     I     |\n");
	printf("------------------------------------------------------------------------------------------\n");
	for(int i=0; i<ff.size(); i++)
	{
		printf(" ff%1i  |",i+1);
		for(int j=0; j<ff[i].size(); j++)
		{
			printf("%10.3f |",ff[i][j]);
		}
		printf("|%10.3f |",I[i]);
		printf("\n");
	}
	printf("------------------------------------------------------------------------------------------\n");
}
