#include "BIWS.h"

//--------------------------------------------------------------
BIWS::BIWS(string fileName)
{
	// setting fitness functions sizes
	I.resize(9,0.);
	ff.resize(9, vector<double>(6,0.));

	// loading leakage data
	leakageEdgeID.resize(nYear);
	leakageEdgeLength.resize(nYear);
	leakageCoefficient.resize(nYear);
	k0.resize(nYear);
	leakageEdgeID[0] = readVectorString("Leakages_edge_id.csv");
	leakageEdgeLength[0] = readVectorDouble("Leakages_edge_length.csv");
	leakageCoefficient[0] = readVectorDouble("Leakages_coefficient.csv");
	for(int i=1; i<nYear; i++)
	{
		leakageEdgeID[i] = leakageEdgeID[0];
		leakageEdgeLength[i] = leakageEdgeLength[0];
		leakageCoefficient[i] = leakageCoefficient[0];
	}

	// resizing STACI vector
	wds.resize(nYear);

	for(int i=0; i<nYear; i++)
	{
		// loading the INP file
		wds[i] = new SeriesHydraulics(fileName);

		// hydraulic model setting
		wds[i]->isPressureDemand = true;
		wds[i]->isLeakage = true;
		for(int j=0; j<wds[i]->nodes.size(); j++)
		{
			wds[i]->nodes[j]->pdExponent = pdExponent;
			wds[i]->nodes[j]->pdDesiredPressure = pdDesiredPressure;
			wds[i]->nodes[j]->pdMinPressure = pdMinPressure;
			
			wds[i]->nodes[j]->leakageExponent = leakageExponent;
			wds[i]->nodes[j]->leakageMinPressure = leakageMinPressure;
			wds[i]->nodes[j]->leakageConstant = 0.0; // zero for now
		}

		// counting the number of nodes with demands
		nDem=0.;
		for(int j=0; j<wds[i]->nodes.size(); j++)
		{
			if(wds[i]->nodes[j]->demand > 0.)
			{
				nDem += 1.;
			}
		}

		// setting original nominal efficiency
		for(int j=0; j<wds[i]->pumpIndex.size(); j++)
		{
			wds[i]->edges[wds[i]->pumpIndex[j]]->setDoubleProperty("efficiencyNominal",oldEta);
		}
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
		cout << " [*] year: " << i << endl;

		// resizing leakage constant vector
		k0[i].clear();
		k0[i].resize(wds[i]->nodes.size(),.0);

		// setting leakage model coefficients
		for(int j=0; j<leakageEdgeID[i].size(); j++)
		{
			int edge_idx = wds[i]->edgeIDtoIndex(leakageEdgeID[i][j]);
			int sn_idx = wds[i]->edges[edge_idx]->startNodeIndex;
			int en_idx = wds[i]->edges[edge_idx]->endNodeIndex;
			double l = wds[i]->edges[edge_idx]->getDoubleProperty("length");

			k0[i][sn_idx] += (l-leakageEdgeLength[i][j])/l*leakageCoefficient[i][j];
			k0[i][en_idx] += leakageEdgeLength[i][j]/l*leakageCoefficient[i][j];
		}

		// leakage growth over years
		for(int j=0; j<wds[i]->nodes.size(); j++)
		{
			wds[i]->nodes[j]->leakageConstant = k0[i][j]*exp(0.25*(double)i*52./260.);
		}

		// hydraulic solver
		wds[i]->seriesSolve();

		// fitness function
		// I1 - service effectiveness
		// I2 - service continouity
		// I6 - full supply continouity
		double nij=0.;
		for(int j=0; j<wds[i]->nodes.size(); j++)
		{
			double wij=1.;
			double dij=1.;
			for(int k=0; k<wds[i]->vectorTime.size(); k++)
			{
				if((int)wds[i]->vectorTime[k]%3600<1e-2)
				{
					if(wds[i]->nodes[j]->vectorHead[k] > 0.0 && wds[i]->nodes[j]->demand > 0.)
					{
						nij += 52.; // as there are 52 weeks
					}
					if(wds[i]->nodes[j]->vectorHead[k] < 0. && wds[i]->nodes[j]->demand > 0.)
					{
						wij = 0.;
					}
					if(wds[i]->nodes[j]->vectorHead[k] < 10. && wds[i]->nodes[j]->demand > 0.)
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
		for(int j=0; j<wds[i]->presIndex.size(); j++)
		{
			for(int k=0; k<wds[i]->vectorTime.size(); k++)
			{
				if((int)wds[i]->vectorTime[k]%3600 < 1e-2)
				{
					if(wds[i]->edges[wds[i]->presIndex[j]]->vectorVolumeFlowRate[k] < 0.)
					{
						Vs[i] += wds[i]->edges[wds[i]->presIndex[j]]->vectorVolumeFlowRate[k];
					}
				}
			}
		}
		for(int j=0; j<wds[i]->poolIndex.size(); j++)
		{
			for(int k=0; k<wds[i]->vectorTime.size(); k++)
			{
				if((int)wds[i]->vectorTime[k]%3600 < 1e-2)
				{
					if(wds[i]->edges[wds[i]->poolIndex[j]]->vectorVolumeFlowRate[k] < 0.)
					{
						Vs[i] += wds[i]->edges[wds[i]->poolIndex[j]]->vectorVolumeFlowRate[k];
					}
				}
			}
		}
		Vs[i] *= -1.;

		for(int j=0; j<wds[i]->nodes.size(); j++)
		{
			double vd=0., vc=0.; // for I9
			for(int k=0; k<wds[i]->vectorTime.size(); k++)
			{
				if((int)wds[i]->vectorTime[k]%3600 < 1e-2)
				{
					Vl[i] += wds[i]->nodes[j]->vectorLeakage[k];
					Vc[i] += wds[i]->nodes[j]->vectorConsumption[k];
					Vd[i] += wds[i]->nodes[j]->vectorDemand[k];
					if(wds[i]->nodes[j]->demand > 0.)
					{
						ff[4][i] += min(wds[i]->nodes[j]->vectorHead[k],pRef);
						vc += wds[i]->nodes[j]->vectorConsumption[k];
						vd += wds[i]->nodes[j]->vectorDemand[k];
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
		vector<double> Lmax(wds[i]->edges.size(),.0);
		for(int j=0; j<wds[i]->nodes.size(); j++)
		{
			for(int k=0; k<wds[i]->vectorTime.size(); k++)
			{
				if((int)wds[i]->vectorTime[k]%3600 < 1e-2)
				{
					if(wds[i]->nodes[j]->vectorHead[k] < 0.)
					{
						for(int l=0; l<wds[i]->nodes[j]->edgeIn.size(); l++)
						{
							int edge_idx = wds[i]->nodes[j]->edgeIn[l];
							int tc = wds[i]->edges[edge_idx]->typeCode;
							if(tc == 0 || tc == 1)
							{
								int si = wds[i]->edges[edge_idx]->startNodeIndex;
								double length = wds[i]->edges[edge_idx]->getDoubleProperty("length");
								if(wds[i]->nodes[si]->vectorHead[k] < 0.)
								{
									if(Lmax[edge_idx]<length)
										Lmax[edge_idx] = length;
								}
								else
								{
									double L = -wds[i]->nodes[j]->vectorHead[k]*length/(wds[i]->nodes[si]->vectorHead[k]-wds[i]->nodes[j]->vectorHead[k]);
									if(Lmax[edge_idx]<L)
										Lmax[edge_idx] = L;
								}
							}
						}
						for(int l=0; l<wds[i]->nodes[j]->edgeOut.size(); l++)
						{
							int edge_idx = wds[i]->nodes[j]->edgeOut[l];
							int tc = wds[i]->edges[edge_idx]->typeCode;
							if(tc == 0 || tc == 1)
							{
								int ei = wds[i]->edges[edge_idx]->endNodeIndex;
								double length = wds[i]->edges[edge_idx]->getDoubleProperty("length");
								if(wds[i]->nodes[ei]->vectorHead[k] < 0.)
								{
									if(Lmax[edge_idx]<length)
										Lmax[edge_idx] = length;
								}
								else
								{
									double L = -wds[i]->nodes[j]->vectorHead[k]*length/(wds[i]->nodes[ei]->vectorHead[k]-wds[i]->nodes[j]->vectorHead[k]);
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
		for(int j=0; j<wds[i]->pumpIndex.size(); j++)
		{
			int si = wds[i]->edges[wds[i]->pumpIndex[j]]->startNodeIndex;
			int ei = wds[i]->edges[wds[i]->pumpIndex[j]]->endNodeIndex;
			vector<double> qCurve = wds[i]->edges[wds[i]->pumpIndex[j]]->getVectorProperty("qCurve");
			double eta0 = wds[i]->edges[wds[i]->pumpIndex[j]]->getDoubleProperty("efficiencyNominal");
			for(int k=0; k<wds[i]->vectorTime.size(); k++)
			{
				if((int)wds[i]->vectorTime[k]%3600 < 1e-2)
				{		
					double dp = wds[i]->nodes[ei]->vectorHead[k]-wds[i]->nodes[si]->vectorHead[k];
					double q = wds[i]->edges[wds[i]->pumpIndex[j]]->vectorVolumeFlowRate[k];
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
