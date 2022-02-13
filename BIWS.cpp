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

	// ----- Saving the original tables ------//
	failsafe_leakageEdgeID = leakageEdgeID;
	failsafe_leakageEdgeLength = leakageEdgeLength;
	failsafe_leakageCoefficient = leakageCoefficient;
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

		///----------- Read cost tables ---------------///
		readCostTable("PipeReplacementCost.csv","PipeReplacementCost");
		readCostTable("ValvePlacementCost.csv","ValvePlacementCost");
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

	I.clear();
	ff.clear();
	I.resize(9,0.);
	ff.resize(9, vector<double>(6,0.));
	// start of main cycle
	for(int i=0; i<nYear; i++)
	{
		//cout << " [*] year: " << i << endl;

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
double BIWS::leakageRepair(int year, string leakageID_local, double leakageEdgeLength_local, bool Activate)
{
	double cost=0., leakageCoefficientLocal = 0., pipeDiameter = 0., C_det = 0., C_repair = 0.;
	int EdgeID;
	for (int i = 0; i < leakageEdgeID.at(year).size(); ++i)
	{
		if(leakageEdgeID[year].at(i) == leakageID_local && leakageEdgeLength_local == leakageEdgeLength[year].at(i))
		{
			leakageCoefficientLocal = leakageCoefficient[year][i];
			EdgeID = wds[year]->edgeIDtoIndex(leakageID_local);
			pipeDiameter = wds[year]->edges.at(EdgeID)->getDoubleProperty("diameter");
			//cout << pipeDiameter << " " << EdgeID << " " << leakageCoefficientLocal << endl;
			if(Activate == true)
			{
				for (int j = year; j < nYear; ++j)
				{
					leakageEdgeLength.at(j).erase(leakageEdgeLength[j].begin()+i);
					leakageCoefficient.at(j).erase(leakageCoefficient[j].begin()+i);
					leakageEdgeID.at(j).erase(leakageEdgeID.at(j).begin()+i);
				}
			}
		}
	}
	C_repair = (94-0.3*pipeDiameter+0.01*pipeDiameter*pipeDiameter)*(1.5+0.11*log10(leakageCoefficientLocal));
	//cout << " C_rep: " << C_repair << endl;
	C_det = 2400*exp(-28*leakageCoefficientLocal);
	//cout << " C_det: " << C_det << endl;
	cost = C_det+C_repair;
	//cout << " a csere költsége: " << cost << endl;
	if(cost == 0.)
		cout << "ZERO COST DETECTED!!!!!!!!!!!!! " << year << " ID: " << leakageID_local << " leakageRepair" << endl;
	return cost;
}

//--------------------------------------------------------------
void BIWS::undo_LeakageRepair(int year, string leakageID_local, double leakageEdgeLength_local, double leakageCoefficient_local, bool Activate)
{
	double cost=0., leakageCoefficientLocal = 0., pipeDiameter = 0., C_det = 0., C_repair = 0.;
	for (int j = year; j < nYear; ++j)
	{
		leakageEdgeLength[j].insert(end(leakageEdgeLength[j]),leakageEdgeLength_local);
		leakageCoefficient[j].insert(end(leakageCoefficient[j]),leakageCoefficient_local);
		leakageEdgeID[j].insert(end(leakageEdgeID[j]),leakageID_local);
	}
}

//--------------------------------------------------------------
double BIWS::replacePipe(int year, string pipeID, double newDiameter, bool Activate)
{
	double cost=0., pipeDiameter = 0., pipeLength = 0., CostPerMeter = 0.;
	double Arp = 13., Brp = 29., Crp = 1200.;
	bool IdentifiedDiameter = false, found = false, getBack = true;
	int pipeIndex = 0, leakageCoefficientLocal = 0, EdgeID = 0;
	pipeIndex = wds[year]->edgeIDtoIndex(pipeID);
	pipeLength = wds[year]->edges.at(pipeIndex)->getDoubleProperty("length");
	roughness_local = wds[year]->edges.at(pipeIndex)->getDoubleProperty("roughness");
	if (newDiameter == -1)
	{
		newDiameter = wds[year]->edges.at(pipeIndex)->getDoubleProperty("diameter")*1000;
	}
	for (int i = 1; i < PipeReplacementCost[0].size(); ++i)
	{
		if (newDiameter == stod(PipeReplacementCost[0][i])/1000.)
		{
			CostPerMeter = stod(PipeReplacementCost[1][i]);
			IdentifiedDiameter = true;
		}
	}
	if (IdentifiedDiameter == false)
	{
		CostPerMeter = Arp+Brp*(newDiameter/1000)+Crp*(newDiameter/1000)*(newDiameter/1000);
	}
	while (getBack == true)
	{
		getBack = false;
		for (int i = 0; i < leakageEdgeID[year].size(); ++i)
		{
			if(leakageEdgeID[year].at(i) == pipeID)
			{
				found = true;
				leakageCoefficientLocal = leakageCoefficient[year][i];
				pipeDiameter = wds[year]->edges.at(pipeIndex)->getDoubleProperty("diameter");
				if(Activate == true)
				{
					for (int j = year; j < nYear; ++j)
					{
						getBack = true;
						//cout << "leakage removed: " << leakageEdgeID.at(j).at(i) << " | " << leakageEdgeLength.at(j).at(i) << endl;
						leakageEdgeLength.at(j).erase(leakageEdgeLength[j].begin()+i);
						leakageCoefficient.at(j).erase(leakageCoefficient[j].begin()+i);
						leakageEdgeID.at(j).erase(leakageEdgeID.at(j).begin()+i);
						wds[j]->edges.at(pipeIndex)->setDoubleProperty("roughness", 120.);
					}
				}
			}
		}
		if (found == false && Activate == true)
		{
			for (int j = year; j < nYear; ++j)
			{
				wds[j]->edges.at(pipeIndex)->setDoubleProperty("roughness", 120.);
			}
		}
	}
	cost = pipeLength*CostPerMeter;
	cout << "Pipe replace cost: " << cost << endl;
	if(cost == 0.)
		cout << "ZERO COST DETECTED!!!!!!!!!!!!! " << year << " ID: " << pipeID << " PIPE REPLACE" << endl;
	return cost;
}
//--------------------------------------------------------------
void BIWS::undo_PipeReplace(int year, string leakageID_local)
{
	leakageEdgeLength = failsafe_leakageEdgeLength;
	leakageCoefficient = failsafe_leakageCoefficient;
	leakageEdgeID = failsafe_leakageEdgeID;
	for (int j = year; j < nYear; ++j)
	{
		wds[j]->edges.at(wds[j]->edgeIDtoIndex(leakageID_local))->setDoubleProperty("roughness", roughness_local);
	}
}

//--------------------------------------------------------------
double BIWS::increaseTank(int year, string tankID, double newVolume, bool Activate)
{
	double cost=0., oldVolume = 0., area = 0., minLevelLocal = 0., maxLevelLocal = 0., Aref_New = 0., initLevelLocal = 0., initLevel_New = 0.;
	int EdgeID = wds[year]->edgeIDtoIndex(tankID);
	minLevelLocal = wds[year]->edges.at(EdgeID)->getDoubleProperty("minLevel");
	maxLevelLocal = wds[year]->edges.at(EdgeID)->getDoubleProperty("maxLevel");
	initLevelLocal = wds[year]->edges.at(EdgeID)->getDoubleProperty("initLevel");
	area = wds[year]->edges.at(EdgeID)->getDoubleProperty("Aref");
	oldVolume = (maxLevelLocal-minLevelLocal)*area;//Mivel elvileg a max level ugyan annyi marad és hengeres a tartály, ezért
	if (Activate == true)
	{
		Aref_New = newVolume/maxLevelLocal;
		initLevel_New = (area*(initLevelLocal-minLevelLocal))/Aref_New;
		for (int j = year; j < nYear; ++j)
		{
			wds[j]->edges.at(EdgeID)->setDoubleProperty("reference_cross_section",Aref_New);
			wds[j]->edges.at(EdgeID)->setDoubleProperty("initLevel",initLevel_New);
		}
	}
	cout << " newVolume: " << newVolume << " oldVolume: " << oldVolume << endl;
	cost = 2000 + 200*(newVolume-oldVolume);
	cout << " increase Tank cost: " << cost << endl;
	if(cost <= 2000.)
		cout << "ZERO COST DETECTED!!!!!!!!!!!!!" << endl;
	return cost;
}

//--------------------------------------------------------------
double BIWS::installValve(int year, string pipeID, bool isStart, string type, double setting, bool Activate)
{
	double cost=0., edgeDiameter=0., Av = 0., Bv = 0.;
	int EdgeID = wds[year]->edgeIDtoIndex(pipeID);
	edgeDiameter = wds[year]->edges.at(EdgeID)->getDoubleProperty("diameter");
	cout << edgeDiameter << endl;
	for (int i = 0; i < ValvePlacementCost.size(); ++i)
	{
		if (ValvePlacementCost[i][0] == type)
		{
			Av = stod(ValvePlacementCost[i][1]);
			Bv = stod(ValvePlacementCost[i][2]);
		}
	}
	if (Activate == true)
	{
		for (int j = year; j < nYear; ++j)
		{
			cout << "majd prv ide..." << endl;
			//wds[j]->addNewPRVValves("new_Valve", pipeID, 0., 1000., 0.1, setting, 0.1, 0.);
		}
	}
	cost = Av*pow(edgeDiameter,Bv);
	cout << "install Valve cost: " << cost << endl;
	return cost;
}

//--------------------------------------------------------------
double BIWS::replacePump(int year, string pumpID, double newQ, double newH, bool Activate)
{
	double cost=0., costOfPipeExtractionPerMeter = 500., baseLvl = 0., pumpLvl = 0., C_ex = 0., C_np = 0., P_BEP = 0.;//->Kilowattban kell
	int SNI = 0, ENI = 0;
	int EdgeID = wds[year]->edgeIDtoIndex(pumpID);
	SNI = wds[year]->edges.at(EdgeID)->startNodeIndex;
	ENI = wds[year]->edges.at(EdgeID)->endNodeIndex;
	P_BEP = newQ*newH*1000*9.81/1000.;
	for (int i = 0; i < wds[year]->presIndex.size(); ++i)
	{
		if ( (wds[year]->edges.at(wds[year]->poolIndex[i])->startNodeIndex == ENI && wds[year]->edges.at(i)->type == "Reservoir") || (wds[year]->edges.at(i)->endNodeIndex == ENI && wds[year]->edges.at(i)->type == "Reservoir"))//Reservoirindexú vector<int> 
		{
			baseLvl = wds[year]->edges.at(i)->getEdgeDoubleProperty("head");
			pumpLvl = wds[year]->nodes.at(SNI)->getProperty("height"); 
		}
		else if(wds[year]->edges.at(i)->startNodeIndex == SNI && wds[year]->edges.at(i)->type == "Reservoir" || wds[year]->edges.at(i)->endNodeIndex == SNI && wds[year]->edges.at(i)->type == "Reservoir")
		{
			baseLvl = wds[year]->edges.at(i)->getEdgeDoubleProperty("head");
			pumpLvl = wds[year]->nodes.at(ENI)->getProperty("height"); 
		}
	}
	if(baseLvl == 0. || pumpLvl == 0.)
		cout << " No reservoir found.... !!! " << pumpID << " replacePump" << endl;
	if (Activate == true)
	{
		for (int j = year; j < nYear; ++j)
		{
			wds[j]->edges.at(EdgeID)->setDoubleProperty("volumeFlowRate",newQ);
			wds[j]->edges.at(EdgeID)->setDoubleProperty("efficiencyNominal",0.8);
		}
	}
	C_ex = (pumpLvl - baseLvl)*costOfPipeExtractionPerMeter;
	C_np = 1475*pow(P_BEP,0.525);
	cost = C_ex+C_np;
	if(cost <= 0.)
		cout << "ZERO COST DETECTED!!!!!!!!!!!!!" << endl;
	cout << " cost: " << cost << endl;
	return cost;
}

//--------------------------------------------------------------
double BIWS::installFrequencyInverter(int year, string pumpID, double RevRate, bool Activate)
{
	double cost=0., P_BEP = 0.;
	int EdgeID = wds[year]->edgeIDtoIndex(pumpID);
	P_BEP = wds[year]->edges.at(EdgeID)->getDoubleProperty("head")*wds[year]->edges.at(EdgeID)->getDoubleProperty("volumeFlowRate")*1000*9.81/1000.;
	cost = 1350+235*P_BEP-1.2*P_BEP*P_BEP;
	if(cost == 0.)
		cout << "ZERO COST DETECTED!!!!!!!!!!!!! ...." << endl;
	if (Activate == true)
	{
		for (int j = year; j < nYear; ++j)
		{
			wds[j]->edges.at(EdgeID)->setDoubleProperty("revolutionNumber",RevRate);
		}
	}
	cout << "cost: " << cost << endl;
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

void BIWS::readCostTable(string fname, string whichCost)
{
 	vector<vector<string>> content;
	vector<string> row;
	string line, word;
 
	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ','))
				row.push_back(word);
			content.push_back(row);
		}
	}
	else
		cout<<"Could not open the file\n";
	if (whichCost == "PipeReplacementCost")
	{
		PipeReplacementCost_Read = true;
		PipeReplacementCost = content;
	}
	if (whichCost == "ValvePlacementCost")
	{
		ValvePlacementCost_Read = true;
		ValvePlacementCost = content;
	}
}