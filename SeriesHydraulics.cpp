#include "SeriesHydraulics.h"

//-------------------------------------------------------------------
SeriesHydraulics::SeriesHydraulics(string fileName) : Sensitivity(fileName)
{
	cout << "1" << endl;
  loadTimeSettings();
  cout << "2" << endl;
  patternIDtoIndex(); // FILLING THE PATTERNS TO THE NODES
  cout << "3" << endl;
}

//-------------------------------------------------------------------
SeriesHydraulics::~SeriesHydraulics(){}

//-------------------------------------------------------------------
void SeriesHydraulics::seriesSolve(ssc seriesSensitivityControl)
{	
	// clearing vectors
	seriesInitialization(); // setting back tanks and clearing outputs

	bool convOk = true;
	time = startTime;
	clockTime = startTime;

	while(time<=endTime && convOk)
	{
		// printing basic info to consol
		if(printLevel > 1)
			seriesInfo();

		updateDemand();
		updatePressurePointPattern();
		updateControl();
		updateRule();

		if(isSensitivityCalculation(seriesSensitivityControl))
		{			
			convOk = calculateSensitivity(seriesSensitivityControl.parameter);
			//saveSensitivityMatrix(seriesSensitivityControl);
			if(printLevel > 1)
				cout << endl << " Sensitivity calculation was done. " << endl;
		}
		else
		{
			convOk = solveSystem();
		}

		// saving head, consumptions for nodes and flows for edges
		saveOutput();

		hydraulicTimeStep = newHydraulicTimeStep();
		time += hydraulicTimeStep;
		updatePool();

		// updating the settings, pools etc.
		//if(time<=endTime)
		//{
		//	updatePool();
		//}

		clockTime = time;
		while(clockTime>=86400.)
			clockTime -= 86400.;
	}

	if(!convOk)
		cout << endl << "!ERROR! Convergence was not achieved in solveSystem!" << endl;
}

//-------------------------------------------------------------------
void SeriesHydraulics::seriesSolve()
{
	ssc seriesSensitivityControl;
	seriesSolve(seriesSensitivityControl);	
}

//-------------------------------------------------------------------
bool SeriesHydraulics::isSensitivityCalculation(ssc seriesSensitivityControl)
{
	bool out = false;
	for(int i=0; i<seriesSensitivityControl.time.size(); i++)
	{
		if(abs(seriesSensitivityControl.time[i] - time)<1e-2)
			out = true;
	}

	return out;
}

//-------------------------------------------------------------------
/*void SeriesHydraulics::saveSensitivityMatrix(ssc seriesSensitivityControl)
{	
	MatrixXd SM;
	int rowSize = seriesSensitivityControl.rowIndex.size();
	int colSize = seriesSensitivityControl.colIndex.size();
	if(seriesSensitivityControl.isPressure)
	{
		if(rowSize == 0 && colSize == 0)
		{
			SM = pressureSensitivity(all,all);
		}
		else if(rowSize != 0 && colSize == 0)
		{
			SM = pressureSensitivity(seriesSensitivityControl.rowIndex,all);
		}
		else if(rowSize == 0 && colSize != 0)
		{
			SM = pressureSensitivity(all,seriesSensitivityControl.colIndex);
		}
		else
		{
			SM = pressureSensitivity(seriesSensitivityControl.rowIndex,seriesSensitivityControl.colIndex);
		}
	}
	else
	{
		if(rowSize == 0 && colSize == 0)
		{
			SM = massFlowRateSensitivity(all,all);
		}
		else if(rowSize != 0 && colSize == 0)
		{
			SM = massFlowRateSensitivity(seriesSensitivityControl.rowIndex,all);
		}
		else if(rowSize == 0 && colSize != 0)
		{
			SM = massFlowRateSensitivity(all,seriesSensitivityControl.colIndex);
		}
		else
		{
			SM = massFlowRateSensitivity(seriesSensitivityControl.rowIndex,seriesSensitivityControl.colIndex);
		}
	}
	seriesSensitivity.push_back(SM);
}
*/
//-------------------------------------------------------------------
void SeriesHydraulics::saveOutput()
{
	vectorTime.push_back(time);
	for(int i=0; i<numberNodes; i++)
	{
		nodes[i]->vectorHead.push_back(nodes[i]->head);
		nodes[i]->vectorConsumption.push_back(nodes[i]->getProperty("consumption"));
		nodes[i]->vectorDemand.push_back(nodes[i]->demand);
		nodes[i]->vectorStatus.push_back(nodes[i]->status);
		nodes[i]->vectorLeakage.push_back(nodes[i]->getProperty("leakage"));
	}
	for(int i=0; i<numberEdges; i++)
	{
		edges[i]->vectorVolumeFlowRate.push_back(edges[i]->volumeFlowRate);
		edges[i]->vectorStatus.push_back(edges[i]->status);
	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::clearOutput()
{
	vectorTime.clear();
	for(int i=0; i<numberNodes; i++)
	{
		nodes[i]->vectorHead.clear();
		nodes[i]->vectorConsumption.clear();
		nodes[i]->vectorDemand.clear();
		nodes[i]->vectorStatus.clear();
		nodes[i]->vectorLeakage.clear();
	}
	for(int i=0; i<numberEdges; i++)
	{
		edges[i]->vectorVolumeFlowRate.clear();
		edges[i]->vectorStatus.clear();
	}
	seriesSensitivity.clear();
}

//-------------------------------------------------------------------
void SeriesHydraulics::seriesInitialization()
{
	clearOutput();
	for(int i=0; i<poolIndex.size(); i++)
	{
		int idx = poolIndex[i];
		double initLevel = edges[idx]->getDoubleProperty("initLevel");
		edges[idx]->setDoubleProperty("waterLevel", initLevel);
	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::saveTimeResult(vector<string> edgeID, vector<string> nodeID, string qUnit, string hUnit)
{
	// saving time
  mkdir(caseName.c_str(),0777);
	FILE* wFile;
	string fn = caseName + "/time.txt";
	wFile = fopen(fn.c_str(),"w");
	fprintf(wFile, " time [s]\n");
	for(int i=0; i<vectorTime.size(); i++)
	{
		fprintf(wFile,"%8.3f\n",vectorTime[i]);
	}
	fclose(wFile);

	// converting units
	double hConvert = 1.;
	if(hUnit == "psi")
		hConvert = 1/0.7032;
	else
		hUnit = "mtr";

	double qConvert = 1000.;
	if(qUnit == "gpm")
		qConvert = 1/0.06309*1000.;
	else
		qUnit = "lps";

	// ID to indecies
	vector<int> idxNode, idxEdge;
	for(int i=0; i<nodeID.size(); i++)
	{	
		int k = nodeIDtoIndex(nodeID[i]);
		if(k!=-1)
	  	idxNode.push_back(k);
	}

	for(int i=0; i<edgeID.size(); i++)
	{	
		int k = edgeIDtoIndex(edgeID[i]);
		if(k!=-1)
	  	idxEdge.push_back(k);
	}

	// saving stuff
	for(int i=0; i<idxNode.size(); i++)
	{
		string fn = caseName + "/nodes";
  	mkdir(fn.c_str(),0777);
		nodes[idxNode[i]]->saveTimeResult(fn, hUnit, qUnit);
	}
	for(int i=0; i<idxEdge.size(); i++)
	{
		string fn = caseName + "/edges";
  	mkdir(fn.c_str(),0777);
		edges[idxEdge[i]]->saveTimeResult(fn, qUnit);
	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::updateDemand()
{	
	for(int i=0; i<numberNodes; i++)
	{
		double demand=0.;
		for(int j=0; j<nodes[i]->vDemand.size(); j++)
		{
			if(nodes[i]->vectorPatternIndex[j] == -1)
			{
				demand += nodes[i]->vDemand[j];
			}
			else
			{
				int timeIndex = time/patternTimeStep;
				int ptIndex = nodes[i]->vectorPatternIndex[j];
				// making patternvalues periodic
				while(timeIndex>=patternValue[ptIndex].size())
					timeIndex -= patternValue[ptIndex].size();
				demand += nodes[i]->vDemand[j]*patternValue[ptIndex][timeIndex];
			}
		}
		nodes[i]->demand = demand;
	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::updatePressurePointPattern()
{
	for(int i=0; i<presIndex.size(); i++)
	{
		int idx = presIndex[i];
		if(edges[idx]->getIntProperty("patternIndex") != -1)
		{
			double headNominal = edges[idx]->getDoubleProperty("headNominal");
			int timeIndex = time/patternTimeStep;
			int ptIndex = edges[idx]->getIntProperty("patternIndex");
			// making patternvalues periodic
			while(timeIndex>=patternValue[ptIndex].size())
				timeIndex -= patternValue[ptIndex].size();

			double head = headNominal * patternValue[ptIndex][timeIndex] - edges[idx]->startHeight;
			edges[idx]->setDoubleProperty("head", head);
		}
	}
}

//-------------------------------------------------------------------
double SeriesHydraulics::newHydraulicTimeStep()
{	
	double dt = hydraulicTimeStepOriginal, dt2 = hydraulicTimeStepOriginal;

	// Checking the [CONTROLS]
	for(int i=0; i<controlEdgeID.size(); i++)
	{
		dt2 = hydraulicTimeStepOriginal;
		double wl = edges[controlNodeIndex[i]]->getDoubleProperty("waterLevel");
		double aref = edges[controlNodeIndex[i]]->referenceCrossSection;
		double set = controlValue[i];
		double vf = edges[controlNodeIndex[i]]->volumeFlowRate;
		if(controlType[i] == "Pool")
		{
			if(vf > 0. && controlAbove[i] && wl < set) // FILLING THE TANK
				dt2 = (set - wl) * aref / vf;
			else if(vf < 0. && !controlAbove[i] && wl > set)
				dt2 = (set - wl) * aref / vf;

			if(dt2<hydraulicTimeStepOriginal/100.)
			{
				if(printLevel > 1)
					cout << endl << "!WARNING! dt is too small at Control " << controlEdgeID[i] << endl;
			}
			if(dt2<dt)
				dt = dt2*1.0;
		}
		else
			cout << endl << "!WARNING! [CONTROLS] control type unkown: " << controlType[i] << endl;
	}

	// Checking the tank fillups
	for(int i=0; i<poolIndex.size(); i++)
	{	
		dt2 = hydraulicTimeStepOriginal;
		int idx = poolIndex[i];
		if(edges[i]->status == 1)
		{
			double vf = edges[idx]->volumeFlowRate;
			double wl = edges[idx]->getDoubleProperty("waterLevel");
			double max_wl = edges[idx]->getDoubleProperty("maxLevel");
			double min_wl = edges[idx]->getDoubleProperty("minLevel");
			double aref = edges[idx]->referenceCrossSection;

			if(vf > 0. && wl < max_wl) // VF > 0. -> FILLING UP
				dt2 = (max_wl - wl) * aref / vf;
			else if(vf < 0. && wl > min_wl) // VF < 0. -> EMPTYING
				dt2 = (min_wl - wl) * aref / vf;

			if(dt2<hydraulicTimeStepOriginal/100.)
			{
				if(printLevel > 1)
					cout << endl << "!WARNING! dt is too small at Tank Fillup: " << i << " th edge, name: " << edges[idx]->name << endl;
			}
			if(dt2<dt)
				dt = dt2*1.0;
		}
	}

	// Checking the [RULES]
	for(int i=0; i<rules.size(); i++)
	{
		for(int j=0; j<rules[i]->conditionID.size(); j++)
		{	
			dt2 = hydraulicTimeStepOriginal;
			if(rules[i]->conditionID[j] == "CLOCKTIME" || rules[i]->conditionID[j] == "TIME")
			{
				double eps = 0;
				string relation = rules[i]->conditionRelation[j];
				if(relation == "<=" || relation == ">")
					eps = hydraulicTimeStepOriginal/10.;

				if(rules[i]->conditionID[j] == "CLOCKTIME")
					dt2 = rules[i]->conditionValue[j] - clockTime;
				else
					dt2 = rules[i]->conditionValue[j] - time;

				bool cond = dt2 < dt && ((dt2 >= 0. && eps > 0.) || (dt2 > 0. && eps == 0.));
				if(cond)
					dt = dt2 + eps;

				if(dt2<0. && cond)
				{
					cout << endl << "!WARNING! Something is wrong, dt is negative at Rule: " << i << ", name: " << rules[i]->ID << endl;
					cin.get();
				}
			}

			if(rules[i]->conditionType[j] == "NODE")
			{	
				int idx = rules[i]->conditionIndex[j];
				double set = rules[i]->conditionValue[j];
				double vf = edges[idx]->volumeFlowRate;
				double wl = edges[idx]->getDoubleProperty("waterLevel");
				double aref = edges[idx]->referenceCrossSection;

				bool above;
				string rel = rules[i]->conditionRelation[j];
				if(rel == "ABOVE" || rel == ">" || rel == ">=")
					above = true;
				else
					above = false;

				if(vf > 0. && above && wl < set) // VF > 0. -> FILLING UP
					dt2 = (set-wl) * aref / vf;
				else if(vf < 0. && !above && wl > set) // VF < 0. -> EMPTYING
					dt2 = (set-wl) * aref / vf;

				bool cond = checkCondition(i, dt2*1.001, dt2*vf/aref);

				if(dt2<0. && cond)
				{
					cout << endl << "!WARNING! Something is wrong, dt is negative at Rule: " << i << ", name: " << rules[i]->ID << endl;
					cin.get();
				}
				else if(dt2<dt && cond)
					dt = dt2;
			}
		}
	}

	// setting back timestep to hydraulicTimeStep*i where i is integer
	dt2 = 0.;
	while(dt2 < time)
		dt2 += hydraulicTimeStepOriginal;
	double res = dt2 - time;
	if(res > 1e-1 && dt == hydraulicTimeStepOriginal) // there was a change in timeStep before
		dt = res;

	// avoiding too small time steps
	//if(dt < hydraulicTimeStepOriginal/100.)
	//{
	//	dt = hydraulicTimeStepOriginal/100.;
	//	cout << endl << "!WARNING! Time step is too small, dt: " << dt << "\n System clock: " << secondsToTime(time) << endl;
	//}

	return dt;
}

//-------------------------------------------------------------------
void SeriesHydraulics::updatePool()
{
	for(int i=0; i<poolIndex.size(); i++)
	{
		int ei = poolIndex[i]; // index of pool in edges list
		int ni = edges[ei]->startNodeIndex; // node index of pool
		double wl = edges[ei]->getDoubleProperty("waterLevel");
		double vf = edges[ei]->volumeFlowRate;
		double aref = edges[ei]->referenceCrossSection;
		double min_wl = edges[ei]->getDoubleProperty("minLevel");
		double max_wl = edges[ei]->getDoubleProperty("maxLevel");

		// new water level
		wl += vf * hydraulicTimeStep / aref;

		// controlling the in and outgoing edges
		if(wl <= min_wl + headTolerance) // gets EMPTY
		{
			edges[ei]->setDoubleProperty("waterLevel",min_wl);
			for(int j=1; j<nodes[ni]->edgeOut.size(); j++)
			{
				int idx = nodes[ni]->edgeOut[j];
				if(edges[idx]->typeCode == 2) // pump
				{
					if(!edges[idx]->status_fix)
						edges[idx]->status = 0; // closing outpointing pump
				}
				else // not pump
				{
					if(edges[idx]->volumeFlowRate > 0.)
					{
						if(!edges[idx]->status_fix)
							deleteFromPipeCVIndex(idx); // opening inflowing pipes/valves
					}
					else
					{
						if(!edges[idx]->status_fix)
							addToPipeCVIndex(-idx); // closing outflowing pipes/valves
					}
				}
			}
			for(int j=0; j<nodes[ni]->edgeIn.size(); j++)
			{
				int idx = nodes[ni]->edgeIn[j];
				if(edges[idx]->typeCode == 2) // pump
				{
					if(!edges[idx]->status_fix)
						edges[idx]->status = 1; // just making sure it is open
				}
				else // not pump
				{
					if(edges[idx]->volumeFlowRate > 0.)
					{
						if(!edges[idx]->status_fix)
							addToPipeCVIndex(idx); // closing outflowing pipes/valves
					}
					else
					{
						if(!edges[idx]->status_fix)
							deleteFromPipeCVIndex(idx); // opening inflowing pipes/valves
					}
				}
			}
		}
		else if(wl >= max_wl - headTolerance) // gets FULL
		{
			edges[ei]->setDoubleProperty("waterLevel",max_wl);
			bool doOverflow = edges[ei]->getBoolProperty("doOverflow");

			if(!doOverflow) // overflow is not allowed (epanet default)
			{
				for(int j=1; j<nodes[ni]->edgeOut.size(); j++)
				{
					int idx = nodes[ni]->edgeOut[j];
					if(edges[idx]->typeCode == 2) // pump
					{
						if(!edges[idx]->status_fix)
							edges[idx]->status = 1; // making sure the pump is operating
					}
					else // not pump
					{
						if(edges[idx]->volumeFlowRate < 0.)
						{
							if(!edges[idx]->status_fix)
								deleteFromPipeCVIndex(idx); // opening outgoing flow
						}
						else
						{
							if(!edges[idx]->status_fix)
								addToPipeCVIndex(idx); // closing ingoing flow
						}
					}
				}
				for(int j=0; j<nodes[ni]->edgeIn.size(); j++)
				{
					int idx = nodes[ni]->edgeIn[j];
					if(edges[idx]->typeCode == 2) // pump
					{
						if(!edges[idx]->status_fix)
							edges[idx]->status = 0; // closing ingoing pump flow
					}
					else
					{
						if(edges[idx]->volumeFlowRate < 0.)
						{
							if(!edges[idx]->status_fix)
								addToPipeCVIndex(-idx); // closing ingoing flow
						}
						else
						{
							if(!edges[idx]->status_fix)
								deleteFromPipeCVIndex(idx); // opening outgoing flow
						}
					}
				}
			}
		}
		else // normal conditions
		{
			edges[ei]->setDoubleProperty("waterLevel",wl);
			// opening up everyting, active valves??
			for(int j=1; j<nodes[ni]->edgeOut.size(); j++)
			{
				int idx = nodes[ni]->edgeOut[j];
				if(!edges[idx]->status_fix)
				{
					if(edges[idx]->typeCode != 2)
					{
						deleteFromPipeCVIndex(idx);
					}
					else
					{
						edges[idx]->status = 1;
					}
				}
			}
			for(int j=0; j<nodes[ni]->edgeIn.size(); j++)
			{
				int idx = nodes[ni]->edgeIn[j];
				if(!edges[idx]->status_fix)
				{
					if(edges[idx]->typeCode != 2)
					{
						deleteFromPipeCVIndex(idx);
					}
					else
					{
						edges[idx]->status = 1;
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::addToPipeCVIndex(int idx)
{
	bool gotIt=false;
	int i=0;
	while(i<pipeCVIndex.size() && !gotIt)
	{
		if(pipeCVIndex[i] == idx)
			gotIt = true;
		i++;
	}
	if(!gotIt)
		pipeCVIndex.push_back(idx);
}

//-------------------------------------------------------------------
void SeriesHydraulics::deleteFromPipeCVIndex(int idx)
{
	int gotIt=-1;
	int i=0;
	while(i<pipeCVIndex.size() && gotIt==-1)
	{
		if(pipeCVIndex[i] == abs(idx))
			gotIt = i;
		i++;
	}
	if(gotIt != -1)
		pipeCVIndex.erase(pipeCVIndex.begin()+gotIt);
}

//-------------------------------------------------------------------
void SeriesHydraulics::updateControl()
{	
	vector<int> changedIndex, oldStatus;
	int status;
	for(int i=0; i<controlEdgeID.size(); i++)
	{	
		status = edges[controlEdgeIndex[i]]->status;
		if(controlType[i] == "Pool")
		{
			double wl = edges[controlNodeIndex[i]]->getDoubleProperty("waterLevel");
			if(controlAbove[i])
			{
				if(wl >= controlValue[i] - headTolerance)
				{
					edges[controlEdgeIndex[i]]->status = controlStatus[i];
				}
			}
			else
			{
				if(wl <= controlValue[i] + headTolerance)
				{
					edges[controlEdgeIndex[i]]->status = controlStatus[i];
				}
			}
		}
		else
			cout << endl << "!WARNING! [CONTROLS] control type unkown: " << controlType[i] << endl;

		if(edges[controlEdgeIndex[i]]->status != status) // i.e. there is a change
		{
			changedIndex.push_back(i);
			oldStatus.push_back(status);
			if(edges[controlEdgeIndex[i]]->status == 0)
			{
				x(controlEdgeIndex[i]) = 0.0;
			}
		}
	}

	if(printLevel > 1)
	{
		if(changedIndex.size() > 0)
		{
	    printf(" *************************** CHANGE BY CONTROLS ***************************\n");
	    printf(" | Edge ID              | old status -> new status | Node ID              | \n");
		  for(int i=0; i<changedIndex.size(); i++)
		  {
		    int idx = changedIndex[i];
		    printf(" | %-20s | %10i -> %-10i | %20s | \n", controlEdgeID[idx].c_str(), oldStatus[i], edges[controlEdgeIndex[idx]]->status, controlNodeID[idx].c_str());
		  }
		  cout << endl;
		}
	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::updateRule()
{
	vector<int> changedIndex, changedRule;
	vector<double> oldStatting, newStatting;
	for(int i=0; i<rules.size(); i++)
	{	
		// HANDLING CONDITIONS
		// ONLY AND IS CONSIDERED, TODO: OR, see ROSSMANN, 2000
		bool cond = checkCondition(i, 0., 0.);

		// EXECUTING THE ACTIONS
		for(int j=0; j<rules[i]->actionID.size(); j++)
		{
			bool act = cond*!rules[i]->actionCounter[j] + !cond*rules[i]->actionCounter[j]; // THEN + ELSE
			if(act)
			{
				if(rules[i]->actionSetting[j]) // SETTING
				{
					if(edges[rules[i]->actionIndex[j]]->setting != rules[i]->actionValue[j])
					{
						oldStatting.push_back(edges[rules[i]->actionIndex[j]]->setting);

						edges[rules[i]->actionIndex[j]]->setting = rules[i]->actionValue[j];
						edges[rules[i]->actionIndex[j]]->status_fix = false;

						newStatting.push_back(edges[rules[i]->actionIndex[j]]->setting);
						changedIndex.push_back(rules[i]->actionIndex[j]);
						changedRule.push_back(i);
					}
				}
				else // STATUS
				{
					if(edges[rules[i]->actionIndex[j]]->status != rules[i]->actionValue[j])
					{
						oldStatting.push_back(edges[rules[i]->actionIndex[j]]->status);

						edges[rules[i]->actionIndex[j]]->status = rules[i]->actionValue[j];
						edges[rules[i]->actionIndex[j]]->status_fix = false;

						newStatting.push_back(edges[rules[i]->actionIndex[j]]->status);
						changedIndex.push_back(rules[i]->actionIndex[j]);
						changedRule.push_back(i);

						// setting new initial conditions for the solver
						if(edges[rules[i]->actionIndex[j]]->status == 0)
							x(rules[i]->actionIndex[j]) = 0.0;
					}
				}
			}
		}
	}

	if(printLevel > 1)
	{
		if(changedIndex.size() > 0)
		{
	    printf(" ********************************* CHANGE BY RULES **********************************\n");
	    printf(" | Edge ID              | old statti -> new statti | Based on RULE                  | \n");
		  for(int i=0; i<changedIndex.size(); i++)
		  {
		    int idx = changedIndex[i];
		    printf(" | %-20s | %10.1f -> %-10.1f | %30s | \n", edges[idx]->name.c_str(), oldStatting[i], newStatting[i], rules[changedRule[i]]->ID.c_str());
		  }
		  cout << endl;
		}
	}
}

//-------------------------------------------------------------------
bool SeriesHydraulics::checkCondition(int i, double dt, double dwl)
{
	bool cond = true;
	for(int j=0; j<rules[i]->conditionID.size(); j++)
	{
		if(rules[i]->conditionType[j] == "SYSTEM")
		{
			if(rules[i]->conditionID[j] == "CLOCKTIME")
			{
				cond *= booleanWithString(rules[i]->conditionRelation[j], clockTime + dt, rules[i]->conditionValue[j]);
			}
			else if(rules[i]->conditionID[j] == "TIME")
			{
				cond *= booleanWithString(rules[i]->conditionRelation[j], time + dt, rules[i]->conditionValue[j]);
			}
		}
		else if(rules[i]->conditionType[j] == "NODE") // actually waterlevel of a tank (that is an edge in staci3)
		{
			cond *= booleanWithString(rules[i]->conditionRelation[j], edges[rules[i]->conditionIndex[j]]->getDoubleProperty("waterLevel") + dwl, rules[i]->conditionValue[j]);
		}
		else
		{
			cout << endl << "!WARNING! [Rules] Condition type is unkown: " << rules[i]->conditionType[j] << " at Rule " << rules[i]->ID << endl;
		}
	}

	return cond;
}

//-------------------------------------------------------------------
bool SeriesHydraulics::booleanWithString(string op, double left, double right)
{
	bool out = false;
	if((op == "<=") && (left <= right))
		out = true;

	if((op == ">=") && (left >= right))
		out = true;

	if((op == "<") && (left < right))
		out = true;

	if((op == ">") && (left > right))
		out = true;

	if((op == "=") && (left == right))
		out = true;

	if((op == "ABOVE") && (left > right))
		out = true;

	if((op == "BELOW") && (left < right))
		out = true;
	return out;
}

//-------------------------------------------------------------------
void SeriesHydraulics::seriesInfo()
{
	cout << endl << "[*] time: " << secondsToTime(time) << endl;
}

//-------------------------------------------------------------------
void SeriesHydraulics::timeTableNode(vector<string> ID, string unit)
{
	vector<int> idx;
	for(int i=0; i<ID.size(); i++)
	{	
		int k = nodeIDtoIndex(ID[i]);
		if(k!=-1)
	  	idx.push_back(k);
	}
	timeTableNode(idx,unit);
}
//-------------------------------------------------------------------
void SeriesHydraulics::timeTableNode(vector<int> idx, string unit)
{
	double convertUnit = 1.;
	if(unit == "psi")
		convertUnit = 1/0.7032;
	else
		unit = " m ";

	cout << endl;
	printf(" |   ID   |");
	for(int i=0; i<idx.size(); i++)
		printf(" %-19s |",nodes[idx[i]]->name.c_str());

	printf("\n |  time  |");
	for(int i=0; i<idx.size(); i++)
		printf(" head [%3s] | status |",unit.c_str());// | status |");

	printf("\n |--------");
	for(int i=0; i<idx.size(); i++)
		printf("+---------------------");
	printf("|");

	for(int i=0; i<vectorTime.size(); i++)
	{
		printf("\n");
		int hour = vectorTime[i]/3600.;
		int minute = (vectorTime[i] - hour*3600.)/60.;
		if(minute < 10)
			printf(" | %3i:0%1i |", hour, minute);
		else
			printf(" | %3i:%2i |", hour, minute);

		for(int j=0; j<idx.size(); j++)
		{
			printf(" %10.3f |", nodes[idx[j]]->vectorHead[i]*convertUnit);
			if(nodes[idx[j]]->vectorStatus[i] == 1)
				printf("  open  |");
			else if(nodes[idx[j]]->vectorStatus[i] == 0)
				printf(" closed |");
		}
	}
	cout << endl;
}
//-------------------------------------------------------------------
void SeriesHydraulics::timeTableEdge(vector<string> ID, string unit)
{
	vector<int> idx;
	for(int i=0; i<ID.size(); i++)
	{	
		int k = edgeIDtoIndex(ID[i]);
		if(k!=-1)
	  	idx.push_back(k);
	}
	timeTableEdge(idx, unit);
}
//-------------------------------------------------------------------
void SeriesHydraulics::timeTableEdge(vector<int> idx, string unit)
{
	double convertUnit = 1000.;
	if(unit == "gpm")
		convertUnit = 1/0.06309*1000.;
	else
		unit = "lps";

	cout << endl;
	printf(" |   ID   |");
	for(int i=0; i<idx.size(); i++)
		printf(" %-19s |",edges[idx[i]]->name.c_str());

	printf("\n |  time  |");
	for(int i=0; i<idx.size(); i++)
		printf(" flow [%3s] | status |",unit.c_str());

	printf("\n |--------");
	for(int i=0; i<idx.size(); i++)
		printf("+---------------------");
	printf("|");

	for(int i=0; i<vectorTime.size(); i++)
	{
		printf("\n");
		int hour = vectorTime[i]/3600.;
		int minute = (vectorTime[i] - hour*3600.)/60.;
		if(minute < 10)
			printf(" | %3i:0%1i |", hour, minute);
		else
			printf(" | %3i:%2i |", hour, minute);

		for(int j=0; j<idx.size(); j++)
		{
			printf(" %10.3f |", edges[idx[j]]->vectorVolumeFlowRate[i]*convertUnit);
			if(edges[idx[j]]->vectorStatus[i] == 2)
				printf(" active |");
			else if(edges[idx[j]]->vectorStatus[i] == 1)
				printf("  open  |");
			else if(edges[idx[j]]->vectorStatus[i] == 0)
				printf(" closed |");
			else if(edges[idx[j]]->vectorStatus[i] == -1)
				printf(" CLOSED |");
		}
	}
	cout << endl;
}

//-------------------------------------------------------------------
string SeriesHydraulics::secondsToTime(double seconds)
{
	string out;
	int hour = seconds/3600.;
	int minute = (seconds - hour*3600.)/60.;
	if(minute < 10)
		out = to_string(hour) + ":0" + to_string(minute);
	else
		out = to_string(hour) + ":" + to_string(minute);
	return out;
}

//-------------------------------------------------------------------
double SeriesHydraulics::timeToSeconds(string s1, string s2)
{
	double time;
	if(s2 == "") // HANDLING 12:00, 14:30, 09:45 etc.
	{
		time = stod(s1.substr(0,s1.find(":")),0) * 3600. + stod(s1.substr(s1.find(":")+1,s1.length()),0) * 60.;
	}
	else // HANDLING 9:45 AM, 2 Days, 4 Hours etc.
	{
		if(s2 == "DAYS" || s2 == "Days" || s2 == "days")
			time = stod(s1,0) * 86400.;
		else if(s2 == "HOURS" || s2 == "Hours" || s2 == "hours")
			time = stod(s1,0) * 3600.;
		else if(s2 == "MIN" || s2 == "Min" || s2 == "min")
			time = stod(s1,0) * 60.;
		else if(s2 == "AM" || s2 == "am")
		{
			double t1 = stod(s1.substr(0,s1.find(":")),0);
			if(t1 >= 12.)
				t1 -= 12.;
			if(s1.find(":") <= s1.length())
				time =  t1* 3600. + stod(s1.substr(s1.find(":")+1,s1.length()),0) * 60.;
			else
				time =  t1* 3600.;
		}
		else if(s2 == "PM" || s2 == "pm")
		{
			double t1 = stod(s1.substr(0,s1.find(":")),0);
			if(t1 != 12)
				t1 += 12.;
			
			if(s1.find(":") <= s1.length())
				time =  t1* 3600. + stod(s1.substr(s1.find(":")+1,s1.length()),0) * 60.;
			else
				time =  t1* 3600.;
		}
		else
		{
			cout << endl << "!WARNING! Unkown input: " << s2 << "\nAvailable inputs: DAYS | HOURS | MIN | AM | PM" << endl;
			time = 0.;
		}
	}
	return time;
}

//-------------------------------------------------------------------
void SeriesHydraulics::patternIDtoIndex()
{
 	for(int i=0; i<numberNodes; i++)
 	{
 		if(nodes[i]->vectorPatternID.size() == 1 && nodes[i]->vectorPatternID[0] == "")
 		{
 			nodes[i]->vectorPatternID[0] = basePattern;
 		}
 		for(int j=0; j<nodes[i]->vectorPatternID.size(); j++)
 		{	
 			nodes[i]->vectorPatternIndex.resize(nodes[i]->vectorPatternID.size(),-1);
 			for(int k=0; k<patternID.size(); k++)
 			{
 				if(nodes[i]->vectorPatternID[j] == patternID[k])
 				{
 					nodes[i]->vectorPatternIndex[j] = k;
 					break;
 				}
 			}
 		}
 	}

 	for(int i=0; i<presIndex.size(); i++)
 	{
 		int idx = presIndex[i];
		for(int k=0; k<patternID.size(); k++)
		{
			if(edges[idx]->getStringProperty("patternID") == patternID[k])
			{
				edges[idx]->setIntProperty("patternIndex",k);
				break;
			}
 		}
 	}
}

//-------------------------------------------------------------------
void SeriesHydraulics::loadTimeSettings()
{
	ifstream fileIn;
  fileIn.open(definitionFile);
  string line;
  if(fileIn.is_open())
  {
	  getline(fileIn,line);
	  while(true)
	  {
	  	while(line[0] != '[')
	  		getline(fileIn,line);

	  	string line_old = line;

	  	if(line.length() >= 6 && line.substr(0,5) == "[END]") // reaching the end of the file
	  		break;

	  	// LOADING THE DEMANDS AND PATTERNS FOR SERIES CALC
	  	if(line.substr(0,11) == "[JUNCTIONS]")
	  	{
	  		int i=0;
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						if(sv.size()==3)
						{
							nodes[i]->vDemand.push_back(demandUnit*stod(sv[2],0));
							nodes[i]->vectorPatternID.push_back("");
						}
						if(sv.size()>3)
						{
							nodes[i]->vDemand.push_back(demandUnit*stod(sv[2],0));
							nodes[i]->vectorPatternID.push_back(sv[3]);
						}
						i++;
					}
				}
			}

			// LOADING THE PATTERNS OF RESERVOIRS FOR SERIES CALC
	  	if(line.substr(0,12) == "[RESERVOIRS]")
	  	{
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						string patternID;
						if(sv.size()>2)
							patternID = sv[2];

						for(int i=0; i<presIndex.size(); i++)
						{
							int idx = presIndex[i];
							if(edges[idx]->name == sv[0])
								edges[idx]->setStringProperty("patternID",patternID);
						}
					}
				}
			}

			// LOADING THE DEMANDS FOR THE NODES
			// Here we load every demand with pattern data
	  	if(line.substr(0,9) == "[DEMANDS]")
	  	{
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						if(sv.size()>=2)
						{
							string id = sv[0];
							for(int i=0; i<numberNodes; i++)
							{
								if(id == nodes[i]->name)
								{
									nodes[i]->vDemand.push_back(demandUnit*stod(sv[1],0));
									if(sv.size()>=3)
										nodes[i]->vectorPatternID.push_back(sv[2]);
									else
										nodes[i]->vectorPatternID.push_back("");
									break;
								}
								if(i==numberNodes-1)
									cout << "\n!WARNING! demand ID " << id << " was not found" << endl;
							}
						}
					}
				}
			}

	  	// LOADING THE DEMANDS AND PATTERNS FOR SERIES CALC
			if(line.substr(0,10) == "[PATTERNS]")
			{
				vector<double> v;
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
					{
						patternValue.push_back(v);
						break;
					}

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						if(patternID.size()==0)
						{
							patternID.push_back(sv[0]);
							for(vector<string>::iterator i=sv.begin()+1; i!=sv.end(); i++)
								v.push_back(stod(*i,0));
						}
						else
						{
							if(sv[0] == patternID.back())
							{
								for(vector<string>::iterator i=sv.begin()+1; i!=sv.end(); i++)
									v.push_back(stod(*i,0));
							}
							else
							{
								patternValue.push_back(v);
								v.clear();
								patternID.push_back(sv[0]);
								for(vector<string>::iterator i=sv.begin()+1; i!=sv.end(); i++)
									v.push_back(stod(*i,0));
							}
						}
					}
				}
			}
			// TIME SETTINGS OF SERIES CALCULATIONS
			if(line.substr(0,7) == "[TIMES]")
			{
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						if(sv.size()>=2)
						{
							if(sv[0] == "Duration" || sv[0] == "DURATION")
							{
								if(sv.size() == 2)
									endTime = timeToSeconds(sv[1],"");
								else
									endTime = timeToSeconds(sv[1],sv[2]);
							}
							if((sv[0] == "HYDRAULIC" && sv[1] == "TIMESTEP") || (sv[0] == "Hydraulic" && sv[1] == "Timestep"))
							{
								if(sv.size() == 3)
									hydraulicTimeStepOriginal = timeToSeconds(sv[2],"");
								else
									hydraulicTimeStepOriginal = timeToSeconds(sv[2],sv[3]);
							}
							if((sv[0] == "PATTERN" && sv[1] == "TIMESTEP") || (sv[0] == "Pattern" && sv[1] == "Timestep"))
							{
								if(sv.size() == 3)
									patternTimeStep = timeToSeconds(sv[2],"");
								else
									patternTimeStep = timeToSeconds(sv[2],sv[3]);
							}
							if((sv[0] == "START" && sv[1] == "CLOCKTIME") || (sv[0] == "Start" && sv[1] == "ClockTime"))
							{
								if(sv.size() == 3)
									startTime = timeToSeconds(sv[2],"");
								else
									startTime = timeToSeconds(sv[2],sv[3]);
							}
						}
					}
				}
			}

			// LOADING THE DEMANDS AND PATTERNS FOR SERIES CALC
	  	if(line.substr(0,9) == "[OPTIONS]")
	  	{
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						if(sv.size()>=2)
						{
							if(sv[0] == "Pattern" || sv[0] == "PATTERN")
							{
								if(sv.size()>1)
									basePattern = sv[1];
							}
						}
					}
				}
			}

	  	// LOADING CONTROLS
  		// TODO: do this in more general way, handling time contorls as well (maybe struct or class)
	  	if(line.substr(0,10) == "[CONTROLS]")
	  	{
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						controlEdgeID.push_back(sv[1]);
						for(int i=0; i<numberEdges; i++)
						{
							if(controlEdgeID.back() == edges[i]->getEdgeStringProperty("name"))
							{
								controlEdgeIndex.push_back(i);
								break;
							}
							if(i == numberEdges-1)
								cout << endl << "!WARNING! Edge in [CONTROLS] " << controlEdgeID.back() << " was not found in edge list." << endl;
						}

						if(sv[2] == "OPEN" || sv[2] == "Open" || sv[2] == "open")
							controlStatus.push_back(true);
						else if(sv[2] == "CLOSED" || sv[2] == "Close" || sv[2] == "close")
							controlStatus.push_back(false);
						else
						{
							controlStatus.push_back(true);
							cout << endl << "!WARNING! [CONTROLS] Unkown input: " << sv[2] << " at " << sv[1] << ". Changed to OPEN." << endl;
						}

						controlNodeID.push_back(sv[5]);

						for(int i=0; i<numberEdges; i++)
						{
							if(edges[i]->getEdgeStringProperty("type") == "Pool")
							{
								if(controlNodeID.back() == edges[i]->getEdgeStringProperty("name"))
								{
									controlNodeIndex.push_back(i);
									controlType.push_back("Pool");
									break;
								}
							}
						}
						if(controlNodeIndex.size() != controlNodeID.size())
							cout << endl << "!WARNING! Node in [CONTROLS] " << controlNodeID.back() << " was not found in node list." << endl;

						if(sv[6] == "ABOVE" || sv[6] == "Above" || sv[6] == "above")
							controlAbove.push_back(true);
						else if(sv[6] == "BELOW" || sv[6] == "Below" || sv[6] == "below")
							controlAbove.push_back(false);
						else
						{
							controlAbove.push_back(true);
							cout << endl << "!WARNING! [CONTROLS] Unkown input: " << sv[6] << " at " << sv[0] << "  " << sv[1] << ". Changed to ABOVE." << endl;
						}

						if(unit == "SI")
							controlValue.push_back(stod(sv[7],0));
						else // US
							controlValue.push_back(0.3048*stod(sv[7],0));
					}
				}
			}

			// RULES
			// TODO: lot of other setting/cases are not considered yet
			if(line.substr(0,7) == "[RULES]"){
				int numberRule=-1;
				while(getline(fileIn,line))
				{
					while(line.length() <= 1) // skipping the empty lines
		  			getline(fileIn,line);

					if(line[0] == '[') // reaching the end of this section
						break;

					if(line[0] != ';')
					{
						vector<string> sv=line2sv(line);
						if(sv[0] == "RULE")
						{
							rules.push_back(new Rule(sv[1]));
							numberRule++;
						}
						if(sv[0] == "IF" || sv[0] == "AND") // TODO: AND / OR handling
						{
							if(sv[1] == "SYSTEM")
							{
								rules[numberRule]->conditionType.push_back(sv[1]); // SYSTEM
								rules[numberRule]->conditionID.push_back(sv[2]); // CLOCKTIME or TIME
								rules[numberRule]->conditionIndex.push_back(-1); // We do not need this in this case
								rules[numberRule]->conditionRelation.push_back(sv[3]);
								rules[numberRule]->conditionValue.push_back(timeToSeconds(sv[4],sv[5]));
							}
							if(sv[1] == "NODE")
							{
								rules[numberRule]->conditionType.push_back(sv[1]);
								rules[numberRule]->conditionID.push_back(sv[2]);
								rules[numberRule]->conditionIndex.push_back(edgeIDtoIndex(sv[2]));
								rules[numberRule]->conditionRelation.push_back(sv[4]);
								double value = stod(sv[5],0);
								if(unit == "US") // converting from feet to meter
									value = 0.3048*value;
								rules[numberRule]->conditionValue.push_back(value);
							}
						}
						if(sv[0] == "THEN" || sv[0] == "ELSE") // only one THEN action and one ELSE action is considered
						{
							rules[numberRule]->actionType.push_back(sv[1]); // LINK, NODE, PUMP
							rules[numberRule]->actionID.push_back(sv[2]); // ID of LINK, NODE, PUMP

							int idx = edgeIDtoIndex(sv[2]);
							if(sv[1] == "LINK" || sv[1] == "PUMP")
								rules[numberRule]->actionIndex.push_back(idx); // Index of LINK, NODE, PUMP
							else
								cout << endl << "!WARNING! Unkown attribute: " << sv[2] << "  at Rule: " << rules[numberRule]->ID << endl;

							if(sv[3] == "SETTING")
							{
								rules[numberRule]->actionSetting.push_back(true);

								if(edges[idx]->typeCode == 6)
									rules[numberRule]->actionValue.push_back(demandUnit*stod(sv[5],0));
								else if(edges[idx]->typeCode == 3)
									rules[numberRule]->actionValue.push_back(headUnit*stod(sv[5],0));
								else
									cout << endl << "!WARNING! Not handeld input in Rules Setting, edge type: " << edges[idx]->type << endl;
							}
							else if(sv[3] == "STATUS")
							{
								rules[numberRule]->actionSetting.push_back(false);
								if(sv[5] == "OPEN")
									rules[numberRule]->actionValue.push_back(1.);
								else
									rules[numberRule]->actionValue.push_back(0.);
							}
							else
								cout << endl << "!WARNING! Unkown attribute: " << sv[3] << "  at Rule: " << rules[numberRule]->ID << endl;

							if(sv[0] == "THEN")
								rules[numberRule]->actionCounter.push_back(false);
							else if(sv[0] == "ELSE")
								rules[numberRule]->actionCounter.push_back(true);
						}
					}
				}
			}
			if(line == line_old) // [XYZ] is not intresting
				if(!getline(fileIn,line)) // if we reached the end of file break
					break;
		}//end of while(getline)
	}//end of if(fileIn.is_open())
} 


//-------------------------------------------------------------------
double SeriesHydraulics::getHydraulicTimeStep()
{
	return hydraulicTimeStep;
}
