#include <fstream>
#include <sstream>
#include <string>
#include "../../Leakage.h"


using namespace std;
using namespace Eigen;
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//typedef pair<int, int> Edge;

Leakage *wds;


int main(int argc, char *argv[])
{
    string caseFolder = "Networks/";
    ofstream wFile("results.txt");

    double nominalPRVDiameter = 0.2;
    double refA = nominalPRVDiameter*nominalPRVDiameter*M_PI/4.;
    double setting = 10.;
    double minorLoss = 0.02;
    string caseName = argv[1];
    string inFileName = argv[2];
    vector<double> settings, minorLosses;
    cout << "[*]Network: " << caseFolder << caseName << ".inp" << endl;
    wds = new Leakage(caseFolder +caseName + ".inp");
    if(inFileName == "NOE")
    {
    	int counter = 0;
    	for (int i = 0; i < wds->edges.size(); ++i)
    	{
    		if(wds->edges.at(i)->getEdgeStringProperty("type") == "Pipe")
    		{
    			counter += 1;
    		}
    	}
	    wFile << to_string(counter) << '\n';
	    wds->calculateLeakage(); 
	    double summarizedLoss = wds->getSummarizedLeakage();
	    cout << summarizedLoss << '\n';
	}
    else
	{		
		vector<string> fileData = readVectorString(inFileName);
		vector<int> addPRVPipe;
		vector<bool> isStart;
		vector<string> addPRVName;
		vector<double> refCros;
		for(int i=0; i<fileData.size(); i++)
		{
		 stringstream ss(fileData[i]);
		 vector<string> sv;

		 while(ss.good())
		 {
		    string substr;
		    getline(ss,substr,',');
		    sv.push_back(substr);
		 }

		 addPRVPipe.push_back(stoi(sv[0]));
		 isStart.push_back(stoi(sv[1]));
		 addPRVName.push_back("PRV_" + to_string(i));
		 refCros.push_back(refA);
		 settings.push_back(setting);
		 minorLosses.push_back(minorLoss);
		}

		// adding the new PRV valves
		wds->addNewPRVValves(addPRVName, addPRVPipe, isStart, 1000., refCros, settings, minorLosses, 0.);
	    wds->calculateLeakage(); 
	    double summarizedLoss = wds->getSummarizedLeakage();
	    double unservedDemands = wds->calculateUnservedDemands();
	    wFile << summarizedLoss << "," << unservedDemands << '\n';
	}
}
