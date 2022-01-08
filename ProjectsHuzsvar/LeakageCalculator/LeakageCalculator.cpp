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
    string inLossConstantFile;
    string CalculateLoss;
    if (argc > 3 && string(argv[3]) == "Calculate_Loss")
    {
    	CalculateLoss = string(argv[3]);
    }
    else if(argc > 3 && string(argv[3]) == "Modify_C")
    {
    	inLossConstantFile = string(argv[3]);
    }
    cout << " cl: " << CalculateLoss << endl;
    bool duplicate = false;
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
    else if(inFileName == "SEL")
    {
    	double summLength = 0.;
    	for (int i = 0; i < wds->edges.size(); ++i)
    	{
    		cout << i << endl;
    		if(wds->edges.at(i)->getEdgeStringProperty("type") == "Pipe")
    		{
    			summLength += wds->edges.at(i)->getDoubleProperty("length");
    		}
    	}
	    wFile << to_string(summLength) << '\n';
	    cout << "Length calculated!" << '\n';
    }
    else if(CalculateLoss == "Calculate_Loss")
    {
		vector<string> fileData = readVectorString(inFileName);
		cout << fileData[0] << endl;
		double C = stod(fileData[0]);
		wds->setLeakageConstant(C);
		bool succesfulConvergence = wds->calculateLeakage();
		double summarizedLoss = wds->getSummarizedLeakage();
		wFile << to_string(summarizedLoss) << '\n';
		cout << summarizedLoss << '\n';
    }
    else
	{	
		else if(argc > 3 && string(argv[3]) == "Modify_C")
	    {
	    	vector<string> fileData = readVectorString(inLossConstantFile);
			cout << fileData[0] << endl;
			double C = stod(fileData[0]);
			wds->setLeakageConstant(C);
	    }
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
			duplicate = false;
			for (int j = 0; j < i; ++j)
			{
				if (addPRVPipe[i] == stoi(sv[0]))
				{
					duplicate = true;
				}
			}
			if (duplicate == true)
			{
				cout << "duplicate cutted" << endl;
			}
			else
			{
				addPRVPipe.push_back(stoi(sv[0]));
				isStart.push_back(stoi(sv[1]));
				addPRVName.push_back("PRV_" + to_string(i));
				refCros.push_back(refA);
				settings.push_back(setting);
				minorLosses.push_back(minorLoss);
			}
		}
		wds->addNewPRVValves(addPRVName, addPRVPipe, isStart, 1000., refCros, settings, minorLosses, 0.);
		bool succesfulConvergence = wds->calculateLeakage();
		cout << "succesfulConvergence: " << succesfulConvergence << endl;
	    if(succesfulConvergence == true)
	    {
		    double summarizedLoss = wds->getSummarizedLeakage();
		    double unservedDemands = wds->calculateUnservedDemands();
		    wFile << summarizedLoss << "," << unservedDemands << '\n';	    	
	    }
	    else
	    {
	    	double summarizedLoss = 99999.;
		    double unservedDemands = 1.;
		    wFile << summarizedLoss << "," << unservedDemands << '\n';
		    cout << "CONVERGENCE FAILURE!!!" << endl;
	    }
	}
}
