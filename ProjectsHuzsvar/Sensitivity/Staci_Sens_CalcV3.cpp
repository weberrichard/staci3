#include "iostream"
#include <fstream>
#include <sstream>
#include "../../Sensitivity.h"
#include <cmath>
#include <algorithm>
#include <complex> 
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <utility>
#include </usr/include/eigen3/Eigen/Eigen>
#include </usr/include/eigen3/Eigen/Dense>
#include <stdio.h>
#include <chrono> 
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>
//-------------------------------------------------------------------------------------//
// for testing on real: /home/namrak/HTamas/halozatok/Sopron/
// VIZ-SOPTVR-J-55-input_mod
//-------------------------------------------------------------------------------------//

using namespace std;
using namespace Eigen;
using namespace chrono;
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//typedef pair<int, int> Edge;

int main(int argc, char *argv[])
{
    cout << "[*]Hydraulic solver started...." << endl;
    //------------------------------------------------------------------------Staci init----------------------------------------------------------------//
    string caseFolder = "../../Networks/Sopron/";
    string caseName = "acsad";
    Sensitivity *wds = new Sensitivity(caseFolder + caseName + ".inp");
    //-----------------------------------------------------------------------Staci init End-------------------------------------------------------------//
    wds->calculateSensitivity("demand");
    MatrixXd a = wds->getPSensMatrix();
    cout << a << endl;
}

