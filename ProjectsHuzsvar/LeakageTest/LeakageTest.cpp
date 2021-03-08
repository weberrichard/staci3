#include <fstream>
#include <sstream>
#include "../../Leakage.h"


using namespace std;
using namespace Eigen;
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//typedef pair<int, int> Edge;

Leakage *wds;


int main(int argc, char *argv[])
{
	string case_folder = "../../Networks/Sopron/";
    string case_name = "tomalom";

    cout << "[*]Network: " << case_folder << case_name << ".inp" << endl;
    wds = new Leakage(case_folder +case_name + ".inp");
    wds->calculateLeakage();
}
