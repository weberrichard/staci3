#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../Vulnerability.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
   // Name of containing folder of staci file
   string caseFolder = "../../Networks/Sopron/";
   //string caseFolder = "../../Networks/Sopron/";

   vector<string> everyCase;
   //everyCase.push_back("villasor_mat_year");
   //everyCase.push_back("ferto_mat_year");
   //everyCase.push_back("sanchegy_mat_year");
   //everyCase.push_back("buk_mat_year");
   //everyCase.push_back("lovo_mat_year");
   //everyCase.push_back("nagycenk_mat_year");
   //everyCase.push_back("vashegy_mat_year");
   //everyCase.push_back("varis_mat_year");
   //everyCase.push_back("becsidomb_mat_year");
   //everyCase.push_back("tomalom_mat_year");
   //everyCase.push_back("szakov_mat_year");
   //everyCase.push_back("kohegy_mat_year");
   //everyCase.push_back("harka_mat_year");
   //everyCase.push_back("pozsonyiut_mat_year");
   //everyCase.push_back("sopronkovesd_mat_year");
   //everyCase.push_back("dudlesz_mat_year");
   //everyCase.push_back("ivan_mat_year");
   //everyCase.push_back("agyagosszergeny_mat_year");
   //everyCase.push_back("kofejto_mat_year");
   //everyCase.push_back("simasag_mat_year");
   //everyCase.push_back("acsad_mat_year");
   //everyCase.push_back("csaford_mat_year");
   everyCase.push_back("nagylozs_mat_year");
   //everyCase.push_back("balf_mat_year");
   //everyCase.push_back("csapod_mat_year");
   //everyCase.push_back("und_mat_year");
   //everyCase.push_back("rojtokmuzsaj_mat_year");
   //everyCase.push_back("brennberg_mat_year");
   //everyCase.push_back("pusztacsalad_mat_year");
   //everyCase.push_back("kutyahegy_mat_year");
   //everyCase.push_back("nyarliget_mat_year");
   //everyCase.push_back("meszlen_mat_year");
   //everyCase.push_back("fertoujlak_mat_year");
   //everyCase.push_back("gorbehalom_mat_year");
   //everyCase.push_back("tozeggyarmajor_mat_year");
   //everyCase.push_back("ebergoc_mat_year");
   //everyCase.push_back("csillahegy_mat_year");
   //everyCase.push_back("jerevan_mat_year");
   //everyCase.push_back("gloriette_mat_year");
   //everyCase.push_back("alomhegy_mat_year");
   //everyCase.push_back("ohermes_mat_year");
   //everyCase.push_back("ujhermes_mat_year");

   int nCases = everyCase.size();
   cout << endl << "   CASES\n***********\n";
   for(int i=0; i<nCases; i++)
      cout << i+1 << "  " << everyCase[i] << endl;
   
   srand( (unsigned)time(NULL) );

   // for writing to files
   ofstream wFile;

   vector<vector<double> > everyLocalGamma(nCases);
   for(int i=0; i<nCases; i++)
   {
      printf("\n[*] %15s\n", everyCase[i].c_str());

      string caseName = everyCase[i];
      Vulnerability *wds = new Vulnerability(caseFolder + caseName + ".inp");

      cout << "load ok " << endl;

      wds->calculateVulnerability();

      double orig_gamma = wds->globalGamma;

      cout << endl << "orig vulner ok " << endl;

      wds->calculateBackupVulnerability();

      vector<double> backup_gamma = wds->backupGamma;

      cout << endl << "orig: " << orig_gamma << endl;
      for(int j=0; j<backup_gamma.size(); j++)
      {
         cout << j << " " << backup_gamma[j]-orig_gamma << endl;
      }

   }
   
   cout << endl << endl;
   return 0;
}