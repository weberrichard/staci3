#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../HydraulicSolver.h"

using namespace std;
using namespace Eigen;

vector<string> separate_lines(string line);
vector<string> separate_lines_semicolon(string line);

int main(int argc, char* argv[])
{
   vector<int> got_material;
   vector<int> number_of_pipes;

   // Name of containing folder of staci file
   string caseFolder = "../../Networks/Sopron/";

   vector<string> everyCase;
   everyCase.push_back("buk_2021");
   /*everyCase.push_back("villasor");
   everyCase.push_back("ferto");
   everyCase.push_back("sanchegy");
   everyCase.push_back("buk");
   everyCase.push_back("lovo");
   everyCase.push_back("nagycenk");
   everyCase.push_back("vashegy");
   everyCase.push_back("varis");
   everyCase.push_back("becsidomb");
   everyCase.push_back("tomalom");
   everyCase.push_back("szakov");
   everyCase.push_back("kohegy");
   everyCase.push_back("harka");
   everyCase.push_back("pozsonyiut");
   everyCase.push_back("sopronkovesd");
   everyCase.push_back("dudlesz");
   everyCase.push_back("ivan");
   everyCase.push_back("agyagosszergeny");
   everyCase.push_back("kofejto");
   everyCase.push_back("simasag");
   everyCase.push_back("acsad");
   everyCase.push_back("csaford");
   everyCase.push_back("nagylozs");
   everyCase.push_back("balf");
   everyCase.push_back("csapod");
   everyCase.push_back("und");
   everyCase.push_back("rojtokmuzsaj");
   everyCase.push_back("brennberg");
   everyCase.push_back("pusztacsalad");
   everyCase.push_back("kutyahegy");
   everyCase.push_back("nyarliget");
   everyCase.push_back("meszlen");
   everyCase.push_back("fertoujlak");
   everyCase.push_back("gorbehalom");
   everyCase.push_back("tozeggyarmajor");
   everyCase.push_back("ebergoc");
   everyCase.push_back("csillahegy");
   everyCase.push_back("jerevan");
   everyCase.push_back("gloriette");
   everyCase.push_back("alomhegy");
   everyCase.push_back("ohermes");
   everyCase.push_back("ujhermes");*/

   int nCases = everyCase.size();
   cout << endl << "   CASES\n***********\n";
   for(int i=0; i<nCases; i++)
      cout << i+1 << "  " << everyCase[i] << endl;
   
   makeDirectory("Network Data");

   // loading pipe material data from csv
   vector<string> raw_data = readVectorString("Pipe_mat_20210223.csv");
   vector<string> pipe_id;
   vector<string> pipe_material;
   vector<int> pipe_visited;
   for(int i=0; i<raw_data.size(); i++)
   {
      vector<string> line = separate_lines(raw_data[i]);
      vector<string> ids = separate_lines_semicolon(line[3]);
      for(int j=0; j<ids.size(); j++)
      {
         pipe_id.push_back(ids[j]);
         string mat = line[4];
         if(mat == "kpe")
            pipe_material.push_back("hdpe");
         else if(mat == "ov")
            pipe_material.push_back("ci");
         else if(mat == "a")
            pipe_material.push_back("s");
         else if(mat == "ko")
            pipe_material.push_back("ss");
         else if(mat == "gov")
            pipe_material.push_back("sgci");
         else if(mat == "olom")
            pipe_material.push_back("lead");
         else
            pipe_material.push_back(mat);
         pipe_visited.push_back(0);
      }
   }

   // loading each WDN
   for(int i=0; i<nCases; i++)
   {
      printf("\n[*] %15s\n", everyCase[i].c_str());

      string caseName = everyCase[i];

      HydraulicSolver *wds = new HydraulicSolver(caseFolder + caseName + ".inp");

      for(int j=0; j<wds->pipeIndex.size(); j++)
      {  
         bool got_it=false;
         int k=0;
         while(!got_it && k<pipe_id.size())
         {  
            int idx = wds->pipeIndex[j];
            if("PIPE_" + pipe_id[k] == wds->edges[idx]->name)
            {
               got_it = true;
               wds->edges[idx]->setStringProperty("material",pipe_material[k]);
               pipe_visited[k]++;
            }
            else
            {
               k++;
            }
         }
      }

      // checking how much material we got
      int counter=0;
      for(int j=0; j<wds->pipeIndex.size(); j++)
      {
         if(wds->edges[wds->pipeIndex[j]]->getStringProperty("material") != "unkown")
         {
            counter++;
         }
      }
      got_material.push_back(counter);
      number_of_pipes.push_back(wds->pipeIndex.size());

      string newFileName = "Network Data/" + everyCase[i] + "_mat.inp";
      wds->saveSystem(newFileName);
   }

   // printing some statistics about wdns
   for(int i=0; i<nCases; i++)
   {
      double per = (double)got_material[i] / (double)number_of_pipes[i] *100.;
      printf("%-15s: %5i/%-5i (%5.1f)\n", everyCase[i].c_str(), got_material[i], number_of_pipes[i], per);
   }

   // printing some stats about material data
   cout << "\n number of all pipe material: " << pipe_id.size() << endl;
   int counter=0;
   int visit_max=0;
   for(int i=0; i<pipe_id.size(); i++)
   {
      if(pipe_visited[i]!=0)
      {
         counter++;
      }
      if(visit_max<pipe_visited[i])
      {
         visit_max=pipe_visited[i];
      }
   }
   cout << " number of all visited pipe material: " << counter << endl;
   cout << " number of max visited pipe material: " << visit_max << endl;

   cout << endl << endl;
   return 0;
}

//--------------------------------------------------
vector<string> separate_lines(string line)
{
   string s="";
   vector<string> sv;
   for(string::iterator i=line.begin(); i!=line.end(); i++)
   {
      if(*i!=' ' && *i!='\t' && *i!='\n' && *i!='\r' && *i!=',')
         s += *i;
      else if(s.length()>0)
      {
         sv.push_back(s);
         s="";
      }
   }
   if(s!="")
      sv.push_back(s);
   return sv;
}

//--------------------------------------------------
vector<string> separate_lines_semicolon(string line)
{
   string s="";
   vector<string> sv;
   for(string::iterator i=line.begin(); i!=line.end(); i++)
   {
      if(*i!=';')
         s += *i;
      else if(s.length()>0)
      {
         sv.push_back(s);
         s="";
      }
   }
   if(s!="")
      sv.push_back(s);
   return sv;
}
