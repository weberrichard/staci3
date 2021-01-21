#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

#include "../../Shutdown.h"
#include "../../HydraulicSolver.h"

using namespace std;
using namespace Eigen;

int main()
{
  // Name of containing folder of staci file
  string caseFolder = "../../Networks/Sopron/";

  vector<string> everyCase;
  everyCase.push_back("tomalom");
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

  // topological stuff
  vector<double> n1n3p(nCases,0.0);
  vector<double> L1L3p(nCases,0.0);
  vector<double> V1V3p(nCases,0.0);
  // one hydraulic simulation stuf
  vector<double> globalAlfa(nCases,0.0);
  vector<double> globalAlfa2(nCases,0.0);
  // multiple hydraulic simulation stuff
  vector<double> globalBeta(nCases,0.0);
  vector<double> globalGamma(nCases,0.0);
  // basic parameters
  vector<double> everyLength(nCases,0.0);
  vector<double> everyCons(nCases,0.0);
  vector<int> everyNodeSize(nCases,0.0);
  vector<int> everyNumberSegment(nCases,0);
  vector<int> everyNumberISOValve(nCases,0);
  vector<int> everyNumberPool(nCases,0);
  vector<int> everyNumberPres(nCases,0);
  vector<vector<int> > everyRank(nCases);
  vector<vector<double> > everySegmentLength(nCases);
  vector<vector<double> > everySegmentVolume(nCases);
  vector<vector<double> > everyAlfaLength(nCases);
  vector<vector<double> > everyAlfaRelLength(nCases);
  vector<vector<double> > everyBetaLength(nCases);
  vector<vector<double> > everyBetaRelLength(nCases);
  vector<vector<double> > everyGammaLength(nCases);
  vector<vector<double> > everyGammaRelLength(nCases);
  // for graph stuff
  vector<double> everyClustering(nCases);
  vector<double> everyPathLength(nCases);
  vector<int> everyDiameter(nCases);
  vector<int> everyNumberArticulation(nCases);

  int absMaxRank=0;

  //========================
  //       MAIN CYCLE       
  //========================
  cout << endl << "   PROGRESS\n**************\n";
  for(int i=0; i<nCases; i++)
  {
    string caseName = everyCase[i];
    Shutdown *sd = new Shutdown(caseFolder + caseName + ".inp");
    int nSegments = sd->getNumberSegment();
    everyNumberSegment[i] = nSegments;
    for(int j=0; j<sd->edges.size(); j++)
    {
      if(sd->edges[j]->typeCode == 9)
        everyNumberISOValve[i]++;
    }

    for(int j=0; j<sd->edges.size(); j++)
    {
      if(sd->edges[j]->typeCode == -1)
      {
        everyNumberPool[i]++;
      }
      else if(sd->edges[j]->typeCode == -2)
      {
        everyNumberPres[i]++;
      }
    }

    // CALCULATING TOPOLOGICAL STUFF
    //*******************************
    vector<int> rank(nSegments,0);
    for(int j=0; j<sd->edges.size(); j++)
    {
      if(sd->edges[j]->typeCode == 9) // iso
      {
        rank[sd->edges[j]->getIntProperty("startSegment")]++;
        rank[sd->edges[j]->getIntProperty("endSegment")]++;
      }
    }

    // find max rank
    int maxRank=0;
    for(int j=0; j<rank.size(); j++)
    {
      if(rank[j]>maxRank)
        maxRank = rank[j];
      if(maxRank>absMaxRank)
        absMaxRank = maxRank;
    }

    // resize vector vector
    for(int j=0; j<i+1; j++)
    {
      everyRank[j].resize(absMaxRank+1,0);
      everySegmentLength[j].resize(absMaxRank+1,0.);
      everySegmentVolume[j].resize(absMaxRank+1,0.);
    }

    // filling ranks in everyRank
    for(int j=0; j<rank.size(); j++)
    {
      everyRank[i][rank[j]]++;
    }
    // filling everySegmentLength
    for(int j=0; j<sd->edges.size(); j++)
    {
      int tc = sd->edges[j]->typeCode;
      if(tc == 0 || tc ==1)
      {
        everySegmentLength[i][rank[sd->edges[j]->segment]] += sd->edges[j]->getDoubleProperty("length");
        everySegmentVolume[i][rank[sd->edges[j]->segment]] += sd->edges[j]->getDoubleProperty("volume");
      }
    }

    int n1=0, n3p=0;
    for(int j=0; j<nSegments; j++)
    {
      if(rank[j]==1)
        n1++;
      if(rank[j]>=3)
        n3p++;
    }
    n1n3p[i] = (double)n1/(double)n3p;

    double L1=0,L3p=0;
    double V1=0,V3p=0;
    for(int j=0; j<sd->edges.size(); j++)
    { 
      int typeCode = sd->edges[j]->typeCode; 
      if(typeCode == 0 || typeCode == 1)
      {
        if(rank[sd->edges[j]->segment] == 1)
        {
          L1 += sd->edges[j]->getDoubleProperty("length");
          V1 += sd->edges[j]->getDoubleProperty("volume");
        }
        if(rank[sd->edges[j]->segment] >= 3)
        {
          L3p += sd->edges[j]->getDoubleProperty("length");
          V3p += sd->edges[j]->getDoubleProperty("volume");
        }
      }
    }
    L1L3p[i] = L1/L3p;
    V1V3p[i] = V1/V3p;

    // CALCULATING ALFA | ONE HYDRAULIC SIMULATION
    //*********************************************
    HydraulicSolver *hs = new HydraulicSolver(caseFolder + caseName + ".inp");

    hs->solveSystem();

    everyNodeSize[i] = hs->nodes.size();

    vector<double> localAlfa(nSegments,0.0);
    vector<double> segmentLength(nSegments,0.0);
    for(int j=0; j<hs->edges.size(); j++)
    {
      int typeCode = hs->edges[j]->typeCode;
      // Calculating the lengths of the segments
      if(typeCode == 1 || typeCode == 0) // pipe, pipeCV
      {
        segmentLength[sd->edges[j]->segment] += hs->edges[j]->getDoubleProperty("length");
      }
      // Calculating the volume flow rate in the segments
      double vf = hs->edges[j]->volumeFlowRate;
      if(typeCode == 9) // iso
      {
        if(vf>0.)
          localAlfa[sd->edges[j]->getIntProperty("endSegment")] += vf;
        if(vf<0.)
          localAlfa[sd->edges[j]->getIntProperty("startSegment")] += -vf;
      }
      if(typeCode == -1 || typeCode == -2) // pool and pres
      { 
        //cout << endl << hs->edges[j]->name << "  vf: " << hs->edges[j]->volumeFlowRate << endl;
        //if(vf>0.)
        //  localAlfa[sd->edges[j]->segment] += vf;
        if(vf<0.)
          localAlfa[sd->edges[j]->segment] += -vf;
      }
    }

    // reranging alfa to [0,1] by deviding with the sum inflow
    double sumDemand=0.0;
    for(int j=0; j<hs->nodes.size(); j++)
        sumDemand += hs->nodes[j]->demand;
    for(int j=0; j<localAlfa.size(); j++)
      localAlfa[j] /= sumDemand;


    double sumLength=0.0;
    for(int j=0; j<nSegments; j++)
      sumLength += segmentLength[j];
    everyLength[i] = sumLength;

    for(int j=0; j<sd->nodes.size(); j++)
      everyCons[i] += sd->nodes[j]->demand;

    for(int j=0; j<nSegments; j++)
      globalAlfa[i] += localAlfa[j]*segmentLength[j];
    globalAlfa[i] /= sumLength;


    everyAlfaLength[i].resize(nSegments,0.0);
    everyAlfaRelLength[i].resize(nSegments,0.0);
    for(int j=0; j<nSegments; j++)
    {
      everyAlfaLength[i][j] = localAlfa[j]*segmentLength[j];
      everyAlfaRelLength[i][j] = localAlfa[j] * segmentLength[j] / sumLength;
    }


    vector<double> localAlfa2(nSegments,0.0);
    vector<double> segmentLength2(nSegments, 0.0);
    double sumLength2=0.0;
    double sumDemand2=0.0;
    for(int j=0; j<hs->nodes.size(); j++)
    {
      localAlfa2[sd->nodes[j]->segment] += hs->nodes[j]->demand;
      sumDemand2 += sd->nodes[j]->demand;
    }


    for(int j=0; j<hs->edges.size(); j++)
    {
      int tc = hs->edges[j]->typeCode;
      double vf = hs->edges[j]->volumeFlowRate;
      if(tc == 9) // iso
      {
        if(vf > 0.)
          localAlfa2[sd->edges[j]->getIntProperty("startSegment")] += vf;
        else
          localAlfa2[sd->edges[j]->getIntProperty("endSegment")] += -vf;
      }
      if(tc == -1 || tc == -2) // pool and pres
      {
        if(vf > 0.)
          localAlfa2[sd->edges[j]->segment] += vf;
      }
      if(tc == 1 || tc == 0)
      {
        double l = sd->edges[j]->getDoubleProperty("length");
        segmentLength2[sd->edges[j]->segment] += l;
        sumLength2 += l;
      }
    }


    for(int j=0; j<nSegments; j++)
      localAlfa2[j] /= sumDemand2;

    for(int j=0; j<nSegments; j++)
    {
      globalAlfa2[i] += segmentLength2[j]*localAlfa2[j];
    }
    globalAlfa2[i] /= sumLength2;


    //for(int j=0; j<sd->edges.size(); j++)
    //{
    //  hs->edges[j]->setEdgeDoubleProperty("userOutput",localAlfa[sd->edges[j]->segment]);
    //}
    //for(int j=0; j<sd->nodes.size(); j++)
    //{
    //  hs->nodes[j]->setProperty("userOutput",localAlfa[sd->nodes[j]->segment]);
    //}

    //hs->saveResult("volumeFlowRate","Pipe");
    //hs->saveResult("volumeFlowRate","Pool");
    //hs->saveResult("volumeFlowRate","PressurePoint");
    //hs->saveResult("head","Node");
    //hs->saveResult("volumeFlowRate","Valve");
    //sd->saveResult("segment","Pipe");
    //sd->saveResult("segment","All");

    // CALCULATING BETA NEW
    //*******************************
    // the demands are depending on the pressure

    printf("%15s\n", everyCase[i].c_str());


    // avoiding the close of the only pressure source
    int nPres = everyNumberPres[i], nPool = everyNumberPool[i];
    int doNotClose = -2;
    if(nPres == 1 && nPool == 0)
      doNotClose = sd->edges[sd->presIndex[0]]->segment;
    else if(nPres == 0 && nPool == 1)
      doNotClose = sd->edges[sd->poolIndex[0]]->segment;


    vector<double> localBeta(nSegments,sumDemand);
    vector<double> localGamma(nSegments,0.0);
    vector<double> segmentDemand(nSegments,0.0);
    for(int j=0; j<nSegments; j++)
    {
      if(j != doNotClose)
      {
        sd->openEverything();
        sd->closeSegment(j);

        for(int k=0; k<hs->edges.size(); k++)
        {
          hs->edges[k]->status = sd->edges[k]->status;
        }
        for(int k=0; k<hs->nodes.size(); k++)
        {
          hs->nodes[k]->status = sd->nodes[k]->status;
        }

        hs->isPressureDemand = false;
        hs->solveSystem();

        for(int k=0; k<hs->nodes.size(); k++)
        {
          localBeta[j] -= hs->nodes[k]->getProperty("consumption");
        }
        //localBeta[j] /= sumDemand;

        hs->isPressureDemand = true;
        hs->solveSystem();

        for(int k=0; k<hs->nodes.size(); k++)
        {
          localGamma[j] += hs->nodes[k]->getProperty("consumption");
        }
        localGamma[j] = sumDemand - localGamma[j];

        //localGamma[j] /= sumDemand;
        //localGamma[j] = 1. - localGamma[j];
      }
      printf(" i / nCase: %3i / %-3i  |  j / nSeg: %3i / %-3i", i, nCases, j, nSegments);
      cout << endl;
    }

    everyBetaLength[i].resize(nSegments,0.0);
    everyBetaRelLength[i].resize(nSegments,0.0);
    for(int j=0; j<nSegments; j++)
    {
      everyBetaLength[i][j] = localBeta[j] * segmentLength[j];
      everyBetaRelLength[i][j] = localBeta[j] * segmentLength[j] / sumLength;
    }

    everyGammaLength[i].resize(nSegments,0.0);
    everyGammaRelLength[i].resize(nSegments,0.0);
    for(int j=0; j<nSegments; j++)
    {
      everyGammaLength[i][j] = localGamma[j] * segmentLength[j];
      everyGammaRelLength[i][j] = localGamma[j]/sumDemand * segmentLength[j]/sumLength ;
    }

    for(int k=0; k<nSegments; k++)
    {
      globalBeta[i] += localBeta[k] * segmentLength[k];
      globalGamma[i] += localGamma[k] * segmentLength[k];
    }
    globalBeta[i] /= sumLength;
    globalGamma[i] /= (sumLength*sumDemand);

    for(int k=0; k<hs->nodes.size(); k++)
    {
      segmentDemand[sd->nodes[k]->segment] += hs->nodes[k]->demand;
    }

    for(int j=0; j<sd->nodes.size(); j++)
    {
      //hs->nodes[j]->setProperty("userOutput",sd->nodes[j]->segment);
      hs->nodes[j]->setProperty("userOutput",everyGammaRelLength[i][sd->nodes[j]->segment]);
      //hs->nodes[j]->setProperty("userOutput",segmentDemand[sd->nodes[j]->segment]);
    }
    hs->saveResult("userOutput","Node");
    //sd->saveResult("segment","Pipe");
    //sd->saveResult("segment","All");

    // saving the segment edge vector for matlab postproc
    mkdir("Network Data",0777);
    mkdir(("Network Data/" + caseName).c_str(),0777);
    remove(("Network Data/" + caseName + "/EdgeVector.txt").c_str());

    vector<int> sev = sd->getSegmentEdgeVector();
    ofstream wfile;
    wfile.open("Network Data/" + caseName + "/EdgeVector.txt");
    for(int i=0; i<sev.size(); i+=2)
      wfile << sev[i] << ',' << sev[i+1] << endl;
    wfile.close();

    // CALCULATING STRUCTURAL GRAPH STUFF
    //everyClustering[i] = globalClustering(sev);
    //everyPathLength[i] = avPathLength(sev);
    //everyDiameter[i] = graphDiameter(sev);
    //everyNumberArticulation[i] = articulationPoints(sev);

    everyClustering[i] = 0.;
    everyPathLength[i] = 0.;
    everyDiameter[i] = 0.;
    everyNumberArticulation[i] = 0.;
  }

  // ********************
  // PRINTING EVERYTHING
  // ********************

  cout << "\n   RESULTS   ";
  cout << "\n*************\n" << endl;
  printf(" %2s | %6s | %6s | %6s | %6s | %5s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s  \n", "i", "n_pres", "n_pool", "n_seg", "n_iso", "nodes", "n1/n3+", "L1/L3+", "V1/V3+", "alfa", "alfa2", "beta", "gamma", "l [km]", "c [cmh]", "apl", "diam", "cluster", "artic");
  cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------M1=[";
  for(int i=0; i<nCases; i++)
  {
    printf("\n %2i , %6i , %6i , %6i , %6i , %5i , %8.4f , %8.4f , %8.4f , %8.4f , %8.4f , %8.4e , %8.4e , %8.4f , %8.4f , %8.4f , %8i , %8.4f , %8i ;", i+1, everyNumberPres[i], everyNumberPool[i], everyNumberSegment[i], everyNumberISOValve[i], everyNodeSize[i], n1n3p[i], L1L3p[i], V1V3p[i], globalAlfa[i], globalAlfa2[i], globalBeta[i], globalGamma[i], everyLength[i]/1000., everyCons[i]*3600., everyPathLength[i], everyDiameter[i], everyClustering[i], everyNumberArticulation[i]);
  }
  printf("];\n");

  cout << "\n    RANKS   ";
  cout << "\n*************\n" << endl;
  printf(" %2s |", "i");
  for(int i=0; i<absMaxRank+1; i++)
    printf(" %4i |", i);
  printf("\n");
  cout << "-----------------------------------------------------------------------------------------------------------M2=[";
  for(int i=0; i<nCases; i++)
  {
    printf("\n %2i ", i+1);
    for(int j=0; j<everyRank[i].size(); j++)
    {
      printf(", %4i ",everyRank[i][j]);
    }
    printf(";");
  }
  printf("];\n");

  cout << "\n WRITING GAMMABAR TO FILE";
  cout << "\n****************************\n" << endl;

  ofstream wFile, wFile2;
  wFile.open("gammaBarLength.txt");
  wFile2.open("gammaBarRelLength.txt");
  for(int i=0; i<nCases; i++)
  {
    wFile << everyGammaLength[i][0];
    wFile2 << everyGammaRelLength[i][0];
    for(int j=0; j<everyGammaLength[i].size(); j++)
    {
      wFile << ", " << everyGammaLength[i][j];
      wFile2 << ", " << everyGammaRelLength[i][j];
    }
    wFile << ";" << endl;
    wFile2 << ";" << endl;
  }
  wFile.close();
  wFile2.close();

  cout << "\n WRITING GAMMABAR TO FILE: OK ";
  cout << "\n********************************\n" << endl;

  cout << endl << endl;
  return 0;
}
