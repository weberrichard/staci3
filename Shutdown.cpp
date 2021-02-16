#include "Shutdown.h"

//--------------------------------------------------------------
Shutdown::Shutdown(string spr_filename) : HydraulicSolver(spr_filename)
{
  buildSegmentGraph();
}
Shutdown::~Shutdown(){}

//--------------------------------------------------------------
void Shutdown::buildSegmentGraph()
{ 
  // clearing node and edge segment variable
  for(int i=0; i<numberNodes; i++)
  {
    nodes[i]->segment = -1;
  }
  for(int i=0; i<numberEdges; i++)
  {
    edges[i]->segment = -1;
  }

  updateEdgeVector();

  // creating the edge vector of the segments
  segmentVector = segmenter(edgeVector);
  numberSegment = segmentVector.size();

  // Setting the segment variables of nodes and edges
  for(int i=0; i<segmentVector.size(); i++)
  {
    for(int j=0; j<segmentVector[i].size(); j+=2)
    {
      nodes[abs(segmentVector[i][j])]->segment = i;
      nodes[abs(segmentVector[i][j+1])]->segment = i;
      for(int k=0; k<edges.size(); k++)
      {
        int typeCode = edges[k]->typeCode;
        if(typeCode != -1 && typeCode != -2 && typeCode != 9) // -1: pp, -2: pool, 9: iso
        {
          if(abs(segmentVector[i][j]) == edges[k]->startNodeIndex && abs(segmentVector[i][j+1]) == edges[k]->endNodeIndex)
            edges[k]->segment = i;
        }
        if(typeCode == -1 || typeCode == -2) // -1: pp, -2: pool
        {
          if(abs(segmentVector[i][j]) == edges[k]->startNodeIndex || abs(segmentVector[i][j+1]) == edges[k]->startNodeIndex)
            edges[k]->segment = i;
        }
      }
    }
  }

  // Handling the nodes connecting only ISO valves
  int counter = numberSegment;
  for(int i=0; i<numberNodes; i++)
  {
    if(nodes[i]->segment == -1)
    { 
      nodes[i]->segment = counter;
      for(unsigned int j=0; j<nodes[i]->edgeIn.size(); j++)
      {
        edges[nodes[i]->edgeIn[j]]->segment = counter;
      }
      for(unsigned int j=0; j<nodes[i]->edgeOut.size(); j++)
      {
        edges[nodes[i]->edgeOut[j]]->segment = counter;
      }
      counter++;
    }
  }
  numberSegment = counter;

  // Creating the edge vector of the segment graph
  segmentEdgeVector.clear();
  for(unsigned int i=0; i<edges.size(); i++)
  {
    if(edges[i]->typeCode == 9) // 9: iso
    {
      int idx_from, idx_to;
      idx_from = edges[i]->startNodeIndex;
      idx_to = edges[i]->endNodeIndex;
      segmentEdgeVector.push_back(nodes[idx_from]->segment);
      segmentEdgeVector.push_back(nodes[idx_to]->segment);

      // Filling up the Valves which segments are they connecting
      edges[i]->setIntProperty("startSegment",nodes[idx_from]->segment);
      edges[i]->setIntProperty("endSegment",nodes[idx_to]->segment);
    }
  }

  segmentRank = calculateRank(segmentEdgeVector);

  // calculating how much pipe lengths are in each segment
  vector<double> Li(numberSegment,0.); // in meter
  vector<double> li(numberSegment,0.); // dimensionless
  relativePipeLength.resize(numberSegment,0.);
  absolutePipeLength.resize(numberSegment,0.);
  double sumL=0.; // in meter

  for(int i=0; i<edges.size(); i++)
  {
      if(edges[i]->typeCode == 1 || edges[i]->typeCode == 0) // pipe, pipeCV
      {
        double L = edges[i]->getDoubleProperty("length");
        Li[edges[i]->segment] += L;
        sumL += L;
      }
  }
  for(int i=0; i<numberSegment; i++)
  {
    li[i] = Li[i] / sumL;
    relativePipeLength[i] = li[i];
    absolutePipeLength[i] = Li[i];
  }
}

//--------------------------------------------------------------
void Shutdown::openEverything()
{
  for(int i=0; i<numberEdges; i++)
  {
    edges[i]->status = 1;
  }
  for(int i=0; i<numberNodes; i++)
  {
    nodes[i]->status = 1;
  }
}

//--------------------------------------------------------------
vector<int> Shutdown::closeSegment(int segmentToClose)
{
  for(int i=0; i<edges.size(); i++)
  {
    if(edges[i]->segment == segmentToClose)
    {
      edges[i]->status = 0;
    }
    if(edges[i]->typeCode == 9) // 9: iso
    {
      if(edges[i]->getIntProperty("startSegment") == segmentToClose || edges[i]->getIntProperty("endSegment") == segmentToClose)
      {
        edges[i]->status = 0;
      }
    }
  }

  for(int i=0; i<nodes.size(); i++)
  {
    if(nodes[i]->segment == segmentToClose)
    {
      nodes[i]->status = 0;
    }
  }

  //vector<int> closedSegment = closeDisconnectedParts();
  vector<int> closedSegment = closeDisconnectedSegments(segmentToClose);

  return closedSegment;
}

//--------------------------------------------------------------
vector<int> Shutdown::closeDisconnectedSegments(int segmentToClose)
{
  vector<int> sac; // segments closed additionally
  

  // building the segment edge vector without the closed segment 
  vector<int> ev; // new segment edge vector
  for(int i=0; i<segmentEdgeVector.size(); i+=2)
  {
    if(segmentEdgeVector[i] != segmentToClose && segmentEdgeVector[i+1] != segmentToClose)
    {
      ev.push_back(segmentEdgeVector[i]);
      ev.push_back(segmentEdgeVector[i+1]);
    }
    if(segmentEdgeVector[i+1] == segmentToClose && segmentRank[segmentEdgeVector[i]] == 1)
    {
      sac.push_back(segmentEdgeVector[i]);
    }
    if(segmentEdgeVector[i] == segmentToClose && segmentRank[segmentEdgeVector[i+1]] == 1)
    {
      sac.push_back(segmentEdgeVector[i+1]);
    }
  }

  // collecting the pressured segments
  vector<int> pi;
  for(int i=0; i<poolIndex.size(); i++)
  {
    pi.push_back(edges[poolIndex[i]]->segment);
  }
  for(int i=0; i<presIndex.size(); i++)
  {
    pi.push_back(edges[presIndex[i]]->segment);
  }

  // there might be more
  pi = unique(pi);

  vector<int> sac2;
  // clearing pressured segemnts from additionally closed segments
  for(int i=0; i<sac.size(); i++)
  {
    bool ispres=false;
    for(int j=0; j<pi.size(); j++)
    {
      if(sac[i] == pi[j])
      {
        ispres = true;
      }
    }
    if(!ispres)
    {
      sac2.push_back(sac[i]);
    }
  }

  vector<vector<int> > ss = segmenter(ev); // segments of new segment vector

  // searching the segments that are not connecting to pressure source
  for(int i=0; i<ss.size(); i++)
  {
    bool ispres = false;
    for(int j=0; j<ss[i].size(); j++)
    {
      for(int k=0; k<pi.size(); k++)
      {
        if(ss[i][j] == pi[k])
        {
          ispres = true;
        }
      }
    }
    if(!ispres)
    {
      for(int j=0; j<ss[i].size(); j++)
      {
        sac2.push_back(ss[i][j]);
      }
    }
  }

  sac2 = unique(sac2);

  // closing nodes and edges
  for(int i=0; i<numberNodes; i++)
  {
    for(int j=0; j<sac2.size(); j++)
    {
      if(nodes[i]->segment == sac2[j])
      {
        nodes[i]->status = 0;
        break;
      }
    }
  }

  for(int i=0; i<numberEdges; i++)
  {
    for(int j=0; j<sac2.size(); j++)
    {
      if(edges[i]->segment == sac2[j])
      {
        edges[i]->status = 0;
        break;
      }
      if(edges[i]->typeCode == 9)
      {
        if(edges[i]->getIntProperty("startSegment") == sac2[j] || edges[i]->getIntProperty("endSegment") == sac2[j])
        {
          edges[i]->status = 0;
        }
      }
    }
  }

  return sac2;
}

//--------------------------------------------------------------
vector<int> Shutdown::closeDisconnectedParts()
{
  // build up edge vector that contains the isolation valves as well
  vector<int> edgeVectorWithValves;
  vector<vector<int> > tmp(2);
  vector<int> pressurePrescribed;
  for(int i=0; i<edges.size(); i++)
  {
    if(edges[i]->status > 0)
    {
      int typeCode = edges[i]->typeCode;
      if(typeCode != -1 && typeCode != -2)  // -2: pp, -1: pool
      {
        edgeVectorWithValves.push_back(edges[i]->startNodeIndex);
        edgeVectorWithValves.push_back(edges[i]->endNodeIndex);
      }
      else
      {
        pressurePrescribed.push_back(edges[i]->startNodeIndex);
      }
    }
  }

  // checking the integrity
  vector<vector<int> > segmentVectorWithValves = segmenter(edgeVectorWithValves);

  vector<int> closedSegment;
  // if we have more than one "part", checking which ones have pressure prescribed
  if(segmentVectorWithValves.size()>1)
  {
    vector<bool> isPressurePrescribed(segmentVectorWithValves.size(),false);
    for(int k=0; k<pressurePrescribed.size(); k++)
    {
      for(int i=0; i<segmentVectorWithValves.size(); i++)
      {
        for(int j=0; j<segmentVectorWithValves[i].size(); j++)
        {
          if(segmentVectorWithValves[i][j] == pressurePrescribed[k])
          {
            isPressurePrescribed[i] = true;
          }
        }
      }
    }

    // closing the nodes those are not connected to pressure node
    for(int i=0; i<isPressurePrescribed.size(); i++)
    {
      if(!isPressurePrescribed[i])
      {
        closedSegment.push_back(nodes[segmentVectorWithValves[i][0]]->segment);
        for(int j=0; j<segmentVectorWithValves[i].size(); j++)
        {
          nodes[segmentVectorWithValves[i][j]]->status = 0;
        }
      }
    }

    // closing the edges as well
    for(int i=0; i<edges.size(); i++)
    {
      int idx_from = edges[i]->startNodeIndex;
      if(nodes[idx_from]->status < 1)
        edges[i]->status = 0;
    }
  }

  return closedSegment;
}

//--------------------------------------------------------------
// Building up the edge vector from pipes, pumps and valves (not ISO)
void Shutdown::updateEdgeVector()
{
  // basic edge vector without ISO valves
  edgeVector.clear();
  // full edge vector
  for(int i=0; i<edges.size(); i++)
  {
    if(edges[i]->status > 0)
    {
      int typeCode = edges[i]->typeCode;
      // pressurepoint: -1 | pool: -2 | iso: 9
      if(typeCode != -1 && typeCode != -2 && typeCode != 9)
      {
        edgeVector.push_back(edges[i]->startNodeIndex);
        edgeVector.push_back(edges[i]->endNodeIndex);
      }
    }
  }
}

//--------------------------------------------------------------
vector<int> Shutdown::findConnectionError(vector<int> connectingNodes)
{ 
  vector<int> connectionErrors;
  vector<int> rank(nodes.size());
  for(int i=0; i<edges.size(); i++)
  { 
    int typeCode = edges[i]->typeCode;
    if(typeCode != -1 && typeCode != -2 && typeCode != 9) // -2: pool, -1: pp, 9: iso
    {
      rank[edges[i]->startNodeIndex]++;
      rank[edges[i]->endNodeIndex]++;
    }
    if(typeCode == -1 || typeCode == -2) // -2: pp, -1: pool
    { 
      rank[edges[i]->startNodeIndex]++;
    }
  }

  for(int i=0; i<connectingNodes.size(); i++)
  {
    if(rank[connectingNodes[i]]!=2 && rank[connectingNodes[i]]!=0)
    {
      connectionErrors.push_back(connectingNodes[i]);
    }
  }

  if(connectionErrors.size()>0)
    cout << endl << "[*] Cleared connection nodes (with  rank): ";
  for(int i=0; i<connectionErrors.size(); i++)
  {
    cout << nodes[connectionErrors[i]]->name << "(" << rank[connectionErrors[i]] << ")" << ", ";
  }

  return connectionErrors;
}

//--------------------------------------------------------------
vector<int> Shutdown::calculateRank(const vector<int> &ev)
{
  vector<int> rank;
  if(ev.size() != 0)
  {
    int max_ev = max(ev);
    rank.resize(max_ev+1);
    for(int i=0; i<ev.size(); i++)
    {
      rank[ev[i]]++;
    }
  }
  else
  {
    rank.resize(1);
    rank[0] = 0; 
  }

  return rank;
}


//--------------------------------------------------------------
//vector<string> Shutdown::shutdownPlan(string pipeID)
//{
//  int idx = edgeIDtoIndex(pipeID);
//
//  vector<string> valvesIdx;
//  if(idx == -1)
//  {
//    cout << endl << "Shutdown::shutdownPlan(pipeID) : pipeID (" << pipeID << ") not found" << endl;
//  }
//  else
//  {
//    int segment = edges[idx]->segment;
//    valvesIdx = closeSegment(segment);
//  }
//
//  return valves;
//}
// Clearing the inner isolation valves
// this was for ShutdownPlan, shoud be updated
/*for(int i=0; i<closedValves.size(); i++)
{
  for(int j=0; j<edges.size(); j++)
  {
    if(edges[j]->typeCode == 9) // 9: iso
    {
      if(edges[j]->name == closedValves[i])
      { 
        for(int k=0; k<closedSegment.size(); k++)
        {
          if(edges[j]->getIntProperty("startSegment") == closedSegment[k] || edges[j]->getIntProperty("endSegment") == closedSegment[k])
          {
            closedValves[i] = "no";
          }
        }
      }
    }
  }
}

for(int i=0; i<closedValves.size(); i++)
  if(closedValves[i] == "no")
    closedValves.erase(closedValves.begin() + i);*/
