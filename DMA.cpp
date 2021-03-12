#include "DMA.h"

//--------------------------------------------------------------
DMA::DMA(string spr_filename) : HydraulicSolver(spr_filename){}
DMA::~DMA(){}

//--------------------------------------------------------------
void DMA::determineDMAZones()
{
	// clearing node and edge DMAZone variable
	for(int i=0; i<numberNodes; i++)
	{
		nodes[i]->DMAZone = -1;
	}
	for(int i=0; i<numberEdges; i++)
	{
		edges[i]->DMAZone = -1;
	}

	updateEdgeVector();
	// creating the edge vector of the DMA zones
	zoneVector = segmenter(edgeVector);
	numberDMAZones = zoneVector.size();

	// Setting the DMAZone variables of nodes and edges
	for(int i=0; i<zoneVector.size(); i++)
	{
		for(int j=0; j<zoneVector[i].size(); j+=2)
		{
			nodes[abs(zoneVector[i][j])]->DMAZone = i;
			nodes[abs(zoneVector[i][j+1])]->DMAZone = i;
			for(int k=0; k<edges.size(); k++)
			{
				int typeCode = edges[k]->typeCode;
				if(typeCode != -1 && typeCode != -2 && typeCode != 10) // -1: pp, -2: pool, 9: iso
				{
					if(abs(zoneVector[i][j]) == edges[k]->startNodeIndex && abs(zoneVector[i][j+1]) == edges[k]->endNodeIndex)
						edges[k]->DMAZone = i;
				}
				if(typeCode == -1 || typeCode == -2) // -1: pp, -2: pool
				{
					if(abs(zoneVector[i][j]) == edges[k]->startNodeIndex || abs(zoneVector[i][j+1]) == edges[k]->startNodeIndex)
				   	edges[k]->DMAZone = i;
				}
			}
		}
	}

	// Handling the nodes connecting only ISO valves
	int counter = numberDMAZones;
	for(int i=0; i<numberNodes; i++)
	{
		if(nodes[i]->DMAZone == -1)
		{ 
			nodes[i]->DMAZone = counter;
			for(unsigned int j=0; j<nodes[i]->edgeIn.size(); j++)
			{
			  edges[nodes[i]->edgeIn[j]]->DMAZone = counter;
			}
			for(unsigned int j=0; j<nodes[i]->edgeOut.size(); j++)
			{
			  edges[nodes[i]->edgeOut[j]]->DMAZone = counter;
			}
			counter++;
		}
	}
	numberDMAZones = counter;

  // Creating the edge vector of the DMA zones
  zoneEdgeVector.clear();
  for(unsigned int i=0; i<edges.size(); i++)
  {
    if(edges[i]->typeCode == 10) // 10: flow meter
    {
      int idx_from, idx_to;
      idx_from = edges[i]->startNodeIndex;
      idx_to = edges[i]->endNodeIndex;
      zoneEdgeVector.push_back(nodes[idx_from]->DMAZone);
      zoneEdgeVector.push_back(nodes[idx_to]->DMAZone);

      // Filling up the Valves which segments are they connecting
      edges[i]->setIntProperty("startDMAZone",nodes[idx_from]->DMAZone);
      edges[i]->setIntProperty("endDMAZone",nodes[idx_to]->DMAZone);
    }
  }

}


//--------------------------------------------------------------
// Building up the edge vector from pipes, pumps and valves (without flow meter)
void DMA::updateEdgeVector()
{
  // basic edge vector without flow meters valves
  edgeVector.clear();
  // full edge vector
  for(int i=0; i<edges.size(); i++)
  {
    if(edges[i]->status > 0)
    {
      int typeCode = edges[i]->typeCode;
      // pressurepoint: -1 | pool: -2 | fm: 10
      if(typeCode != -1 && typeCode != -2 && typeCode != 10)
      {
        edgeVector.push_back(edges[i]->startNodeIndex);
        edgeVector.push_back(edges[i]->endNodeIndex);
      }
    }
  }
}

//--------------------------------------------------------------
void DMA::addNewFlowMeter(vector<string> flowMeterName, vector<string> pipeName, vector<bool> isStart, double density, vector<double> referenceCrossSection, double volumeFlowRate)
{ 
  vector<int> pipeISO;
  for(unsigned int i=0; i<pipeName.size(); i++)
  {
    int pi = edgeIDtoIndex(pipeName[i]); // pipe index
    pipeISO.push_back(pi);
  }

  addNewFlowMeter(flowMeterName, pipeISO, isStart, density, referenceCrossSection, volumeFlowRate);
}

//--------------------------------------------------------------
void DMA::addNewFlowMeter(vector<string> flowMeterName, vector<int> pipeISO, vector<bool> isStart, double density, vector<double> referenceCrossSection, double volumeFlowRate)
{
  for(unsigned int i=0; i<flowMeterName.size(); i++)
  {
    int pi = pipeISO[i]; // pipe index
    int ns = nodeIDtoIndex(edges[pi]->startNodeName); // start node index
    int ne = nodeIDtoIndex(edges[pi]->startNodeName); // end node index
    double x = nodes[ns]->xPosition + nodes[ne]->xPosition; // coordinates of additional node
    double y = nodes[ns]->yPosition + nodes[ne]->yPosition;
    string newNode = flowMeterName[i] + "_N";
    if(isStart[i])
    {
      string sn = edges[pi]->startNodeName;
      edges[pi]->startNodeName = newNode;
      nodes.push_back(new Node(newNode, x, y, edges[pi]->startHeight, 0.0, 0.0, density));
      edges.push_back(new FlowMeter(flowMeterName[i], sn, newNode, density, referenceCrossSection[i], volumeFlowRate));
      int n=edges.size();
      edges[n-1]->startHeight = edges[pi]->startHeight;
      edges[n-1]->endHeight = edges[pi]->startHeight;
      numberEdges++;
      numberNodes++;
    }
    else
    {
      string en = edges[pi]->endNodeName;
      edges[pi]->endNodeName = newNode;
      nodes.push_back(new Node(newNode, x, y, edges[pi]->endHeight, 0.0, 0.0, density));
      edges.push_back(new FlowMeter(flowMeterName[i], en, newNode, density, referenceCrossSection[i], volumeFlowRate));
      int n=edges.size();
      edges[n-1]->startHeight = edges[pi]->endHeight;
      edges[n-1]->endHeight = edges[pi]->endHeight;
      numberEdges++;
      numberNodes++;
    }
  }

  cout << edges[numberEdges-1]->name << "  " << edges[numberEdges-1]->startNodeName << "  " << edges[numberEdges-1]->endNodeName << "  " << edges[numberEdges-1]->typeCode << endl;
  cin.get();

  numberEquations = numberEdges + numberNodes;

  buildSystem();
  buildIndexing();

  // giving initial values to head and volume flow rates
  initialization();

  // resizing Eigen vectors
  x.resize(numberEquations);
  f.resize(numberEquations);

  // Setting initial conditions to x vector
  for(int i=0; i<numberEdges; i++)
    x(i) = edges[i]->volumeFlowRate;
  for(int i=0; i<numberNodes; i++)
    x(numberEdges + i) = nodes[i]->head;

  buildJacobian(); // building the Jacobian matrix
}