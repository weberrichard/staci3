#include "IOxml.h"

//--------------------------------------------------------------------------------
IOxml::IOxml(const char *a_xmlFileName) {
  xmlFileName = a_xmlFileName;

  xmlMain = XMLNode::openFileHelper(xmlFileName, "staci");
  nodeSetting = xmlMain.getChildNode("settings");
  nodeNodes = xmlMain.getChildNode("nodes");
  nodesEdges = xmlMain.getChildNode("edges");

  string cpp_xml_debug = xmlMain.getChildNode("settings").getChildNode("cpp_xml_debug").getText();
  if(cpp_xml_debug == "true")
    debug = true;
  else
    debug = false;
}

//--------------------------------------------------------------------------------
void IOxml::loadSystem(vector<Node *> &nodes, vector<Edge *> &edges) {

  XMLNode nodeSetting = xmlMain.getChildNode("settings");
  XMLNode nodeNodes = xmlMain.getChildNode("nodes");
  XMLNode nodesEdges = xmlMain.getChildNode("edges");

  // Reading the number of nodes and edges
  nodeNumber = nodeNodes.nChildNode("node");
  edgeNumber = nodesEdges.nChildNode("edge");
  if(debug)
    cout << endl << "FileName: " << xmlFileName << endl << "\n Number of nodes: " << nodeNumber << endl << "\tNumber of edges: " << edgeNumber;
  // Reading the nodes in detail
  if(debug)
    cout << endl << endl << "Reading the NODES in detail: " << endl << "-----------------------------------------" << endl;

  string id, is_endnode;
  double height, demand, pressure, density, xcoord, ycoord;

  for(int i = 0; i < nodeNumber; i++){
    is_endnode = nodeNodes.getChildNode("node", i).getChildNode("endnode").getText();
    if(is_endnode == "false"){
      id = nodeNodes.getChildNode("node", i).getChildNode("id").getText();
      xcoord = stod(nodeNodes.getChildNode("node", i).getChildNode("xcoord").getText());
      ycoord = stod(nodeNodes.getChildNode("node", i).getChildNode("ycoord").getText());
      height = stod(nodeNodes.getChildNode("node", i).getChildNode("height").getText());
      demand = stod(nodeNodes.getChildNode("node", i).getChildNode("demand").getText());
      pressure = stod(nodeNodes.getChildNode("node", i).getChildNode("pressure").getText());
      density = stod(nodeNodes.getChildNode("node", i).getChildNode("density").getText());
      //demand = demand / 3.6; // SPR stores in m^3/h, handeld in l/s in STACI
      demand = demand / 3600.; // SPR stores in m^3/h, handeld in m3/s in STACI
      nodes.push_back(new Node(id, xcoord, ycoord, height, demand, pressure, density));

      if(debug)
        cout << nodes.size() << ". node, " << id << ", " << xcoord << ", " << ycoord << ", " << height << ", " << demand << endl;
    }
  }

  // Reading EDGES in detail
  if(debug)
    cout << endl << endl << "Reading EDGES in detail:" << endl << "-----------------------------------------" << endl;

  string node_from, node_to;
  double aref, mass_flow_rate;

  for(int i = 0; i < edgeNumber; i++){
    id = nodesEdges.getChildNode("edge", i).getChildNode("id").getText();
    aref = stod(nodesEdges.getChildNode("edge", i).getChildNode("aref").getText());
    node_from = nodesEdges.getChildNode("edge", i).getChildNode("node_from").getText();
    node_to = nodesEdges.getChildNode("edge", i).getChildNode("node_to").getText();
    density = stod(nodesEdges.getChildNode("edge", i).getChildNode("density").getText());
    mass_flow_rate = stod(nodesEdges.getChildNode("edge", i).getChildNode("mass_flow_rate").getText());

    if(debug)
      cout << endl << i << ". edge: " << id << ", node_from = " << node_from << ", node_to = " << node_to << ", density = " << density << ", aref = " << aref;

    string edgeType = nodesEdges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(0).getName();
    XMLNode nodeEdgeSpec = nodesEdges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(edgeType.c_str());
    if(edgeType == "press") // PRESSUREPOINT
    {
      double pres = stod(nodeEdgeSpec.getChildNode("pressure").getText());
      edges.push_back(new PressurePoint(id, aref, node_from, density, pres/9.81/1000., 0.0, mass_flow_rate));
    }
    else if(edgeType == "pipe") // PIPELINE
    {
      double L = stod(nodeEdgeSpec.getChildNode("length").getText());
      double D = stod(nodeEdgeSpec.getChildNode("diameter").getText());
      double roughness = stod(nodeEdgeSpec.getChildNode("roughness").getText());
      edges.push_back(new Pipe(id, node_from, node_to, density, L, D, roughness, mass_flow_rate, false, 2, 1.0));
    }
    else if(edgeType == "pool") // POOL
    {
      double Hb = stod(nodeEdgeSpec.getChildNode("bottom_level").getText());
      double Hw = stod(nodeEdgeSpec.getChildNode("water_level").getText());
      edges.push_back(new Pool(id, node_from, density, aref, Hb, Hw, mass_flow_rate));
    }
    else if(edgeType == "valve") // VALVE
    {
      edges.push_back(new ValveISO(id, node_from, node_to, density, aref, mass_flow_rate));
    }
    else
    {
      cout << endl << " !!!WARNING!!! " << endl << " Reading from XML file, unknown edge type: " << edgeType << "  ID: " << id << endl << " Continouing..." << endl;
    }
  }
}

//--------------------------------------------------------------------------------
void IOxml::loadInitialValue(vector<Node *> &nodes, vector<Edge *> &edges) {

  nodeNumber = nodeNodes.nChildNode("node");
  edgeNumber = nodesEdges.nChildNode("edge");

  string id;
  double head;
  for(int i = 0; i < nodeNumber; i++){
    id = nodeNodes.getChildNode("node", i).getChildNode("id").getText();
    for(int j = 0; j < nodes.size(); j++){
      if(nodes.at(j)->name == id){
        head = stod(nodeNodes.getChildNode("node", i).getChildNode("pressure").getText());
        nodes[j]->setProperty("head",head);
      }
    }
  }

  double mass_flow_rate;
  for (int i = 0; i < edgeNumber; i++) {
    id = nodesEdges.getChildNode("edge", i).getChildNode("id").getText();
    for (unsigned int j = 0; j < edges.size(); j++) {
      if (edges.at(j)->getEdgeStringProperty("name") == id) {
        mass_flow_rate = stod(nodesEdges.getChildNode("edge", i).getChildNode("mass_flow_rate").getText());
        edges[j]->setDoubleProperty("massFlowRate",mass_flow_rate);
      }
    }
  }
}

//--------------------------------------------------------------------------------
string IOxml::readSetting(string which) {

  XMLNode nodeSetting = xmlMain.getChildNode("settings");

  string out;
  int number = nodeSetting.nChildNode(which.c_str());

  if(number == 1)
    out = xmlMain.getChildNode("settings").getChildNode(which.c_str()).getText();
  else if(number >= 2){
    cout << endl << "!!WARNING!! IOXML: setting (" + which + ") found multiple times!" << endl;
    out = "0";
  }
  else{
    cout << endl << "!!WARNING!! IOXML: setting (" + which + ") not found!" << endl;
    out = "0";
  }

  return out;
}

/*! WR: function to change the value in a tag inside on the second level + handles the existince of the tag (if there is a new tag it will be added)
*/
void IOxml::writingTagValue(XMLNode node, string tag_name, double tag_value){
    ostringstream os;
    if(node.getChildNode(tag_name.c_str()).isEmpty())
    {
        node.addChild(XMLNode::parseString(("<" + tag_name + "> </" + tag_name + ">").c_str()));
        if(debug)
            cout << "\n The following child has been added: " << tag_name << '\n';
    }
    node.getChildNode(tag_name.c_str()).deleteText();
    os.str("");
    os << scientific << setprecision(5) << tag_value;
    node.getChildNode(tag_name.c_str()).addText(os.str().c_str());
}

//--------------------------------------------------------------------------------
void IOxml::curveReader(const string id, const XMLNode nodeEdgeSpec, vector<double> &px, vector<double> &py) {

  string curve_id = nodeEdgeSpec.getChildNode("id").getText();
  string x_val = nodeEdgeSpec.getChildNode("x_val").getText();
  string y_val = nodeEdgeSpec.getChildNode("y_val").getText();
  string x_dim = nodeEdgeSpec.getChildNode("x_dim").getText();
  string y_dim = nodeEdgeSpec.getChildNode("y_dim").getText();
  int numberX = nodeEdgeSpec.getChildNode("points").nChildNode("point_x");
  int numberY = nodeEdgeSpec.getChildNode("points").nChildNode("point_y");

  if(numberX == numberY)
  {
    double x, y;
    for (int i = 0; i < numberX; i++) {
      x = stod(nodeEdgeSpec.getChildNode("points").getChildNode("point_x", i).getText());
      y = stod(nodeEdgeSpec.getChildNode("points").getChildNode("point_y", i).getText());

      px.push_back(x);
      py.push_back(y);
    }
  }
  else
    cout << endl << "xml reading: ERROR! number of X =" << numberX << " <-> number of Y =" << numberY << endl << "Moving on " << id << " element..." << endl << endl;
}