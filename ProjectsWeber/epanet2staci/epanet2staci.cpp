/*===================================================================*\
                              epanet2staci  
                            ----------------

  This code converts an EPANET file (.net) to a STACI file (.spr).
  Importing nodes, pipelines, valves, pumps, pools, pressure points.
  Converts every units to SI.
  
  HOW TO USE: 

      - CASE_NAME is the name of the file without .net extension

      - SET:  nothing :)

      - OUTPUT: CASE_NAME.spr            File for STACI

      - RUN: ./epanet2staci CASE_NAME

      - COMPILE: clang++ epanet2staci.cpp -o epanet2staci.out

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

vector<string> line2sv(string line);
double VectorMax(const vector<double> &v);
double VectorMin(const vector<double> &v);

int main(int argc, char* argv[]){

// ###################################
// SETTINGS SETTINGS SETTINGS SETTINGS
// ###################################
	bool do_print = false; // printing to terminal
	double dem_pat_scal = 1.0; // scaling the demands 

	string case_name;

	time_t ido=time(0);
	if(argc == 1){
		case_name = "Net1";
  }else if(argc == 2){
    case_name = argv[1];
  }

	int warning_counter=0;

	cout << "\nEPAENT ==>> STACI Converter by WR" << endl;
	cout << "\nFile in: " << case_name << endl;

	// FOR NODES
	vector<string> node, node_dem_pat;
	vector<double> elev, dem, x, y;

	// FOR PIPES
	vector<string> pipe, node_from, node_to;
	vector<double> l, D, roughness;

	// FOR PRES
	vector<string> pres_id, pres_node_from, pres_node_to;
	vector<double> pres_p;

	// FOR POOLS
	vector<string> pool, pool_node_from, pool_node_to;
	vector<double> pool_botlev, pool_watlev, pool_aref;

	// FOR PUMPS
	vector<string> pump, pump_node_from, pump_node_to, pump_par;
	vector<vector<double> > pump_cv_x, pump_cv_y;

	// FOR VALVES
	vector<string> valve, valve_node_from, valve_node_to, valve_type;
	vector<double> valve_d, valve_set;

	// FOR DEMAND PATTERNS
	vector<string> dem_pat_id;
	vector<vector<double> > dem_pat;

	// FOR CURVES
	vector<string> curve;
	vector<vector<double> > curve_x, curve_y;

	// SETTINGS
	string flow_unit;
	string friction_model;

// ################
// READING THE FILE
// ################
	ifstream file_in;
  file_in.open(case_name+".inp");
  string line;
  if(file_in.is_open()){
	  while(getline(file_in,line)){
			if(line.substr(0,11) == "[JUNCTIONS]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						node.push_back("NODE_"+sv[0]);
						elev.push_back(stod(sv[1],0));
						dem.push_back(stod(sv[2],0));
						x.push_back(0.);y.push_back(0.);
						if(sv.size()==4)
							node_dem_pat.push_back("DEM_PAT"+sv[3]);
						else
							node_dem_pat.push_back("DEFAULT");
					}
				}
			}
			if(line.substr(0,12) == "[RESERVOIRS]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						pres_id.push_back("PRES_"+sv[0]);
						pres_p.push_back(stod(sv[1],0));
						// Creating the node
						node.push_back("NODE_"+sv[0]);
						elev.push_back(-1.); // no data in inp, will be equal to other end's
						dem.push_back(0.);
						x.push_back(0.);y.push_back(0.);
						node_dem_pat.push_back("DEFAULT");
						pres_node_to.push_back(node.back()+"_END");
						pres_node_from.push_back(node.back());
						node.push_back("NODE_"+sv[0]+"_END");
						elev.push_back(-1.);
						dem.push_back(0.);
						x.push_back(0.);y.push_back(0.);
						node_dem_pat.push_back("DEFAULT");
					}
				}
			}

			if(line.substr(0,7) == "[TANKS]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						pool.push_back("POOL_"+sv[0]);
						pool_botlev.push_back(stod(sv[1],0));
						pool_watlev.push_back(stod(sv[2],0));
						double d = stod(sv[5],0);
						pool_aref.push_back(d*d*M_PI/4);
						// TODO: minlevel, maxlevel etc.
						// Creating the node
						node.push_back("NODE_"+sv[0]);
						elev.push_back(stod(sv[1],0));
						dem.push_back(0.);
						x.push_back(0.);y.push_back(0.);
						node_dem_pat.push_back("DEFAULT");
						pool_node_to.push_back(node.back()+"_END");
						pool_node_from.push_back(node.back());
						node.push_back("NODE_"+sv[0]+"_END");
						elev.push_back(stod(sv[1],0));
						dem.push_back(0.);
						x.push_back(0.);y.push_back(0.);
						node_dem_pat.push_back("DEFAULT");
					}
				}
			}

			if(line.substr(0,7) == "[PIPES]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						pipe.push_back("PIPE_"+sv[0]);
						node_from.push_back("NODE_"+sv[1]);
						node_to.push_back("NODE_"+sv[2]);
						l.push_back(stod(sv[3],0));
						D.push_back(stod(sv[4],0));
						roughness.push_back(stod(sv[5],0));
					}
				}
			}

			if(line.substr(0,7) == "[PUMPS]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						pump.push_back("PUMP_"+sv[0]);
						pump_node_from.push_back("NODE_"+sv[1]);
						pump_node_to.push_back("NODE_"+sv[2]);
						vector<double> v;
						pump_cv_x.push_back(v);
						pump_cv_y.push_back(v);
						if(sv.size()<=4){
							cout << "\n !WARNING! No pump curve, or power is defined for " << "PUMP_"+sv[0] << endl;
							warning_counter++;
						}
						if(sv.size()==5){
							if(sv[3]=="HEAD")
								pump_par.push_back("CV_"+sv[4]);
						} // TODO later: add POWER, SPEED etc.
					}
				}
			}

			if(line.substr(0,8) == "[VALVES]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						valve.push_back("VALVE_"+sv[0]);
						valve_node_from.push_back("NODE_"+sv[1]);
						valve_node_to.push_back("NODE_"+sv[2]);
						valve_d.push_back(stod(sv[3],0));
						valve_type.push_back(sv[4]);
						valve_set.push_back(stod(sv[5],0));
					}
				}
			}

			if(line.substr(0,10) == "[PATTERNS]"){
				vector<double> v;
				while(true){
					getline(file_in,line);
					if(line[0] != ';'){
						if(line.length()<=1 || line.substr(0,5) == "[END]"){
							dem_pat.push_back(v);
							break;
						}
						vector<string> sv=line2sv(line);
						if(dem_pat_id.size()==0){
							dem_pat_id.push_back("DP_"+sv[0]);
							for(vector<string>::iterator i=sv.begin()+1; i!=sv.end(); i++)
								v.push_back(stod(*i,0));
						}else{
							if("DP_"+sv[0] == dem_pat_id.back()){
								for(vector<string>::iterator i=sv.begin()+1; i!=sv.end(); i++)
									v.push_back(stod(*i,0));
							}else{
								dem_pat.push_back(v);
								v.clear();
								dem_pat_id.push_back("DP_"+sv[0]);
								for(vector<string>::iterator i=sv.begin()+1; i!=sv.end(); i++)
									v.push_back(stod(*i,0));
							}
						}
					}
				}
			}

			// PUMP CURVES
			if(line.substr(0,8) == "[CURVES]"){
				vector<double> x,y;
				while(true){
					getline(file_in,line);
					if(line[0] != ';'){
						if(line.length()<=1 || line.substr(0,5) == "[END]"){
							curve_x.push_back(x);
							curve_y.push_back(y);
							break;
						}
						vector<string> sv=line2sv(line);
						if(curve.size()==0){
							curve.push_back("CV_"+sv[0]);
							x.push_back(stod(sv[1],0));
							y.push_back(stod(sv[2],0));
						}else{
							if("CV_"+sv[0] == curve.back()){
								x.push_back(stod(sv[1],0));
								y.push_back(stod(sv[2],0));
							}else{
								curve_x.push_back(x);
								curve_y.push_back(y);
								x.clear(); y.clear();
								curve.push_back("CV_"+sv[0]);
								x.push_back(stod(sv[1],0));
								y.push_back(stod(sv[2],0));
							}
						}
					}
				}
			}

			// BASIC SETTINGS
			if(line.substr(0,9) == "[OPTIONS]"){
			//TODO There are more options that might be intresting later
				vector<double> x,y;
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						if(sv.size()>=2){
							if(sv[0] == "Units"){
								flow_unit = sv[1];
							}
							if(sv[0] == "Headloss"){ //TODO D-W, C-W
								if(sv[1] == "H-W")
									friction_model =  "HW";
							}
						}
					}
				}
			}

			if(line.substr(0,13) == "[COORDINATES]"){
				while(true){
					getline(file_in,line);
					if(line.length()<=1 || line.substr(0,5) == "[END]")
						break;
					if(line[0] != ';'){
						vector<string> sv=line2sv(line);
						int idx=-1;
						for(int i=0; i<node.size(); i++)
							if(node[i]=="NODE_"+sv[0]){
								idx = i;
								break;
							}
						if(idx==-1){
							cout << "\n !WARNING! node not found for coordinate " << "NODE_"+sv[0] << endl;
							warning_counter++;
						}
						else{
							x[idx] = stod(sv[1],0);
							y[idx] = stod(sv[2],0);
						}
					}
				}
			}
	  }
  }else{
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl << "epanet2staci(), File is not open when calling epanet2staci() function!!! file: " << case_name << endl;
		exit(-1);
  }
  file_in.close();
  cout << "\n[*]Reading EPANET file: OK" << endl;

// #######################
// POSTPROCESSING THE DATA
// #######################	
  //Setting the elevation for PRES nodes
  for(int i=0; i<pres_id.size(); i++){
  	for(int j=0; j<pipe.size(); j++){
 			if(node_from[j].substr(5,node_from[j].length()-5)==pres_id[i].substr(5,pres_id[i].length()-5)){
 				int idx_to=-1, idx_from=-1;
 				for(int k=0; k<node.size(); k++){
 					if(node_to[j]==node[k])
 						idx_to = k;
 					if(node_from[j]==node[k])
 						idx_from = k;
 					if(idx_from != -1 && idx_to != -1)
 						break;
 				}
 				elev[idx_from] = elev[idx_to];
 				elev[idx_from+1] = elev[idx_to];
 				break;
 			}
 			if(node_to[j].substr(5,node_to[j].length()-5)==pres_id[i].substr(5,pres_id[i].length()-5)){
 				int idx_to=-1, idx_from=-1;
 				for(int k=0; k<node.size(); k++){
 					if(node_to[j]==node[k])
 						idx_to = k;
 					if(node_from[j]==node[k])
 						idx_from = k;
 					if(idx_from != -1 && idx_to != -1)
 						break;
 				}
 				elev[idx_to] = elev[idx_from];
 				elev[idx_to+1] = elev[idx_from];
 				break;
 			}
  	}
  	for(int j=0; j<pump.size(); j++){
 			if(pump_node_from[j].substr(5,pump_node_from[j].length()-5)==pres_id[i].substr(5,pres_id[i].length()-5)){
 				int idx_to=-1, idx_from=-1;
 				for(int k=0; k<node.size(); k++){
 					if(pump_node_to[j]==node[k])
 						idx_to = k;
 					if(pump_node_from[j]==node[k])
 						idx_from = k;
 					if(idx_from != -1 && idx_to != -1)
 						break;
 				}
 				elev[idx_from] = elev[idx_to];
 				elev[idx_from+1] = elev[idx_to];
 				break;
 			}
 			if(pump_node_to[j].substr(5,pump_node_to[j].length()-5)==pres_id[i].substr(5,pres_id[i].length()-5)){
 				int idx_to=-1, idx_from=-1;
 				for(int k=0; k<node.size(); k++){
 					if(pump_node_to[j]==node[k])
 						idx_to = k;
 					if(pump_node_from[j]==node[k])
 						idx_from = k;
 					if(idx_from != -1 && idx_to != -1)
 						break;
 				}
 				elev[idx_to] = elev[idx_from];
 				elev[idx_to+1] = elev[idx_from];
 				break;
 			}
 		}
  }

  // Setting x-y coordinates for dummy nodes of pools and press
  for(int i=0; i<node.size(); i++){
  	if(x[i]==0){
  		x[i] = x[i-1];
  		y[i] = y[i-1];
  	}
  }

  // Matching curves to pumps
  for(int i=0; i<pump.size(); i++){
  	for(int j=0; j<curve.size(); j++){
  		if(pump_par[i]==curve[j]){
  			for(int k=0; k<curve_x[j].size(); k++){
  				pump_cv_x[i].push_back(curve_x[j][k]);
  				pump_cv_y[i].push_back(curve_y[j][k]);
  			}
  		}
  	}
  }

  // Modifying pump curves, STACI must have at least 3 points
  for(int i=0; i<pump.size(); i++)
  	if(pump_cv_x[i].size()<=2){
  		cout << "\n !WARNING! Pump " << pump[i] << "does not have enough points. STACI needs at least 3 points\n";
  		warning_counter++;
  	}

  cout << "\nflow_unit: " << flow_unit << endl;
  // Converting units from EPANET to STACI (m3/h)
  double dem_unit;
  string unit;
  // US units
  if(flow_unit == "CFS"){ // cubic feet seconds
  	dem_unit = 101.940648;
  	unit = "US";
  }
  if(flow_unit == "GPM"){ // gallon(US) per minute
  	dem_unit = 0.227125;
  	unit = "US";
  }
  if(flow_unit == "MGD"){ // million gallon(US) per day
  	dem_unit = 157.72549166667;
  	unit = "US";
  }
  if(flow_unit == "IMGD"){ // imperial million gallon per day
  	dem_unit = 189.42041666667;
  	unit = "US";
  }
  if(flow_unit == "AFD"){ // acre-feet / day
  	dem_unit = 51.4;
  	unit = "US";
  }
  // NORMAL units
  if(flow_unit == "LPS"){ // litre per sec
  	dem_unit = 3.6;
  	unit = "SI";
  }
  if(flow_unit == "LPM"){ // litre per minute
  	dem_unit = 0.06;
  	unit = "SI";
  }
  if(flow_unit == "MLD"){ // million litre per day
  	dem_unit = 1000./24.;
  	unit = "SI";
  }
  if(flow_unit == "CMH"){ // cubic meter per hour
  	dem_unit = 1.;
  	unit = "SI";
  }
  if(flow_unit == "CMD"){ // cubic meter per day
  	dem_unit = 1./24.;
  	unit = "SI";
  }

  for(int i=0; i<node.size(); i++)
  	dem[i] = dem[i]*dem_unit*dem_pat_scal;
  if(unit=="US"){
		for(int i=0; i<node.size(); i++)
	  		elev[i] = elev[i]*0.3048; // feet to meter
	  for(int i=0; i<pipe.size(); i++){
  		l[i] = l[i]*0.3048; // feet to meter
  		D[i] = D[i]*0.0254; // inches to meter
  		 // TODO Roughness in D-W
  	}
		for(int i=0; i<pres_id.size(); i++)
		  pres_p[i] = pres_p[i]*0.3048; // feet to meter
		for(int i=0; i<pool.size(); i++){
			pool_botlev[i] = pool_botlev[i]*0.3048; // feet to meter
			pool_watlev[i] = pool_watlev[i]*0.3048; // feet to meter
			pool_aref[i] = pool_aref[i]*0.3048*0.3048; // square feet to square meter
	  }
	  for(int i=0; i<valve.size(); i++)
  		valve_d[i] = valve_d[i]*0.0254; // inches to meter
	  for(int i=0; i<pump.size(); i++){
  		for(int j=0; j<pump_cv_x[i].size(); j++){
  			pump_cv_x[i][j] = pump_cv_x[i][j]*dem_unit; // feet to meter
  			pump_cv_y[i][j] = pump_cv_y[i][j]*0.3048; // feet to meter
  		}
	  }
  }

	if(unit=="SI"){
		for(int i=0; i<pipe.size(); i++)
	  	D[i] = D[i]/1000.; // millimeter to meter
	  for(int i=0; i<valve.size(); i++)
  		valve_d[i] = valve_d[i]/1000.; // millimeter to meter
	}

	for(int i=0; i<pres_id.size(); i++){
		int idx=-1;
		for(int j=0; j<node.size(); j++)
			if(node[j]==pres_node_from[i])
				idx = j;

		pres_p[i] = (pres_p[i]-elev[idx])*9.81*1000; // meter to Pascal minus elevation
	}

  cout << "\n[*]Processing data: OK" << endl;

// ############################################
// PRINTING THE RESULTS FROM MEMORY TO TERMINAL
// ############################################
  if(do_print){
	  cout << "\nCase name: "+case_name+"\n";
	  printf("====================================\n");
	  if(friction_model=="HW")
		  printf("Friction model: Hazen-Williams\n");

	  printf("\n\n[*] NODES");
	  printf("\n %10s | %10s | %10s | %10s | %11s | %9s |", "ID", "Elevation", "Demand", "Demand_pat", "X", "Y");
	  for(int i=0; i<node.size(); i++)
		  printf("\n %10s | %10.2f | %10.2f | %10s | %10.2f|%10.2f |", node[i].c_str(), elev[i], dem[i], node_dem_pat[i].c_str(),x[i],y[i]);

		printf("\n\n[*] PIPES");
	  printf("\n %10s | %10s | %10s | %10s | %10s | %10s |", "ID", "Node_from", "Node_to", "Length", "Diameter", "Roughness");
	  for(int i=0; i<pipe.size(); i++)
		  printf("\n %10s | %10s | %10s | %10.2f | %10.2f | %10.2f |", pipe[i].c_str(), node_from[i].c_str(), node_to[i].c_str(), l[i], D[i], roughness[i]);

		printf("\n\n[*] PRES");
	  printf("\n %10s | %10s | %10s | %10s |", "ID", "Node_from", "Node_to", "Head");
	  for(int i=0; i<pres_id.size(); i++)
		  printf("\n %10s | %10s | %10s | %10.2f |", pres_id[i].c_str(), pres_node_from[i].c_str(), pres_node_to[i].c_str(), pres_p[i]);

		printf("\n\n[*] POOLS");
	  printf("\n %10s | %10s | %10s | %10s | %10s | %10s |", "ID", "Node_from", "Node_to", "BotLev", "WatLev", "Aref");
	  for(int i=0; i<pres_id.size(); i++)
		  printf("\n %10s | %10s | %10s | %10.2f | %10.2f | %10.2f |", pool[i].c_str(), pool_node_from[i].c_str(), pool_node_to[i].c_str(), pool_botlev[i], pool_watlev[i], pool_aref[i]);

		printf("\n\n[*] VALVES");
	  printf("\n %10s | %10s | %10s | %10s | %10s | %10s |", "ID", "Node_from", "Node_to", "Diameter", "Type", "Setting");
	  for(int i=0; i<valve.size(); i++)
		  printf("\n %10s | %10s | %10s | %10.2f | %10s | %10.2f |", valve[i].c_str(), valve_node_from[i].c_str(), valve_node_to[i].c_str(), valve_d[i],valve_type[i].c_str(), valve_set[i]);

		printf("\n\n[*] DEMAND PATTERNS");
	  printf("\n %10s | %10s |", "ID", "Multipl");
	  for(int i=0; i<dem_pat_id.size(); i++){
		  printf("\n %10s |", dem_pat_id[i].c_str());
		  for(int j=0; j<dem_pat[i].size(); j++)
		  	printf(" %10.4f |", dem_pat[i][j]);
	  }

	  printf("\n\n[*] PUMPS");
	  printf("\n %10s | %10s | %10s | %10s |", "ID", "Node_from", "Node_to", "Curve");
	  for(int i=0; i<pump.size(); i++){
		  printf("\n %10s | %10s | %10s |", pump[i].c_str(), pump_node_from[i].c_str(), pump_node_to[i].c_str());
		  for(int j=0; j<pump_cv_x[i].size(); j++)
		  	printf(" %10.2f-%10.2f |", pump_cv_x[i][j], pump_cv_y[i][j]);
	  }

	  printf("\n\n[*] CURVES");
	  printf("\n %10s | %10s |", "ID", "X - Y");
	  for(int i=0; i<curve.size(); i++){
		  printf("\n %10s |", curve[i].c_str());
		  for(int j=0; j<curve_x[i].size(); j++)
		  	printf(" %10.2f-%10.2f |", curve_x[i][j], curve_y[i][j]);
	  }
	}

// ######################
// WRITING THE STACI FILE
// ######################
  int case_width = 800, case_height = 800;
  ofstream file_out;
  file_out.open(case_name + ".spr");
  file_out << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\r\n<staci>\r\n";
  file_out << "\t<settings>\r\n\t\t<id>\r\n\t\t\t"+case_name+"\r\n\t\t</id>\r\n\t\t<name>\r\n\t\t\t"+case_name+"\r\n\t\t</name>\r\n\t\t<color/>\r\n\t\t<counter>\r\n\t\t\t1\r\n\t\t</counter>\r\n\t\t<user1>\r\n\t\t\t0\r\n\t\t</user1>\r\n\t\t<user2>\r\n\t\t\t0\r\n\t\t</user2>\r\n\t\t<zone>\r\n\t\t\t0\r\n\t\t</zone>\r\n\t\t<cpp_xml_debug>\r\n\t\t\tfalse\r\n\t\t</cpp_xml_debug>\r\n\t\t<solution_exists>\r\n\t\t\ttrue\r\n\t\t</solution_exists>\r\n\t\t<tsolution_exists>\r\n\t\t\tfalse\r\n\t\t</tsolution_exists>\r\n\t\t<csolution_exists>\r\n\t\t\tfalse\r\n\t\t</csolution_exists>\r\n\t\t<serial_exists>\r\n\t\t\tfalse\r\n\t\t</serial_exists>\r\n\t\t<out_file>\r\n\t\t\tstaci.out\r\n\t\t</out_file>\r\n\t\t<debug_level>\r\n\t\t\t2\r\n\t\t</debug_level>\r\n\t\t<friction_model>\r\n\t\t\t"+friction_model+"\r\n\t\t</friction_model>\r\n\t\t<version>\r\n\t\t\t1\r\n\t\t</version>\r\n\t\t<xstart>\r\n\t\t\t0\r\n\t\t</xstart>\r\n\t\t<ystart>\r\n\t\t\t0\r\n\t\t</ystart>\r\n\t\t<xend>\r\n\t\t\t100\r\n\t\t</xend>\r\n\t\t<yend>\r\n\t\t\t100\r\n\t\t</yend>\r\n\t\t<owner>\r\n\t\t\trweber\r\n\t\t</owner>\r\n\t\t<width>\r\n\t\t\t"+to_string(case_width)+"\r\n\t\t</width>\r\n\t\t<height>\r\n\t\t\t"+to_string(case_height)+"\r\n\t\t</height>\r\n\t\t<hasmap>\r\n\t\t\tfalse\r\n\t\t</hasmap>\r\n\t\t<image/>\r\n\t\t<locked>\r\n\t\t\tfalse\r\n\t\t</locked>\r\n\t\t<iter_max>\r\n\t\t\t100\r\n\t\t</iter_max>\r\n\t\t<e_p_max>\r\n\t\t\t0.0001\r\n\t\t</e_p_max>\r\n\t\t<e_mp_max>\r\n\t\t\t0.01\r\n\t\t</e_mp_max>\r\n\t\t<relax>\r\n\t\t\t1\r\n\t\t</relax>\r\n\t\t<relax_mul>\r\n\t\t\t1.2\r\n\t\t</relax_mul>\r\n\t\t<tt_length>\r\n\t\t\t1.0\r\n\t\t</tt_length>\r\n\t\t<cl_length>\r\n\t\t\t1.0\r\n\t\t</cl_length>\r\n\t\t<mp_init>\r\n\t\t\t100.0\r\n\t\t</mp_init>\r\n\t\t<p_init>\r\n\t\t\t5.0\r\n\t\t</p_init>\r\n\t\t<sum_of_inflow>\r\n\t\t\t0.0\r\n\t\t</sum_of_inflow>\r\n\t\t<sum_of_demand>\r\n\t\t\t0.0\r\n\t\t</sum_of_demand>\r\n\t\t<fluid_volume>\r\n\t\t\t1.41372e+02\r\n\t\t</fluid_volume>\r\n\t\t<reservoir_capacity>\r\n\t\t\t0.0\r\n\t\t</reservoir_capacity>\r\n\t\t<start_date>\r\n\t\t\t2007-01-01\r\n\t\t</start_date>\r\n\t\t<start_time>\r\n\t\t\t0.0\r\n\t\t</start_time>\r\n\t\t<time_step>\r\n\t\t\t1.0\r\n\t\t</time_step>\r\n\t\t<step_count>\r\n\t\t\t24\r\n\t\t</step_count>\r\n\t\t<play_step>\r\n\t\t\t10.0\r\n\t\t</play_step>\r\n\t\t<muls>\r\n\t\t\t<mul>\r\n\t\t\t\t0.0\r\n\t\t\t</mul>\r\n\t\t</muls>\r\n\t</settings>\r\n";
  file_out << "\t<texts/>\r\n";

  // NODES
  file_out << "\t<nodes>\r\n";
  double xmin = VectorMin(x), xmax = VectorMax(x), ymin = VectorMin(y), ymax = VectorMax(y);
  for(int i=0;i<node.size();i++){
    double xpos,ypos;
    double delta = max(ymax-ymin,xmax-xmin);
    xpos = (x[i]-xmin)/delta*(case_width-20)+10;
    ypos = -((y[i]-ymin)/delta*(case_height-20)+10)+case_height;
    int color = -8684677;
    string endnode="false";
    if(node[i].substr(node[i].length()-3) == "END")
    	endnode="true";
    file_out << "\t\t<node>\r\n\t\t\t<id>\r\n\t\t\t\t"+node[i]+"\r\n\t\t\t</id>\r\n\t\t\t<name>\r\n\t\t\t\t"+node[i]+"\r\n\t\t\t</name>\r\n\t\t\t<color>\r\n\t\t\t\t"+to_string(color)+"\r\n\t\t\t</color>\r\n\t\t\t<counter>\r\n\t\t\t\t14\r\n\t\t\t</counter>\r\n\t\t\t<user1>\r\n\t\t\t\t0\r\n\t\t\t</user1>\r\n\t\t\t<user2>\r\n\t\t\t\t0\r\n\t\t\t</user2>\r\n\t\t\t<zone>\r\n\t\t\t\t0\r\n\t\t\t</zone>\r\n\t\t\t<xcoord>\r\n\t\t\t\t"+to_string(x[i])+"\r\n\t\t\t</xcoord>\r\n\t\t\t<ycoord>\r\n\t\t\t\t"+to_string(y[i])+"\r\n\t\t\t</ycoord>\r\n\t\t\t<xpos>\r\n\t\t\t\t"+to_string(xpos)+"\r\n\t\t\t</xpos>\r\n\t\t\t<ypos>\r\n\t\t\t\t"+to_string(ypos)+"\r\n\t\t\t</ypos>\r\n\t\t\t<height>\r\n\t\t\t\t"+to_string(elev[i])+"\r\n\t\t\t</height>\r\n\t\t\t<demand>\r\n\t\t\t\t"+to_string(dem[i])+"\r\n\t\t\t</demand>\r\n\t\t\t<density>\r\n\t\t\t\t1000.0\r\n\t\t\t</density>\r\n\t\t\t<cl_input>\r\n\t\t\t\t0.0\r\n\t\t\t</cl_input>\r\n\t\t\t<pressure>\r\n\t\t\t\t0.00000\r\n\t\t\t</pressure>\r\n\t\t\t<head>\r\n\t\t\t\t0.00000\r\n\t\t\t</head>\r\n\t\t\t<endnode>\r\n\t\t\t\t"+endnode+"\r\n\t\t\t</endnode>\r\n\t\t\t<travel_time>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</travel_time>\r\n\t\t\t<concentration>\r\n\t\t\t\t0.0000\r\n\t\t\t</concentration>\r\n\t\t\t<curve>\r\n\t\t\t\t<id>\r\n\t\t\t\t\tmuls\r\n\t\t\t\t</id>\r\n\t\t\t\t<x_val/>\r\n\t\t\t\t<x_dim/>\r\n\t\t\t\t<y_val/>\r\n\t\t\t\t<y_dim/>\r\n\t\t\t\t<points>\r\n\t\t\t\t\t<point_x>\r\n\t\t\t\t\t\t0,0000\r\n\t\t\t\t\t</point_x>\r\n\t\t\t\t\t<point_y>\r\n\t\t\t\t\t\t0,0000\r\n\t\t\t\t\t</point_y>\r\n\t\t\t\t</points>\r\n\t\t\t</curve>\r\n\t\t</node>\r\n";
  }
  file_out << "\t</nodes>\r\n";

  // EDGES
  file_out << "\t<edges>\r\n";
  // PIPES
  for(int i=0; i<pipe.size(); i++){
    string edge_name = "PIPE_" + node_from[i] + "_to_" + node_to[i];
    double aref = D[i]*D[i]/4*M_PI;
    file_out << "\t\t<edge>\r\n\t\t\t<id>\r\n\t\t\t\t"+pipe[i]+"\r\n\t\t\t</id>\r\n\t\t\t<name>\r\n\t\t\t\t"+edge_name+"\r\n\t\t\t</name>\r\n\t\t\t<color>\r\n\t\t\t\t0\r\n\t\t\t</color>\r\n\t\t\t<counter>\r\n\t\t\t\t19\r\n\t\t\t</counter>\r\n\t\t\t<user1>\r\n\t\t\t\t0\r\n\t\t\t</user1>\r\n\t\t\t<user2>\r\n\t\t\t\t0\r\n\t\t\t</user2>\r\n\t\t\t<zone>\r\n\t\t\t\t0\r\n\t\t\t</zone>\r\n\t\t\t<pipe_type>\r\n\t\t\t\tpipe\r\n\t\t\t</pipe_type>\r\n\t\t\t<node_from>\r\n\t\t\t\t"+node_from[i]+"\r\n\t\t\t</node_from>\r\n\t\t\t<node_to>\r\n\t\t\t\t"+node_to[i]+"\r\n\t\t\t</node_to>\r\n\t\t\t<density>\r\n\t\t\t\t1000\r\n\t\t\t</density>\r\n\t\t\t<aref>\r\n\t\t\t\t"+to_string(aref)+"\r\n\t\t\t</aref>\r\n\t\t\t<mass_flow_rate>\r\n\t\t\t\t0.00\r\n\t\t\t</mass_flow_rate>\r\n\t\t\t<volume_flow_rate>\r\n\t\t\t\t0.00\r\n\t\t\t</volume_flow_rate>\r\n\t\t\t<velocity>\r\n\t\t\t\t0.00\r\n\t\t\t</velocity>\r\n\t\t\t<headloss>\r\n\t\t\t\t0.00\r\n\t\t\t</headloss>\r\n\t\t\t<head_loss_per_unit_length>\r\n\t\t\t\t0.00\r\n\t\t\t</head_loss_per_unit_length>\r\n\t\t\t<travel_time>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</travel_time>\r\n\t\t\t<concentration>\r\n\t\t\t\t0.0000\r\n\t\t\t</concentration>\r\n\t\t\t<edge_spec>\r\n\t\t\t\t<pipe>\r\n\t\t\t\t\t<length>\r\n\t\t\t\t\t\t"+to_string(l[i])+"\r\n\t\t\t\t\t</length>\r\n\t\t\t\t\t<diameter>\r\n\t\t\t\t\t\t"+to_string(D[i])+"\r\n\t\t\t\t\t</diameter>\r\n\t\t\t\t\t<roughness>\r\n\t\t\t\t\t\t"+to_string(roughness[i])+"\r\n\t\t\t\t\t</roughness>\r\n\t\t\t\t\t<friction_coeff>\r\n\t\t\t\t\t\t"+to_string(-0.02)+"\r\n\t\t\t\t\t</friction_coeff>\r\n\t\t\t\t\t<cl_k>\r\n\t\t\t\t\t\t0.016\r\n\t\t\t\t\t</cl_k>\r\n\t\t\t\t\t<cl_w>\r\n\t\t\t\t\t\t0.0\r\n\t\t\t\t\t</cl_w>\r\n\t\t\t\t</pipe>\r\n\t\t\t</edge_spec>\r\n\t\t</edge>\r\n";
  }
  // VALVES
  for(int i=0; i<valve.size(); i++){
	  int color = -16711936;
		string valve_name = "VALVE_" + valve_node_from[i] + "_to_" + valve_node_to[i];
		double aref = valve_d[i]*valve_d[i]/4.*M_PI;
	  file_out << "\t\t<edge>\r\n\t\t\t<id>\r\n\t\t\t\t"+valve[i]+"\r\n\t\t\t</id>\r\n\t\t\t<name>\r\n\t\t\t\t"+valve_name+"\r\n\t\t\t</name>\r\n\t\t\t<color>\r\n\t\t\t\t"+to_string(color)+"\r\n\t\t\t</color>\r\n\t\t\t<counter>\r\n\t\t\t\t0\r\n\t\t\t</counter>\r\n\t\t\t<pipe_type>\r\n\t\t\t\tvalve\r\n\t\t\t</pipe_type>\r\n\t\t\t<node_from>\r\n\t\t\t\t"+valve_node_from[i]+"\r\n\t\t\t</node_from>\r\n\t\t\t<node_to>\r\n\t\t\t\t"+valve_node_to[i]+"\r\n\t\t\t</node_to>\r\n\t\t\t<density>\r\n\t\t\t\t1000\r\n\t\t\t</density>\r\n\t\t\t<aref>\r\n\t\t\t\t3,1416e-04\r\n\t\t\t</aref>\r\n\t\t\t<mass_flow_rate>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</mass_flow_rate>\r\n\t\t\t<volume_flow_rate>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</volume_flow_rate>\r\n\t\t\t<velocity>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</velocity>\r\n\t\t\t<headloss>\r\n\t\t\t\t0.00\r\n\t\t\t</headloss>\r\n\t\t\t<head_loss_per_unit_length>\r\n\t\t\t\t0.00\r\n\t\t\t</head_loss_per_unit_length>\r\n\t\t\t<travel_time>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</travel_time>\r\n\t\t\t<concentration>\r\n\t\t\t\t0\r\n\t\t\t</concentration>\r\n\t\t\t<edge_spec>\r\n\t\t\t\t<valve>\r\n\t\t\t\t\t<position>\r\n\t\t\t\t\t\t"+to_string(0.)+"\r\n\t\t\t\t\t</position>\r\n\t\t\t\t\t<adzeta>\r\n\t\t\t\t\t\t"+to_string(0.)+"\r\n\t\t\t\t\t</adzeta>\r\n\t\t\t\t\t<curve>\r\n\t\t\t\t\t\t<id>\r\n\t\t\t\t\t\t\tcurve\r\n\t\t\t\t\t\t</id>\r\n\t\t\t\t\t\t<x_val>\r\n\t\t\t\t\t\t\te\r\n\t\t\t\t\t\t</x_val>\r\n\t\t\t\t\t\t<x_dim>\r\n\t\t\t\t\t\t\tpercent\r\n\t\t\t\t\t\t</x_dim>\r\n\t\t\t\t\t\t<y_val>\r\n\t\t\t\t\t\t\tdzeta\r\n\t\t\t\t\t\t</y_val>\r\n\t\t\t\t\t\t<y_dim>\r\n\t\t\t\t\t\t\tPa*s2/kg2\r\n\t\t\t\t\t\t</y_dim>\r\n\t\t\t\t\t\t<points>\r\n\t\t\t\t\t\t\t<point_x>\r\n\t\t\t\t\t\t\t\t0,0000\r\n\t\t\t\t\t\t\t</point_x>\r\n\t\t\t\t\t\t\t<point_y>\r\n\t\t\t\t\t\t\t\t0,0000\r\n\t\t\t\t\t\t\t</point_y>\r\n\t\t\t\t\t\t\t<point_x>\r\n\t\t\t\t\t\t\t\t100,0000\r\n\t\t\t\t\t\t\t</point_x>\r\n\t\t\t\t\t\t\t<point_y>\r\n\t\t\t\t\t\t\t\t1e10\r\n\t\t\t\t\t\t\t</point_y>\r\n\t\t\t\t\t\t</points>\r\n\t\t\t\t\t</curve>\r\n\t\t\t\t</valve>\r\n\t\t\t</edge_spec>\r\n\t\t</edge>\r\n";
  }

  // PUMPS
	for(int i=0; i<pump.size(); i++){
  	int color = -65536;
  	string pump_name = "PUMP_" + pump_node_from[i] + "_to_" + pump_node_to[i];
  	file_out << "\t\t<edge>\r\n\t\t\t<id>\r\n\t\t\t\t"+pump[i]+"\r\n\t\t\t</id>\r\n\t\t\t<name>\r\n\t\t\t\t"+pump_name+"\r\n\t\t\t</name>\r\n\t\t\t<color>\r\n\t\t\t\t"+to_string(color)+"\r\n\t\t\t</color>\r\n\t\t\t<counter>\r\n\t\t\t\t3\r\n\t\t\t</counter>\r\n\t\t\t<pipe_type>\r\n\t\t\t\tpump\r\n\t\t\t</pipe_type>\r\n\t\t\t<node_from>\r\n\t\t\t\t"+pump_node_from[i]+"\r\n\t\t\t</node_from>\r\n\t\t\t<node_to>\r\n\t\t\t\t"+pump_node_to[i]+"\r\n\t\t\t</node_to>\r\n\t\t\t<density>\r\n\t\t\t\t1000\r\n\t\t\t</density>\r\n\t\t\t<aref>\r\n\t\t\t\t1.0\r\n\t\t\t</aref>\r\n\t\t\t<mass_flow_rate>\r\n\t\t\t\t0.0\r\n\t\t\t</mass_flow_rate>\r\n\t\t\t<volume_flow_rate>\r\n\t\t\t\t0.0\r\n\t\t\t</volume_flow_rate>\r\n\t\t\t<velocity>\r\n\t\t\t\t0.0\r\n\t\t\t</velocity>\r\n\t\t\t<headloss>\r\n\t\t\t\t0.0\r\n\t\t\t</headloss>\r\n\t\t\t<head_loss_per_unit_length>\r\n\t\t\t\t0.0\r\n\t\t\t</head_loss_per_unit_length>\r\n\t\t\t<travel_time>\r\n\t\t\t\t0.0\r\n\t\t\t</travel_time>\r\n\t\t\t<concentration>\r\n\t\t\t\t0.0\r\n\t\t\t</concentration>\r\n\t\t\t<edge_spec>\r\n\t\t\t\t<pump>\r\n\t\t\t\t\t<curve>\r\n\t\t\t\t\t\t<id>\r\n\t\t\t\t\t\t\tcurve\r\n\t\t\t\t\t\t</id>\r\n\t\t\t\t\t\t<x_val>\r\n\t\t\t\t\t\t\tQ\r\n\t\t\t\t\t\t</x_val>\r\n\t\t\t\t\t\t<x_dim>\r\n\t\t\t\t\t\t\tm3/h\r\n\t\t\t\t\t\t</x_dim>\r\n\t\t\t\t\t\t<y_val>\r\n\t\t\t\t\t\t\tH\r\n\t\t\t\t\t\t</y_val>\r\n\t\t\t\t\t\t<y_dim>\r\n\t\t\t\t\t\t\tm\r\n\t\t\t\t\t\t</y_dim>\r\n\t\t\t\t\t\t<points>\r\n";
  	for(int j=0; j<pump_cv_x[i].size(); j++)
	  	file_out << "\t\t\t\t\t\t\t<point_x>\r\n\t\t\t\t\t\t\t\t"+to_string(pump_cv_x[i][j])+"\r\n\t\t\t\t\t\t\t</point_x>\r\n\t\t\t\t\t\t\t<point_y>\r\n\t\t\t\t\t\t\t\t"+to_string(pump_cv_y[i][j])+"\r\n\t\t\t\t\t\t\t</point_y>\r\n";

	  file_out << "\t\t\t\t\t\t</points>\r\n\t\t\t\t\t</curve>\r\n\t\t\t\t</pump>\r\n\t\t\t</edge_spec>\r\n\t\t\t<user1>\r\n\t\t\t\t0.0\r\n\t\t\t</user1>\r\n\t\t\t<user2>\r\n\t\t\t\t0.0\r\n\t\t\t</user2>\r\n\t\t</edge>\r\n";
	}

	for(int i=0; i<pool.size(); i++){
  	string pool_name = "POOL_" + pool_node_from[i] + "_to_" + pool_node_to[i];
  	int color = -16776961;
		file_out << "\t\t<edge>\r\n\t\t\t<id>\r\n\t\t\t\t"+pool[i]+"\r\n\t\t\t</id>\r\n\t\t\t<name>\r\n\t\t\t\t"+pool_name+"\r\n\t\t\t</name>\r\n\t\t\t<color>\r\n\t\t\t\t"+to_string(color)+"\r\n\t\t\t</color>\r\n\t\t\t<counter>\r\n\t\t\t\t1\r\n\t\t\t</counter>\r\n\t\t\t<pipe_type>\r\n\t\t\t\tpool\r\n\t\t\t</pipe_type>\r\n\t\t\t<node_from>\r\n\t\t\t\t"+pool_node_from[i]+"\r\n\t\t\t</node_from>\r\n\t\t\t<node_to>\r\n\t\t\t\t"+pool_node_to[i]+"\r\n\t\t\t</node_to>\r\n\t\t\t<density>\r\n\t\t\t\t1000\r\n\t\t\t</density>\r\n\t\t\t<aref>\r\n\t\t\t\t"+to_string(pool_aref[i])+"\r\n\t\t\t</aref>\r\n\t\t\t<mass_flow_rate>\r\n\t\t\t\t0.0\r\n\t\t\t</mass_flow_rate>\r\n\t\t\t<volume_flow_rate>\r\n\t\t\t\t0.0\r\n\t\t\t</volume_flow_rate>\r\n\t\t\t<velocity>\r\n\t\t\t\t0.0\r\n\t\t\t</velocity>\r\n\t\t\t<headloss>\r\n\t\t\t\t0.0\r\n\t\t\t</headloss>\r\n\t\t\t<head_loss_per_unit_length>\r\n\t\t\t\t0.0\r\n\t\t\t</head_loss_per_unit_length>\r\n\t\t\t<travel_time>\r\n\t\t\t\t0.0\r\n\t\t\t</travel_time>\r\n\t\t\t<concentration>\r\n\t\t\t\t0,0000\r\n\t\t\t</concentration>\r\n\t\t\t<edge_spec>\r\n\t\t\t\t<pool>\r\n\t\t\t\t\t<water_level>\r\n\t\t\t\t\t\t"+to_string(pool_watlev[i])+"\r\n\t\t\t\t\t</water_level>\r\n\t\t\t\t\t<bottom_level>\r\n\t\t\t\t\t\t"+to_string(pool_botlev[i])+"\r\n\t\t\t\t\t</bottom_level>\r\n\t\t\t\t</pool>\r\n\t\t\t</edge_spec>\r\n\t\t\t<user1>\r\n\t\t\t\t0.0\r\n\t\t\t</user1>\r\n\t\t\t<user2>\r\n\t\t\t\t0.00000e+00\r\n\t\t\t</user2>\r\n\t\t</edge>\r\n";
	}

	for(int i=0; i<pres_id.size(); i++){
		string pres_name = "PRES_" + pres_node_from[i] + "_to_" + pres_node_to[i];
  	int color = -16776961;
		file_out << "\t\t<edge>\r\n\t\t\t<id>\r\n\t\t\t\t"+pres_id[i]+"\r\n\t\t\t</id>\r\n\t\t\t<name>\r\n\t\t\t\t"+pres_name+"\r\n\t\t\t</name>\r\n\t\t\t<color>\r\n\t\t\t\t"+to_string(color)+"\r\n\t\t\t</color>\r\n\t\t\t<counter>\r\n\t\t\t\t2\r\n\t\t\t</counter>\r\n\t\t\t<pipe_type>\r\n\t\t\t\tpress\r\n\t\t\t</pipe_type>\r\n\t\t\t<node_from>\r\n\t\t\t\t"+pres_node_from[i]+"\r\n\t\t\t</node_from>\r\n\t\t\t<node_to>\r\n\t\t\t\t"+pres_node_to[i]+"\r\n\t\t\t</node_to>\r\n\t\t\t<density>\r\n\t\t\t\t1000\r\n\t\t\t</density>\r\n\t\t\t<aref>\r\n\t\t\t\t1.0\r\n\t\t\t</aref>\r\n\t\t\t<mass_flow_rate>\r\n\t\t\t\t0.0\r\n\t\t\t</mass_flow_rate>\r\n\t\t\t<volume_flow_rate>\r\n\t\t\t\t0.0\r\n\t\t\t</volume_flow_rate>\r\n\t\t\t<velocity>\r\n\t\t\t\t0.0\r\n\t\t\t</velocity>\r\n\t\t\t<headloss>\r\n\t\t\t\t0.0\r\n\t\t\t</headloss>\r\n\t\t\t<head_loss_per_unit_length>\r\n\t\t\t\t0.0\r\n\t\t\t</head_loss_per_unit_length>\r\n\t\t\t<travel_time>\r\n\t\t\t\t0.0\r\n\t\t\t</travel_time>\r\n\t\t\t<concentration>\r\n\t\t\t\t0\r\n\t\t\t</concentration>\r\n\t\t\t<edge_spec>\r\n\t\t\t\t<press>\r\n\t\t\t\t\t<pressure>\r\n\t\t\t\t\t\t"+to_string(pres_p[i])+"\r\n\t\t\t\t\t</pressure>\r\n\t\t\t\t</press>\r\n\t\t\t</edge_spec>\r\n\t\t\t<user1>\r\n\t\t\t\t0.0\r\n\t\t\t</user1>\r\n\t\t\t<user2>\r\n\t\t\t\t0.0\r\n\t\t\t</user2>\r\n\t\t</edge>\r\n";
	}

  file_out << "\t</edges>\r\n";

  file_out << "</staci>\r\n";
  file_out.close();
  cout << "\n[*]Writing to file: OK" << endl;
  cout << endl << "[*]File has been created successfully, name: " << case_name + ".spr" << endl;

  cout <<"\n    =========\n    |SUMMARY|\n    =========\n";
  cout << "\nNodes:         " << node.size(); 
  cout << "\nPipes:         " << pipe.size(); 
  cout << "\nValves:        " << valve.size();
  cout << "\nPumps:         " << pump.size();
  cout << "\nPres:          " << pres_id.size();
  cout << "\nPool:          " << pool.size();

  cout << "\n\nOverall CPU time:     " << time(0)-ido << " s";
  cout << "\nWarnings: " << warning_counter << endl;

  cout << endl << endl;
	return 0;
}

vector<string> line2sv(string line){
	string s="";
	vector<string> sv;
	for(string::iterator i=line.begin(); i!=line.end(); i++){
		if(*i!=' ' && *i!='\t' && *i!=';' && *i!='\n' && *i!='\r')
			s += *i;
		else if(s.length()>0){
			sv.push_back(s);
			s="";
		}
	}
	if(s!="")
		sv.push_back(s);
	return sv;
}

double VectorMax(const vector<double> &v){
  double max;
  if(v.size()>0){
    max=v[0];
    for(int i=0; i<v.size(); i++)
      if(v[i]>max)
        max = v[i];
  }else{
    cout << "\nERROR! VectorMax: vector is empty! Exiting..." << endl;
    exit(-1);
  }
  return max;
}

double VectorMin(const vector<double> &v){
  double max;
  if(v.size()>0){
    max=v[0];
    for(int i=0; i<v.size(); i++)
      if(v[i]<max)
        max = v[i];
  }else{
    cout << "\nERROR! VectorMax: vector is empty! Exiting..." << endl;
    exit(-1);
  }
  return max;
}