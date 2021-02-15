#include <iostream>
#include "/home/rweber/Dropbox/0_PhD/staci/SVDCalibration.h"
//#include "/home/rweber/Dropbox/0_PhD/staci/Staci.h"
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

SVDCalibration *wds;

vector<string> LoadMeasIds(string file_name);

int main(int argc, char* argv[]){

  double tol = 1.0;
  // Name of containing folder of staci file
  string case_name;
  if(argc == 1){
      //case_name = "VIZ-SOPTVR-B-input_mod";
      //case_name = "Anytown_DC1";
      //case_name = "C-town_mod_10";
      //case_name = "VIZ-SOPTVR-M-input_mod";
      //case_name = "VIZ-SOPTVR-Q-input_mod";
      //case_name = "lin_tank_10";
      case_name = "lin_alap_9";
      //case_name = "racs_alap_9";
      //case_name = "kor_alap_6";
  }else if(argc == 2){
      case_name = argv[1];
  }
  cout << endl << "Case: " << case_name << endl;
  cout << setprecision(6) << fixed;
  srand((unsigned int) time(0));
  //string case_folder = "/home/rweber/0_PhD/Halozatok/sopron_halozatok/";
  //string case_folder = "/home/rweber/0_PhD/Halozatok/general_alap/";
  //string case_folder = "/home/rweber/0_PhD/Halozatok/anytown/";
  //string case_folder = "/home/rweber/0_PhD/Halozatok/";
  string case_folder = "";

  string sens_par_name = "friction_coeff";
  // Creating the case
  //wds = new Staci(case_folder + case_name + ".spr");
  wds = new SVDCalibration(case_folder + case_name + ".spr");
  wds->Set_debug_level(1);
  wds->build_system();
  wds->ini();
  int n_cspok = wds->cspok.size();
  int n_agelemek = wds->agelemek.size();
  int n_cso=0;
  for(int i=0; i<n_agelemek; i++)
    if(wds->agelemek.at(i)->GetType() == "Cso")
      n_cso++;
  
  vector<double> fric_real, p_real;
  cout << "\nfric_real:\n";
  for(int i=0; i<n_agelemek; i++)
    if(wds->agelemek.at(i)->GetType() == "Cso"){
      //wds->agelemek.at(i)->Set_dprop("erdesseg",-0.03);
      //wds->agelemek.at(i)->Set_friction_model("HW");
      //wds->agelemek.at(i)->Set_dprop("erdesseg",100.*(.1*sin(i)+1));
      //wds->agelemek.at(i)->Set_dprop("erdesseg",-0.04*(.1*sin(i)+1));
      fric_real.push_back(wds->agelemek.at(i)->Get_dprop("erdesseg"));
      cout << "f at " << i << "\t" << fric_real.at(i) << endl;
    }

  cout << "\ncons:";
  double cons_nom = wds->get_sum_of_pos_consumption();
  if(case_folder == "/home/rweber/0_PhD/Halozatok/general_alap/" || case_folder == ""){
    wds->cspok.at(n_cspok-1)->Set_fogy(-cons_nom);
    //wds->cspok.at(5)->Set_fogy(cons_nom/2.);
    for(int i=0; i<n_cspok; i++){
      if(wds->cspok.at(i)->Get_fogy()>0.)
        wds->cspok.at(i)->Set_fogy(cons_nom/(double)(n_cspok-1.0));
      cout << endl << wds->cspok.at(i)->Get_fogy();
    }
  }

  wds->solve_system();
  cout << "\nPressure dist:";
  for(int i=0; i<n_cspok; i++){
    cout << endl << wds->cspok.at(i)->Get_p();
    p_real.push_back(wds->cspok.at(i)->Get_p());
  }
  cout << "\nFlow meters:";
  for(int i=0; i<n_agelemek;i++)
    if(wds->agelemek.at(i)->GetType() == "Cso")
      cout << endl << wds->agelemek.at(i)->Get_v();

  vector<int> fogy_meas_idx;
  //fogy_meas_idx.push_back(5);
  //fogy_meas_idx.push_back(8);
  vector<string> fogy_meas_id;
  for(int i=0; i<fogy_meas_idx.size(); i++)
    fogy_meas_id.push_back(wds->cspok.at(fogy_meas_idx.at(i))->Get_nev());
  cout << "\nfogy_meas_id.size: " << fogy_meas_idx.size() << endl << endl << endl;

  vector<double> fogy_pert; fogy_pert.push_back(1.0); fogy_pert.push_back(0.1);
  string fogy_pert_mode = "None";

  //vector<vector<double> > p_meas_all;
  //vector<vector<double> > p_meas = wds->GenMeas(fric_real,fogy_sum,meas_node_id,fogy_pert,fogy_pert_mode,fogy_meas_id,fogy_meas,p_meas_all);
  //wds->Set_p_meas_idx(meas_node_id);
  //wds->Set_d_meas_idx(fogy_meas_id);
  //wds->Set_d_meas_sum(fogy_sum);
  wds->GenMeas(fric_real,fogy_pert,fogy_pert_mode);
  vector<vector<double> > p_meas = wds->Get_p_meas();
  vector<vector<double> > p_meas_all = wds->Get_p_meas_all();
  cout << "\np_meas: " << p_meas.size() << endl;
  cout << "\np_meas_all: " << p_meas_all.size() << 'x' << p_meas_all[0].size() << endl;

  vector<double> fogy_sum = wds->Get_d_meas_sum();
  cout << "\nMeasured values:";
  printf("\n a_fogy | p ");
  for(int i=0; i<fogy_sum.size(); i++){
    printf("\n %6.3f | ",fogy_sum.at(i));
    for(int j=0; j<p_meas[i].size(); j++)
      printf(" %6.3f",p_meas.at(i).at(j));
  }

  //for(int i=0; i<fogy_sum.size(); i++)
  //  fogy_sum.at(i) = fogy_sum.at(i)*wds->stat_tools->NormalDist(1.0,0.02);

  //cout << "\nfogy_meas size: " << fogy_meas.size() << "x" << fogy_meas.at(0).size();
  cout << "\nfogy_sum size: " << fogy_sum.size();
  cout << "\nfogy_meas_id size: " << fogy_meas_id.size() << endl;

  //double tol = -1e3;
  vector<double> fric_est;
  for(int i=0; i<n_agelemek; i++)
    if(wds->agelemek.at(i)->GetType() == "Cso"){ 
      //fric_est.push_back(-0.02);
      //fric_est.push_back(fric_real.at(i)*wds->stat_tools->NormalDist(1.0,0.2));
      fric_est.push_back(fric_real.at(i)*1.5);
      wds->agelemek.at(i)->Set_dprop("erdesseg",fric_est.at(i));
    }
  wds->solve_system_sens(sens_par_name,0.);

  //vector<vector<double> > p_est_all;
  //vector<vector<double> > p_est = wds->GenMeas(fric_est,fogy_sum,meas_node_id,fogy_pert,fogy_pert_mode,fogy_meas_id,fogy_meas,p_est_all);

  //cout << "\nSM_row.size: " << wds->SM_row_sum_Pressures.size();
  //cout << "\nSM_col.size: " << wds->SM_col_sum_Pressures.size();

  int konv = wds->Calibrate(fric_est, tol);
  //int konv = wds->Staci::Calibration(fric_est, p_meas, meas_node_id, fogy_sum, fogy_meas_id, fogy_meas, tol);

  wds->GenMeas(fric_est,fogy_pert,fogy_pert_mode);
  vector<vector<double> > p_est = wds->Get_p_meas();
  vector<vector<double> > p_est_all = wds->Get_p_meas_all();

  //cout << "\nfrics after Calibration: \n";
  //printf(" # |    calib   |    real    |   estima   |\n");
  //printf("---+------------+------------+------------|\n");
  //for(int i=0; i<n_agelemek; i++)
  //  if(wds->agelemek.at(i)->GetType() == "Cso")
  //    printf("%2i | %10.3f | %10.3f | %10.3f |\n", i, wds->agelemek.at(i)->Get_dprop("erdesseg"), fric_real.at(i),fric_est.at(i));

  //cout << "\npres after Calibration: \n";
  //printf(" # |    calib   |    real    |   estima   |\n");
  //printf("---+------------+------------+------------|\n");
  //for(int i=0; i<n_cspok; i++)
  //  printf("%2i | %10.6f | %10.6f | %10.6f |\n", i, wds->cspok[i]->Get_p(), p_real.at(i),p_est_all.at(i).at(0));

  double e_fric=0.;
  for(int i=0; i<fric_real.size(); i++)
    e_fric += pow(fric_real[i] - fric_est.at(i), 2.);
  e_fric = pow(e_fric/n_cso,.5);
  cout << "\n\n Fric difference at estimation: " << e_fric;

  e_fric=0.;
  for(int i=0; i<fric_real.size(); i++)
    e_fric += pow(wds->agelemek.at(i)->Get_dprop("erdesseg") - fric_real.at(i), 2.);
  e_fric = pow(e_fric/n_cso,.5);
  cout << "\n Fric difference after calib:   " << e_fric;

  vector<double> fric_calib;
  for(int i=0; i<n_agelemek; i++)
    if(wds->agelemek[i]->GetType() == "Cso")
      fric_calib.push_back(wds->agelemek[i]->Get_dprop("erdesseg"));

  //vector<vector<double> > p_calib_all;
  //vector<vector<double> > p_calib = wds->GenMeas(fric_calib,fogy_sum,meas_node_id,fogy_pert,fogy_pert_mode,fogy_meas_id,fogy_meas,p_calib_all);
  wds->GenMeas(fric_calib,fogy_pert,fogy_pert_mode);
  vector<vector<double> > p_calib = wds->Get_p_meas();
  vector<vector<double> > p_calib_all = wds->Get_p_meas_all();


  double e_p=0.;
  for(int i=0; i<fogy_sum.size(); i++)
    for(int j=0; j<n_cspok; j++)
      e_p += pow(p_est_all[i][j] - p_meas_all[i][j],2.);
  e_p = pow(e_p/(double)n_cspok,.5);
  printf("\n\n Pres difference at estimation: %8.5e",e_p);

  e_p=0.;
  for(int i=0; i<fogy_sum.size(); i++)
    for(int j=0; j<n_cspok; j++)
      e_p += pow(p_calib_all[i][j] - p_meas_all[i][j],2.);
  e_p = pow(e_p/(double)n_cspok,.5);
  printf("\n Pres difference after calib:   %8.5e",e_p); 

  printf("\n\n Konv: ");
  switch(konv) {
    case 0 : printf("NO :(\n"); break;
    case 1 : printf("FULL :)\n"); break;
    case 2 : printf("LOCAL :|\n"); break;
}

  //string new_case_name = case_name +"_c";
  //wds->new_def_file = case_folder + new_case_name + ".spr";
  //wds->copy_file(wds->get_def_file(), wds->new_def_file);
  //wds->set_res_file(wds->new_def_file);
  //wds->save_results(true);

  cout << endl;
  return 0;
}

vector<string> LoadMeasIds(string file_name){
  ifstream file;
  file.open(file_name);
  string line;
  vector<string> names;
  if(file.is_open()){
  while(getline(file,line)){
    if(line[0] != '/'){
      string s;
      for(string::iterator j=line.begin(); j!=line.end(); j++){
        if(*j!=' ')
          s +=*j;
      }
    names.push_back(s);
    }
  }
  }else
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl << "LoadMeasIds(), File is not open when calling LoadMeasIds() function!!!" << endl;
  file.close();

  return names;
}