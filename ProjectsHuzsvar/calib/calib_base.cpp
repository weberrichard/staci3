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

int main(int argc, char* argv[]){

  // Name of containing folder of staci file
  string case_name = "lin_alap_9";
  cout << endl << "Case: " << case_name << endl;
  string sens_par_name = "friction_coeff";

  // Creating the case
  wds = new SVDCalibration(case_name + ".spr");
  wds->Set_debug_level(1);
  wds->build_system();
  wds->ini();
  int n_cspok = wds->cspok.size();
  int n_agelemek = wds->agelemek.size();
  
  // Setting the real friction and the estimations
  vector<double> fric_real, fric_est;
  cout << "\nfric_real:\n";
  for(int i=0; i<n_agelemek; i++)
    if(wds->agelemek.at(i)->GetType() == "Cso"){
      fric_real.push_back(wds->agelemek.at(i)->Get_dprop("erdesseg"));
      fric_est.push_back(fric_real.at(i)*1.5);
      cout << "f at " << i << "\t" << fric_real.at(i) << endl;
    }

  // Setting the demand perturbation values and mode
  vector<double> d_noise; d_noise.push_back(1.0); d_noise.push_back(0.1);
  string d_noise_mode = "None";

  // Evaluating the "measurements"
  wds->GenMeas(fric_real,d_noise,d_noise_mode);

  // Printing the loaded measurement data to consol
  vector<double> d_meas_sum = wds->Get_d_meas_sum();
  vector<vector<double> > p_meas = wds->Get_p_meas();
  cout << "\nMeasured values:";
  printf("\n a_fogy | p ");
  for(int i=0; i<d_meas_sum.size(); i++){
    printf("\n %6.3f | ",d_meas_sum.at(i));
    for(int j=0; j<p_meas[0].size(); j++)
      printf(" %6.3f",p_meas.at(i).at(j));
  }

  // Perfroming the calibration
  double tol = .9;
  int konv = wds->Calibrate(fric_est, tol);

  // Printing the frictions after, before the calibration and the "real"
  cout << "\nfrics after Calibration: \n";
  printf(" # |    calib   |    real    |   estima   |\n");
  printf("---+------------+------------+------------|\n");
  for(int i=0; i<n_agelemek; i++)
    if(wds->agelemek.at(i)->GetType() == "Cso")
      printf("%2i | %10.3f | %10.3f | %10.3f |\n", i, wds->agelemek.at(i)->Get_dprop("erdesseg"), fric_real.at(i),fric_est.at(i));

  // Calculating the norm of the differences
  double e_fric=0.;
  for(int i=0; i<fric_real.size(); i++)
    e_fric += pow(fric_real[i] - fric_est.at(i), 2.);
  e_fric = pow(e_fric/(n_agelemek-1),.5);
  cout << "\n\n Fric difference at estimation: " << e_fric;

  e_fric=0.;
  for(int i=0; i<fric_real.size(); i++)
    e_fric += pow(wds->agelemek.at(i)->Get_dprop("erdesseg") - fric_real.at(i), 2.);
  e_fric = pow(e_fric/(n_agelemek-1),.5);
  cout << "\n Fric difference after calib:   " << e_fric;

  // Printing whether Calibration converged
  printf("\n\n Konv: ");
  switch(konv) {
    case 0 : printf("NO :(\n"); break;
    case 1 : printf("FULL :)\n"); break;
    case 2 : printf("LOCAL :|\n"); break;
  }

  cout << endl;
  return 0;
}