/*===================================================================*\
                              QualitySolver
                            -----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef QUALITYSOLVER_H
#define QUALITYSOLVER_H

#include "WaterAge.h"
#include "chlorine.h"
#include "SeriesHydraulics.h"
#include "Biofilm.h"

class QualitySolver : public SeriesHydraulics
{
public:
  QualitySolver(string spr_filename);
  ~QualitySolver();

  // setting adaptive parameters for dt handling
  void setAdaptiveParameters(double e_max, double e_min, double dt_inc);
  // setting spatial parameters
  void setAdaptiveParameters(double e_min, double e_max, int modAttempts, double dt_inc, double timeMax);

  int solveQuality(double t, double h, double tmax, double tolerance, double DX, string whichSolver);

  vector<int> generateNumberOfIntervals(double DX);

  vector< vector< vector<double> > > initializeWithStartingCondition(vector<int> numberOfIntervals);

  vector< vector<double> > updateInitializeTimestep(vector< vector<double> > out, double initTimestep);

  vector< vector<double> > updateInitializeTimestep(vector< vector<double> > out, double initTimestep, vector<int> numberOfIntervals);

  void updateNodeVariables(vector < vector< vector<double> > > &x, string WhichSolver, double h);

  void solveQuality(string WhichSolver);

  vector<double> source(double t, vector<double> x_actual, vector<double> x_upwind, double flow, double DX, string Which);

  double getSmallestTimestep(vector< vector<double> > input);

  double time_start = 0.;
  double time_max = 2*3600;

protected:

  vector<double> Runge_Kutta(vector< vector<double> >& x, double t, double h, vector<double> x_previous, double flow, double DX, string WhichSolver);
  // general ode45 solver with S at the right side, i.e. dx/dt = S
  

private:
  WaterAge *wa = new WaterAge();

  chlorine *ch = new chlorine();

  biofilm *bf = new biofilm();

  // for adaptive timestep handling
  double minScaleFactor = 0.125;
  double maxScaleFactor = 25.0;
  int SourceDimension = 1;
  int attempts = 12;
  string solverType = "waterAge";
  // actual timestep
  double dt;
  double step_time = 300;
  double length_max = 0.;
  double tolerance = 0.1;
  // spacial division
  double step_X = 10.;

  vector<double> startingCondition;
  vector<string> ParameterList;
};

#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     This Runge-Kutta-Prince-Dormand method is an adaptive procedure for    //
//     approximating the solution of the differential equation x'(t) = f(t,x) //
//     with initial condition x(t0) = c.  This implementation evaluates       //
//     f(t,x) six times per step using embedded fourth order and fifth order  //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        x[i+1] = x[i] +  h * ( 31 / 540 * k1 + 190 / 297 * k3               //
//                           - 145 / 108 * k4 + 351 / 220 * k5 + 1/20 * k6 )  //
//     where                                                                  //
//     k1 = f( t[i],x[i] ),                                                   //
//     k2 = f( t[i]+h/5, x[i] + h*k1/5 ),                                     //
//     k3 = f( t[i]+3h/10, x[i]+(h/40)*(3 k1 + 9 k2) ),                       //
//     k4 = f( t[i]+3h/5, x[i]+(h/10)*(3 k1 - 9 k2 + 12 k3) ),                //
//     k5 = f( t[i]+2h/3, x[i]+(h/729)*(226 k1 - 675 k2 + 880 k3 + 55 k4) )   //
//     k6 = f( t[i]+h, x[i]+(h/2970)*(-1991 k1 + 7425 k2 - 2660 k3            //
//                                                  - 10010 k4 + 10206 k5) )  //
//     t[i+1] = t[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( 77 k1 - 400 k3 + 1925 k4 - 1701 k5 + 99 k6 ) / 2520       //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * x[i] / [err * (tmax - t[0])] | ^ 1/4     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
// int Embedded_Prince_Dormand_v1_4_5( double (*f)(double, double),           //
//       double x[], double t, double h, double tmax, double h_next,          //
//                                                        double tolerance )  //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation x'=f(t,x) with the      //
//     initial condition x(t) = x[0].  The value at tmax is returned in x[1]. //
//     The function returns 0 if successful or -1 if it fails.                //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (t,x) of //
//                integral curve of the differential equation x' = f(t,x)     //
//                which passes through the point (t0,x0) corresponding to the //
//                initial condition x(t0) = x0.                               //
//     double x[] On input x[0] is the initial value of x at t, on output     //
//                x[1] is the solution at tmax.                               //
//     double t   The initial value of t.                                     //
//     double h   Initial step size.                                          //
//     double tmax The endpoint of t.                                         //
//     double h_next   A pointer to the estimated step size for successive    //
//                      calls to Runge_Kutta_Prince_Dormand_v1_4_5.           //
//     double tolerance The tolerance of x(tmax), i.e. a solution is sought   //
//                so that the relative error < tolerance.                     //
//                                                                            //
//  Return Values:                                                            //
//     0   The solution of x' = f(t,x) from t to tmax is stored x[1] and      //
//        h_next has the value to the next size to try.                       //
//    -1   The solution of x' = f(t,x) from t to tmax failed.                 //
//    -2   Failed because either tmax < t or the step size h <= 0.            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Prince-Dormand's embedded 4th and 5th order methods  //
//     to approximate the solution of the differential equation y'=f(x,y)     //
//     with the initial condition y = y[0] at x = x0.  The value at x + h is  //
//     returned in y[1].  The function returns err / h ( the absolute error   //
//     per step size ).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y[0]).                   //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at x + h.                              //
//     double x   Initial value of x.                                         //
//     double h   Step size                                                   //
//                                                                            //
//  Return Values:                                                            //
//     This routine returns the err / h.  The solution of y(x) at x + h is    //
//     returned in y[1].                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////