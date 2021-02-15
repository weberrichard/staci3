/*===================================================================*\
                               Statistic
                            ---------------

  Statistical calculations, e.g. average, deviation, normal dist.

  Can be used independently from staci3.
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef STATISTIC_H
#define STATISTIC_H

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <iostream>

#include </usr/include/eigen3/Eigen/Eigen>

using namespace std;
using namespace Eigen;

// WR average of a vector
double average(const vector<double> &x);
double average(const vector<int> &x);

// WR maximum of a vector in absolute value
double absoluteMax(const vector<double> &x);

// WR standard deviation of a vector
double standardDeviation(const vector<double> &x);

// Normal distribution double with average of m and deviation of szig
double normalDistribution(double m, double szig);

// Uniform distribution double betweein min and max
double uniformDistribution(double min, double max);

// Chi squrad test to x vector with expected value of m and deviation of szig, only normal distribution
bool chiSquaredTest(const vector<double> &x, double m, double szig);

// Discrete homogenity test with chi square test where f1 f2 are the frequencies in case of 95%
bool discreteHomogenityChi2(vector<int> f1, vector<int> f2);

// 2-sample U test wit averages, standardDeviations, sizes and p is the significanse
bool twoSampleUTest(double m1, double s1, int n1, double m2, double s2, int n2, double p);

// Approximation of the inverse function of the error function //really good approximation
double erfinv(double x);

// Sample correlation coefficient between x and y, x.size() must be equal to y.size()
double correlCoefficient(const vector<double> &x, const vector<double> &y);

// Spearman's rho
double spearmanCoefficient(const vector<double> &x, const vector<double> &y);

// Simple Grubbs test, returns a vector without the element under x_mean-3*x_dev and abov x_mean+3*x_dev
vector<double> grubbsTest(const vector<double> &x);

/// Linear interpolation using x-y points between the two closest points	
double linearInterpolation(const vector<double> &x, const vector<double> &y, double xp);

// Interpolation using x-y points, with two points linear, with three or above parabole with the closest (in x) points
double interpolate(vector<double> x, vector<double> y, double xp);
	
// Finding the unique elements in an integer vector
vector<int> unique(const vector<int> &x);

// Selection sorting algorithm, pretty slow for large vectors 
void selectionSort(vector<double> &x);

// Dynamic sorting algorithm from index "low" until index "high" and following the changes and saving in idx
void quickSortWithIndicies(vector<double> &x, vector<int> &idx, int low, int high);
void quickSortWithIndicies(vector<int> &x, vector<int> &idx, int low, int high);

// Dynamic sorting algorithm from index "low" until index "high"
void quickSort(vector<double> &x, vector<int> &idx, int low, int high);
void quickSort(vector<int> &x, vector<int> &idx, int low, int high);

// Function for quickSort
int partition(vector<double> &x, vector<int> &idx, int low, int high);
int partition(vector<int> &x, vector<int> &idx, int low, int high);

// Changing two elements in vector double/integer
void swap(vector<double> &x, int i, int j);
void swap(vector<int> &x, int i, int j);

// Finding the min element and index in vector
double min(vector<double> x, int &idx);
double min(vector<double> x); // without index tracing

// Finding the max element and index in vector
double max(vector<double> x, int &idx);
int max(vector<int> x, int &idx);
double max(vector<double> x); // without index tracing
int max(vector<int> x); // without index tracing
	
// Fitting polynomial function on points with X-Y coordinates
// returns the coefficients of the polynomail in a0 + a1*x^1 + a2*x^2 + ... order
VectorXd leastSquaresPolynomial(const VectorXd &x, const VectorXd &y, int &order);

// Powering Eigen VectorXd by elements to order
VectorXd eigenVectorXdPow(const VectorXd &x, int order);

// Finding the n number of max values in an Eigen VectorXd, also tracking the indicies
// if av is true, it takes the absolute value of f
vector<double> findMaxValues(const VectorXd &f, vector<int> &vidx, int n, bool av);

#endif