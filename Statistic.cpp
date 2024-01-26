#include "Statistic.h"

using namespace std;
using namespace Eigen;

//--------------------------------------------------------------
double average(const vector<double> &x)
{
  double a=0.;
  for(int i=0; i<x.size() ;i++)
    a += x[i];
  return a/(double)x.size();
}
double average(const vector<int> &x)
{
  double a=0.;
  for(int i=0; i<x.size() ;i++)
    a += x[i];
  return a/(double)x.size();
}

//--------------------------------------------------------------
double sum(const vector<double> &x)
{
  double s=0.;
  for(int i=0; i<x.size(); i++)
    s += x[i];
  return s;
}

//--------------------------------------------------------------
double absoluteMax(const vector<double> &x)
{
  double m=0.;
  for(int i=0; i<x.size(); i++)
    if((abs(x[i])>m))
      m = abs(x[i]);
  return m;
}

//--------------------------------------------------------------
double standardDeviation(const vector<double> &x)
{
  double s=0.;
  double m = average(x);
  for(int i=0; i<x.size(); i++)
    s += pow(x[i]-m,2.);
  return pow(s/((double)x.size()-1),.5);
}

//--------------------------------------------------------------
double normalDistribution(double m, double szig)
{
  random_device rd;
  default_random_engine generator( rd() );
  normal_distribution<double> dist(m,szig);
  return dist(generator);
}

//--------------------------------------------------------------
double uniformDistribution(double min, double max)
{
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}

//--------------------------------------------------------------
bool chiSquaredTest(const vector<double> &x, double m, double szig)
{
  bool r=false;
  int no_intervals = 6;
  int n_data = x.size();
  vector<double> limits;
  double sqrt2 = sqrt(2.);
  // boundaries between the intervals using standard normal distribution function: FI^-1(0.1),FI^-1(0.2)...FI^-1(0.9)
  for(int i=0; i<no_intervals-1;i++){
    limits.push_back( -sqrt2*erfinv(1-2*(i+1)/(double)no_intervals) );
    //cout << "limits: " << limits.at(i) << "\t FI: " << .5*erfc(-limits.at(i)/sqrt2) << endl;
  }

  vector<int> counters(no_intervals);
  for(int k=0; k<n_data; k++){
    double x_st = (x.at(k)-m)/szig; // Standardizing the x vector with m and szig
      for(int i=0; i<no_intervals-2; i++){
          if(limits.at(i)<x_st && limits.at(i+1)>x_st) 
              counters.at(i+1)++;
      }
      if(x_st<limits.at(0))
              counters.at(0)++;
      if(x_st>limits.at(no_intervals-2))
              counters.at(no_intervals-1)++;
  }

    double khi=0;
    for(int i=0; i<no_intervals; i++){
      //cout << "counters: " << counters.at(i);
      khi += pow(counters.at(i)-n_data/no_intervals,2)/(n_data/no_intervals);
      //cout << "\t khi: " << khi << endl;
    }
    /*p: probability(0.95,0.98,0.99)
      k: degrees of freedom (no. of intervals(6,10) - 1) TODO: using vector,vector; then search for p and k, save i and j, if not found ERROR
      |------------------------------|
      |1-p|  0.05  |  0.02  |  0.01  |
      |k -+--------+--------+--------|
      |5  | 11.070 | 13.388 | 15.086 |
      |9  | 16.919 | 19.679 | 21.666 |
      |------------------------------|*/
    if(khi<=11.070)
      r = true;

  return r;
}

//--------------------------------------------------------------
bool discreteHomogenityChi2(vector<int> f1, vector<int> f2)
{
  int n1=0,n2=0;
  vector<double> r1(f1.size()), r2(f2.size());
  for(int i=0; i<f1.size(); i++)
    n1 += f1[i];
  for(int i=0; i<f1.size(); i++)
    r1[i] = (double)f1[i] / n1;

  for(int i=0; i<f2.size(); i++)
    n2 += f2[i];
  for(int i=0; i<f2.size(); i++)
    r2[i] = (double)f2[i] / n2;

  int m;
  if(f1.size()<f2.size())
  {
    m = f2.size();
    r1.resize(m,0.);
    f1.resize(m,0);
  }
  else
  {
    m = f1.size();
    r2.resize(m,0.);
    f2.resize(m,0);
  }

  double x2obs=0.;
  int r=m;
  for(int i=0; i<m; i++)
  {
    if(f1[i] + f2[i] != 0)
      x2obs += pow((r1[i]-r2[i]),2) / (f1[i] + f2[i]);
    else
      r--;
  }
  x2obs *= n1;
  x2obs *= n2;

  // critical values from 0 to 20 95%
  //vector<double> x2crit{0.,0.,3.841458821,5.991464547,7.814727903,9.487729037,11.07049769,12.59158724,14.06714045,15.50731306,16.9189776,18.30703805,19.67513757,21.02606982,22.36203249,23.6847913,24.99579014,26.2962276,27.58711164,28.86929943,30.14352721};
  // critical values from 0 to 20 98%
  vector<double> x2crit{0.,0.,5.411894431,7.824046011,9.837409311,11.6678434,13.3882226,15.03320775,16.62242187,18.16823076,19.67901609,21.16076754,22.61794081,24.05395669,25.47150914,26.87276464,28.25949634,29.63317731,30.99504721,32.34616093,33.68742507,};

  if(x2crit[r] == 0)
    cout << endl << "!WARNING! x2crit = 0, r (" << r << ") is not handeld." << endl;

  if(x2obs < x2crit[r])
    return true;
  else
    return false;
}


//--------------------------------------------------------------
bool twoSampleUTest(double m1, double s1, int n1, double m2, double s2, int n2, double p){
  bool test_ok = false;
  double u_akt, u_krit;

  u_akt = (m1-m2)/pow(pow(s1,2.)/(double)n1 + pow(s2,2.)/(double)n2,.5);
  u_krit = sqrt(2)*erfinv(p);
  if(abs(u_akt)<u_krit)
    test_ok = true;
  return test_ok;
}

//--------------------------------------------------------------
double erfinv(double x)
{
  double a3 = -0.140543331, a2 = 0.914624893, a1 = -1.645349621, a0 = 0.886226899;
  double b4 = 0.012229801, b3 = -0.329097515, b2 = 1.442710462, b1 = -2.118377725, b0 = 1;
  double c3 = 1.641345311, c2 = 3.429567803, c1 = -1.62490649, c0 = -1.97084045;
  double d2 = 1.637067800, d1 = 3.543889200, d0 = 1;
  double x2, r, y;
  int sign_x;

  if (x < -1 || x > 1){
    cout << "ERROR: erfinv(x) function called with fail x, x must be >=-1 AND <=1" << endl;
    return -1;
  }

  if (x == 0)
    return 0;

  if (x > 0)
    sign_x = 1;
  else {
    sign_x = -1;
    x = -x;
  }
 
  if (x <= 0.7) {
    x2 = x * x;
    r = x * (((a3 * x2 + a2) * x2 + a1) * x2 + a0);
    r /= (((b4 * x2 + b3) * x2 + b2) * x2 + b1) * x2 + b0;
  }
  else {
    y = sqrt(-log ((1 - x) / 2));
    r = (((c3 * y + c2) * y + c1) * y + c0);
    r /= ((d2 * y + d1) * y + d0);
  }
 
  r = r*sign_x;
  x = x*sign_x;
 
  r -= (erf(r) - x) / (2 / sqrt(M_PI) * exp(-r * r));
  r -= (erf(r) - x) / (2 / sqrt(M_PI) * exp(-r * r));
 
  return r;
}

//--------------------------------------------------------------
double correlCoefficient(const vector<double> &x, const vector<double> &y)
{
  if(x.size() != y.size()){
    cout << "\n\n !!!!! ERROR !!!!!\n Size of x and y are not matching in CorrelCoef, Statistic.cpp\n";
    return 0;
  }
  int n = x.size();
  double szam=0.0,sx=0.0,sy=0.0,ax=0.0,ay=0.0;
  for(int i=0;i<n;i++){
    ax += x.at(i);
    ay += y.at(i);
  }
  ax /= n; ay /= n;
  for(int i=0;i<n;i++)
    szam += (x.at(i)-ax)*(y.at(i)-ay);
  for(int i=0;i<n;i++){
    sx += pow(x.at(i)-ax,2.);
    sy += pow(y.at(i)-ay,2.);
  }
  sx = pow(sx,0.5);
  sy = pow(sy,0.5);
  return szam/sx/sy;
}

//--------------------------------------------------------------
double spearmanCoefficient(const vector<double> &xo, const vector<double> &yo)
{
  double sc;
  if(xo.size() == yo.size())
  {
    int nv = xo.size();
    vector<double> xs = xo, ys = yo;
    selectionSort(xs);
    vector<int> rx(xo.size()),ry(yo.size());
    for(int i=0; i<xs.size(); i++)
    {
      for(int j=0; j<xo.size(); j++)
      {
        if(xo[i] == xs[j])
        {
          rx[i] = j;
          break;
        }
      }
    }

    selectionSort(ys);
    for(int i=0; i<ys.size(); i++)
    {
      for(int j=0; j<yo.size(); j++)
      {
        if(yo[i] == ys[j])
        {
          ry[i] = j;
          break;
        }
      }
    }

    double d=0.;
    for(int i=0; i<rx.size(); i++)
      d += pow(rx[i]-ry[i],2);
    
    sc = 1.-6.*d/nv/(nv*nv-1.);
  }
  else
  {
    cout << endl << "!WARNING! Something is wrong. xo size " << xo.size() << " does not match with yo size " << yo.size() << " in spearmanCoefficient function. " << endl;
    exit(-1);
  }

  return sc;
}

//--------------------------------------------------------------
vector<double> grubbsTest(const vector<double> &x)
{
  double a,m,s;
  a = average(x);
  m = absoluteMax(x);
  s = standardDeviation(x);
  vector<double> y;
  for(int i=0; i<x.size(); i++)
    if(x.at(i) > m-3.*s && x.at(i) < m+3.*s)
      y.push_back(x.at(i));
  return y;
}

//--------------------------------------------------------------
double linearInterpolation(const vector<double> &x, const vector<double> &y, double xp)
{
  double yp = 0.0;
  int nx=x.size(), ny=y.size();
  if(nx == 0 || ny == 0 || nx == 1 || ny == 1)
  {
    cout << endl << " !!!ERROR!!! interpolate: Number of x/y coordinates (" << nx << "/" << ny << " is equal to 0 or 1!" << endl << " Exiting..." << endl;
    exit(0);
  }
  else if(nx!=ny)
  {
    cout << endl << " !!!ERROR!!! Number of x coordinates (" << nx << ") NOT equal to number of y coordinates(" << ny << ") !" << endl;
    exit(-1);
  }
  else
  {
    int idx;
    if(xp<min(x,idx)){
      //cout << endl << " !!!WARNING!! interpolate: xp = " << xp << " is smaller than min(x) = " << min(x,idx) << endl << " Overriding xp = min(x), then Continouing..." << endl;
      xp = min(x,idx);
      yp = y[idx];
    }
    else if(xp>max(x,idx))
    {
      //cout << endl << " !!!WARNING!! xp = " << xp << " is larger than max(x) = " << max(x,idx) << endl << " Overriding xp = max(x), then Continouing..." << endl;
      xp = max(x,idx);
      yp = y[idx];
    }
    else
    {
      // Finding the closest points
      int idx;
      double x1=min(x,idx);
      double y1=y[idx];
      double x2=max(x,idx);
      double y2=y[idx];
      for(int i=0; i<nx; i++)
      {
        if(x[i]<xp && x[i]>x1){
          x1 = x[i];
          y1 = y[i];
        }
        if(x[i]>xp && x[i]<x2){
          x2 = x[i];
          y2 = y[i];
        }
      }
      double a = (y2-y1)/(x2-x1);
      double b = y2 - a*x2;
      yp = a*xp + b;
    }
  }

  return yp;
}

//--------------------------------------------------------------
double interpolate(vector<double> x, vector<double> y, double xp){

  double yp = 0.0;
  int nx=x.size(), ny=y.size();
  if(nx == 0 || ny == 0 || nx == 1 || ny ==1)
  {
    cout << endl << " !!!ERROR!!! interpolate: Number of x/y coordinates (" << nx << "/" << ny << " is equal to 0 or 1!" << endl << " Exiting..." << endl;
    exit(0);
  }
  if(nx == ny)
  {
    int idx;
    if(xp<min(x,idx)){
      cout << endl << " !!!WARNING!! interpolate: xp = " << xp << " is smaller than min(x) = " << min(x,idx) << endl << " Overriding xp = min(x), then Continouing..." << endl;
      xp = min(x,idx);
      yp = y[idx];
    }
    else if(xp>max(x,idx))
    {
      cout << endl << " !!!WARNING!! xp = " << xp << " is larger than max(x) = " << max(x,idx) << endl << " Overriding xp = max(x), then Continouing..." << endl;
      xp = max(x,idx);
      yp = y[idx];
    }else
    {
      if(nx == 2) // Fitting a linear
      {
        double a = (y[1]-y[0])/(x[1]-x[0]);
        double b = y[1] - a*x[1];
        yp = a*xp + b;
      }
      if(nx>=3)// Above 3 points, fitting a parabole with the three closest
      {
        vector<double> xdist(nx);
        vector<int> idx;
        for(int i=0; i<nx; i++)
          xdist[i] = abs(x[i] - xp);
        quickSortWithIndicies(xdist,idx,0,nx-1);
        yp = (xp-x[idx[1]])*(xp-x[idx[2]])/(x[idx[0]]-x[idx[1]])/(x[idx[0]]-x[idx[2]])*y[idx[0]] + (xp-x[idx[0]])*(xp-x[idx[2]])/(x[idx[1]]-x[idx[0]])/(x[idx[1]]-x[idx[2]])*y[idx[1]] + (xp-x[idx[0]])*(xp-x[idx[1]])/(x[idx[2]]-x[idx[0]])/(x[idx[2]]-x[idx[1]])*y[idx[2]]; // Lagrange polinomial
      }
    }
  }
  else
  {
    cout << endl << " !!!ERROR!!! Number of x coordinates (" << nx << ") NOT equal to number of y coordinates(" << ny << ") !" << endl;
    exit(0);
  }

  return yp;
}

//-------------------------------------------------------------- 
vector<int> unique(const vector<int> &x)
{
  vector<int> out;
  for(int i=0; i<x.size(); i++)
  {
    bool unique = true;
    for(int j=0; j<out.size(); j++)
    {
      if(x[i] == out[j])
      {
        unique = false;
        break;
      }
    }
    if(unique)
    {
      out.push_back(x[i]);
    }
  }

  return out;
}

//-------------------------------------------------------------- 
void selectionSort(vector<double> &x)  
{  
  int i, j, min_idx, n = x.size();

  // One by one move boundary of unsorted subarray  
  for(i = 0; i < n-1; i++)  
  {  
    // Find the minimum element in unsorted array  
    min_idx = i;  
    for(j = i+1; j < n; j++)  
    if(x[j] < x[min_idx])  
        min_idx = j;  

    // Swap the found minimum element with the first element  
    swap(x, min_idx, i);  
  }  
}  

//-------------------------------------------------------------- 
void quickSortWithIndicies(vector<double> &x, vector<int> &idx, int low, int high){

  idx.clear();
  idx.resize(x.size());
  for(int i=0; i<idx.size(); i++)
    idx[i] = i;
  quickSort(x,idx,low,high);
}
//-------------------------------------------------------------- 
void quickSortWithIndicies(vector<int> &x, vector<int> &idx, int low, int high){

  idx.clear();
  idx.resize(x.size());
  for(int i=0; i<idx.size(); i++)
    idx[i] = i;
  quickSort(x,idx,low,high);
}

/* The main function that implements QuickSort 
 arr[] --> Array to be sorted, 
  low  --> Starting index, 
  high  --> Ending index */
void quickSort(vector<double> &x, vector<int> &idx, int low, int high) 
{
  if (low < high) 
  { 
    /* pi is partitioning index, arr[p] is now 
       at right place */
    int pi = partition(x, idx, low, high); 

    // Separately sort elements before 
    // partition and after partition 
    quickSort(x, idx, low, pi - 1); 
    quickSort(x, idx, pi + 1, high); 
  } 
}
void quickSort(vector<int> &x, vector<int> &idx, int low, int high) 
{
  if (low < high) 
  { 
    /* pi is partitioning index, arr[p] is now 
       at right place */
    int pi = partition(x, idx, low, high); 

    // Separately sort elements before 
    // partition and after partition 
    quickSort(x, idx, low, pi - 1); 
    quickSort(x, idx, pi + 1, high); 
  } 
} 

/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition(vector<double> &x, vector<int> &idx, int low, int high) 
{ 
  int pivot = x[high];    // pivot 
  int i = (low - 1);  // Index of smaller element 

  for (int j = low; j <= high- 1; j++) 
  { 
    // If current element is smaller than or 
    // equal to pivot 
    if (x[j] <= pivot) 
    { 
      i++;    // increment index of smaller element 
      swap(x,i,j);
      if(idx.size() == x.size())
        swap(idx,i,j);
    } 
  } 
  swap(x,i+1,high); 
  if(idx.size() == x.size())
    swap(idx,i+1,high);
  return (i + 1); 
}
int partition(vector<int> &x, vector<int> &idx, int low, int high) 
{ 
  int pivot = x[high];    // pivot 
  int i = (low - 1);  // Index of smaller element 

  for (int j = low; j <= high- 1; j++) 
  { 
    // If current element is smaller than or 
    // equal to pivot 
    if (x[j] <= pivot) 
    { 
      i++;    // increment index of smaller element 
      swap(x,i,j);
      if(idx.size() == x.size())
        swap(idx,i,j);
    } 
  } 
  swap(x,i+1,high); 
  if(idx.size() == x.size())
    swap(idx,i+1,high);
  return (i + 1); 
}

//--------------------------------------------------------------
void swap(vector<double> &x, int i, int j)
{
  double temp = x[i];
  x[i] = x[j];
  x[j] = temp;
}

//--------------------------------------------------------------
void swap(vector<int> &x, int i, int j)
{
  int temp = x[i];
  x[i] = x[j];
  x[j] = temp;
}

//--------------------------------------------------------------
double min(vector<double> x, int &idx)
{
  double min=x[0];
  idx = 0;
  for(int i=1; i<x.size(); i++)
  {
    if(min > x[i])
    {
      min = x[i];
      idx = i;
    }
  }
  return min;
}
//--------------------------------------------------------------
double min(vector<double> x)
{
  int idx;
  return min(x,idx);
}

//--------------------------------------------------------------
double max(vector<double> x, int &idx)
{
  double max=x[0];
  idx = 0;
  for(int i=1; i<x.size(); i++)
  {
    if(max < x[i])
    {
      max = x[i];
      idx = i;
    }
  }
  return max;
}
//--------------------------------------------------------------
int max(vector<int> x, int &idx)
{
  int max=x[0];
  idx = 0;
  for(int i=1; i<x.size(); i++)
  {
    if(max < x[i])
    {
      max = x[i];
      idx = i;
    }
  }
  return max;
}
//--------------------------------------------------------------
double max(vector<double> x)
{
  int idx;
  return max(x,idx);
}
//--------------------------------------------------------------
int max(vector<int> x)
{
  int idx;
  return max(x,idx);
}


//--------------------------------------------------------------
VectorXd leastSquaresPolynomial(const VectorXd &x, const VectorXd &y, int &order)
{

  if(x.rows() != y.rows())
  {
    cout << endl << "!!!ERROR!!! Statistic::leastSquaresPolynomial(): size of X is not equal to size of Y" << endl << "Exiting..." << endl;
    exit(0);
  }
  if(order>x.rows()-1)
  {
    cout << endl << "!WARNING! Statistic::leastSquaresPolynomial(): order(" << order << ") is larger or equal than number of points(" << x.rows() << ")!" << endl << "order is decreased to x.rows()-1 (" << x.rows()-1 << "), then continouing..." << endl;
    order = x.rows()-1;
  }

  VectorXd out = VectorXd::Zero(order+1);
  MatrixXd Gram = MatrixXd::Zero(order+1,order+1);

  for(int i=0; i<Gram.rows(); i++)
  {
    for(int j=i; j<Gram.cols(); j++)
    {
      Gram(i,j) = (eigenVectorXdPow(x,i).transpose()*eigenVectorXdPow(x,j)).sum(); // <fi,fj> = [x1^i, ... xm^i]*[x1^j, ... xm^j]^T
      if(j!=i)
        Gram(j,i) = Gram(i,j); // Since Gram matrix is symmetric
    }
  }

  VectorXd b = VectorXd::Zero(order+1);
  for(int i=0; i<b.rows(); i++)
    b(i) = (y.transpose()*eigenVectorXdPow(x,i)).sum();

  out = Gram.ldlt().solve(b); //LDL decomposition is suitable, since Gram is positive definite (also symmetric)
  //out = Gram.fullPivLu().solve(b);

  return out;
}

//--------------------------------------------------------------
VectorXd eigenVectorXdPow(const VectorXd &x, int order)
{
  int n=x.rows();
  VectorXd out;
  if(order==0)
    out = VectorXd::Constant(n,1.);
  else
    out = x;

  if(order>=2 || order <=-2)
  {
    for(int i=0; i<n; i++)
      for(int j=0; j<order-1; j++)
        out(i) = out(i)*out(i);
  }

  if(order<0)
  {
    for(int i=0; i<n; i++)
      out(i) = 1./out(i);
  }
  return out;
}

//--------------------------------------------------------------
vector<double> findMaxValues(const VectorXd &f, vector<int> &vidx, int n, bool av)
{
  vector<double> out;
  vidx.clear();
  for(int j=0; j<n; j++)
  {
    double max = -1e20;
    int idx = 0;
    for(int i=0; i<f.size(); i++)
    {
      bool gotit = false;
      for(int k=0; k<vidx.size(); k++)
      {
        if(i==vidx[k])
          gotit = true;
      }
      if(av) // with absolute value
      {
        if(abs(f(i))>max && !gotit)
        {
          max = abs(f(i));
          idx = i;
        }
      }
      else // without abosolute value
      {
        if(f(i)>max && !gotit)
        {
          max = f(i);
          idx = i;
        }
      }
    }
    out.push_back(max);
    vidx.push_back(idx);
  }
  return out;
}