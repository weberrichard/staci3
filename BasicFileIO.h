/*===================================================================*\
                              BasicFileIO
                            ---------------

  Reading and writing data from and to files.

  Can be used independently from staci3.
  
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef BASICFILEIO_H
#define BASICFILEIO_H

#include <vector>
#include <iostream>
#include <fstream>
#include </usr/include/eigen3/Eigen/Eigen>
#include <sys/stat.h> // mkdir

using namespace std;
using namespace Eigen;

// Reading doubles from file, each line will be an element
vector<string> readVectorString(string fileName);

// Reading doubles from file to Eigen MatrixXd, separeted with "separator"
MatrixXd readMatrixXdDouble(string fileName, char separator);

// Reading doubles from file to vector<double>
vector<double> readVectorDouble(string fileName);
vector<int> readVectorInt(string fileName);

// Counting the rows for Eigen Vectors
int countRows(string fileName);

// Counting the rows and cols for Eigen Vectors
vector<int> countRowsCols(string fileName, char separator);

// Writing vector double to file
void writeVectorDouble(string filename, vector<double> v);

// Make new directory, works for windows and linux
void makeDirectory(string name);

#endif
