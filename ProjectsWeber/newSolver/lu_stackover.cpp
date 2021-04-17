#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>

 #include "LUDecomp.h"

 using namespace std;

 LUDecomp::LUDecomp()
 {
 }

 void LUDecomp::h_pivot_decomp(int MAT1, double a[], int p[], int q[])
 {
int i = 0, j = 0, k = 0;
int n = MAT1;
int pi = 0, pj = 0, tmp = 0;
double max = 0.0;
double ftmp = 0.0;
//Stores the scaling of each row or column.
double* vv = new double[MAT1 + 5];

//Loop over rows toget the implicit scaling information.
max = 0.0;
for (i = 0; i < n; i++)
{
    for (j = 0; j < n; j++)
    {

        if ((ftmp = fabs(a(i,j))) > max)
        {
            max=ftmp;
        }
    }
    //No nonzero largest element.
                if (max == 0.0)
                {
                    throw("Singular matrix in LUdcmp");
                }
                //Save the scaling.
                vv[i]=1.0/max;

            }

        // The k element determines which pivot element you are in thereby
        // determining the submatrix starting at the upper left corner of the matrix.
for (k = 0; k < n; k++)
{

    // pi: stores row needing to be swapped.
    // pj: stores column needing to be swapped.
    // max: makes a zero element in the matrix into a very tiny number.
    pi = -1, pj = -1, max = TINY;

    //find pivot in submatrix a(k:n,k:n) by finding the absolute value of the biggest element.
    for (i = k; i < n; i++)
    {
        for (j = k; j < n; j++)
        {
            //j = k;
            ftmp = vv[i] * fabs(a(i,j));
            // Decides if current max is bigger than current element.
                    if (ftmp>max)
                    {
                        max = ftmp;
                        // Index of row being swapped.
                        pi=i;
                        // Index of column being swapped.
                        pj=j;
                    }
                }
            }

    {
        // Stores the permutation of row swaps.
        tmp = p[k];
        p[k] = p[pi];
        p[pi] = tmp;
    }

    //Swaps the scalling factor if needed.
    if (k != pi)
    {
        vv[pi] = vv[k];
        cout << "Scaling factor: " << vv[pi] << endl;
    }

    // Swaps the indicated rows to move the max pivot
    // element of the submatrix k into place.
    for (j = 0; j < n; j++)
    {
        // The k and pi index stays the same so the row
        // number stays the same, the j changes to iterate threw the row.
        ftmp = a(k,j);
        a(k,j)=a(pi,j);
        a(pi,j)=ftmp;
        //cout << a(k,j) << " , " << a(pi,j) << endl;
            }

    {
        // Stores the permutation of column swaps.
        tmp = q[k];
        q[k] = q[pj];
        q[pj] = tmp;
        //cout << q[k] << " , " << q[pj] << endl;
    }

    // Swaps the indicated columns to move the max pivot
    // element of the submatrix k into place.
    for (i = 0; i < n; i++)
    {
        // The k and pj index stays the same so the column
        // number stays the same, the i changes to iterate threw the column.
        ftmp = a(i,k);
        a(i,k)=a(i,pj);
        a(i,pj)=ftmp;
        //cout << a(i,k) << " , " << a(i,pj) << endl;
            }
        // END PIVOT
    cout << fixed << showpoint;
    cout << setprecision(20);
    // Check pivot size and decompose
    if ((fabs(a(k,k))>TINY))
    {
        for (i=k+1;i<n;i++)
        {
            // Column normalisation, Does first element under pivot k row i.

            ftmp=a(i,k)/=a(k,k);

            cout << "k,k " <<a(k,k) << " , " << endl;
            // Does the rest of row i.
            for (j=k+1;j<n;j++)
            {
                //a(ik)*a(kj) subtracted from lower right submatrix elements
                a(i,j)-=(ftmp*a(k,j));
                //cout <<"i,j "<< a(i,j) << endl;
            }
        }
    }

}
    //END DECOMPOSE
for (i = 0; i < n; i++)
{
    for (j = 0; j < n; j++)
    {

        cout << a(i,j)<<" ";
    }
    cout << endl;
}
 }

 void LUDecomp::h_solve(int MAT1, double a[], double x[], int p[], int q[])
 {
// Forward substitution; see  Golub, Van Loan 96
// And see http://www.cs.rutgers.edu/~richter/cs510/completePivoting.pdf
int i = 0, ii = 0, j = 0;
double ftmp = 0.0;
double* xtmp = new double[MAT1 + 5];

cout << fixed << showpoint;
cout << setprecision(4);

// Swap rows
// Put be vector back like it should be by using the permutations from the row swapping.
for (i = 0; i < MAT1; i++)
{
    xtmp[i] = x[p[i]]; //value that should be here
    //cout << xtmp[i] << endl;
}

// Ly=b
for (i = 0; i < MAT1; i++)
{
    ftmp = xtmp[i];
    if (ii != 0)
        for (j = ii - 1; j < i; j++)
            ftmp -= a(i,j)*xtmp[j];

            else if (ftmp!=0.0)
            ii=i+1;

    xtmp[i] = ftmp;
    //cout << xtmp[i] << endl;
}

// Backward substitution
// Partially taken from Sourcebook on Parallel Computing p577
// Solves Ux=y
cout << "xtmp " << xtmp[MAT1 - 1] << " a " << a(MAT1-1,MAT1-1)<< endl;
xtmp[MAT1 - 1] /= a(MAT1-1,MAT1-1);
//cout << xtmp[MAT1 - 1] << endl;
for (i = MAT1 - 2; i >= 0; i--)
{
    ftmp = xtmp[i];
    //cout << "ftmp " << ftmp << endl;
    for (j = i + 1; j < MAT1; j++)
    {
        ftmp -= a(i,j)*xtmp[j];
        //cout << "ftmp in "<<ftmp << endl;
    }

    xtmp[i] = (ftmp) / a(i,i);

}

    // Last bit
    // Swap columns
    // Takes the final answer and puts it back into its proper order by
    // using the permutations from the column swapping.
for (i = 0; i < MAT1; i++)
{
    x[q[i]] = xtmp[i];
}

delete xtmp;
 }

 // Method to get output from the LU Decomposition.
 void LUDecomp::output(unsigned int MAT1, double a[], double b[])
 {

// Pivot array's for the permutation vectors.
int* p_pivot = new int[MAT1 + 5];
int* q_pivot = new int[MAT1 + 5];

// Sets the elements in the permutation vectors up to receive permutations.
// p_pivot is for row permutations and is initialized to {0,1,...,r};
// q_pivot is for column permutations and is initialized to {0,1,...,r};
for (unsigned int i = 0; i < MAT1; i++)
{
    p_pivot[i] = i;
    q_pivot[i] = i;
}

// Call to decomposition method passing (size,matrix to be decomposed, not used,   not used).
h_pivot_decomp(MAT1, a, p_pivot, q_pivot);

// Call to solve passing (size, matrix in LU form, b vector, not used, not used).
h_solve(MAT1, a, b, p_pivot, q_pivot);

// Have solution.
// Used for file output.
ofstream outFile;
// Allow for appenending to a file already created.
outFile.open("outSolMatrix.txt");

// Sets the precision of the output to the file.
outFile << fixed << showpoint;
outFile << setprecision(4);

// Output results to file answer is {0,1,...,n}.
for (unsigned int i = 0; i < MAT1; i++)
{
    outFile << i << " " << b[i] << endl;
}

outFile << "End" << endl;

delete p_pivot;
delete q_pivot;

outFile.close();
 }