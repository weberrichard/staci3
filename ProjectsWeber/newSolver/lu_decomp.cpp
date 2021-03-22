#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

class jacMatrix
{
public:
	// elements along the diag
	vector<double> d;

	// row and col index, COL major
	vector<int> ic;
	vector<int> jc;
	vector<double> xc;

	int len=0;

	jacMatrix(int l)
	{
		len = l;
	}
	void insertDiag(double x)
	{
		d.push_back(x);
	}
	void insertNonDiag(int i, int j, double x)
	{
		ic.push_back(i);
		jc.push_back(j);
		xc.push_back(x);
	}
	void updateDiag(int i, double x)
	{
		d[i] = x;
	}
	void print()
	{
		cout << "JAC: " << endl;
		printf("diag:\n");
		for(int k=0; k<d.size(); k++)
		{
			printf("[%2i] : %6.4f\n",k,d[k]);
		}
		printf("------------------------\n");
		printf("[%2s] %2s-%-2s: %6s\n","k","ic","jc","xc");
		for(int k=0; k<ic.size(); k++)
		{
			printf("[%2i] %2i-%-2i: %6.4f\n",k,ic[k],jc[k],xc[k]);
		}

	}
};

class sparseL
{
public:
	// bottom left
	vector<int> i; // row index
	vector<int> j; // col index
	vector<double> x; // actual value
	int n=0; // number of non-zero at bottom left part

	// bottom right
	vector<int> i2; // row index
	vector<int> j2; // col index
	vector<double> x2; // actual value
	sparseL(){};
	void insert(int ii, int jj, double xx)
	{
		i.push_back(ii);
		j.push_back(jj);
		x.push_back(xx);
		n++;
	}
	int insert2(int ii, int jj, double xx)
	{
		int out;
		int k=0;
		if(k==i2.size())
		{
			i2.push_back(ii);
			j2.push_back(jj);
			x2.push_back(xx);
			out = -1;
		}
		else
		{
			bool gotIt = false;
			while(!gotIt && k<i2.size()-1)
			{	
				bool b1 = (ii==i2[k] && jj>j2[k]);
				bool b2 = (i2[k]==i2[k+1] && jj<j2[k+1]);
				bool b3 = (i2[k]<i2[k+1]);
				if(b1 && (b2|b3))
				{
					i2.insert(i2.begin()+k+1,ii);
					j2.insert(j2.begin()+k+1,jj);
					x2.insert(x2.begin()+k+1,xx);
					gotIt = true;
					out = k+1;
				}
				else
				{
					k++;
				}
			}
			if(!gotIt)
			{
				i2.push_back(ii);
				j2.push_back(jj);
				x2.push_back(xx);
				out = -1;
			}
		}
		return out;
	}
	void update(int k, double xx)
	{
		x[k] = xx;
	}
	void update2(int k, double xx)
	{
		x2[k] = xx;
	}
	void clear()
	{
		i.clear();
		j.clear();
		x.clear();
		i2.clear();
		j2.clear();
		x2.clear();
	}
	void print()
	{
		cout << "L: " << endl;
		printf("[%2s] %2s-%-2s: %6s\n","k","i","j","x");
		for(int k=0; k<i.size(); k++)
		{
			printf("[%2i] %2i-%-2i: %6.4f\n",k,i[k],j[k],x[k]);
		}
		printf("------------------------\n");
		for(int k=0; k<i2.size(); k++)
		{
			printf("[%2i] %2i-%-2i: %6.4f\n",k,i2[k],j2[k],x2[k]);
		}
	}
};

class sparseU
{
public:
	vector<double> d; // for diag elements

	// for non-diag elements at upper right
	vector<int> i; 
	vector<int> j;
	vector<double> x;

	// for non-diag elements at upper right
	vector<int> i2; 
	vector<int> j2; 
	vector<double> x2;

	sparseU(){};
	void insertDiag(double xx)
	{
		d.push_back(xx);
	}
	void insert(int ii, int jj, double xx)
	{
		i.push_back(ii);
		j.push_back(jj);
		x.push_back(xx);
	}
	void insert2(int kk, int ii, int jj, double xx)
	{
		if(kk==-1)
		{
			i2.push_back(ii);
			j2.push_back(jj);
			x2.push_back(xx);
		}
		else
		{
			i2.insert(i2.begin()+kk,ii);
			j2.insert(j2.begin()+kk,jj);
			x2.insert(x2.begin()+kk,xx);
		}
	}
	void updateDiag(int k, double xx)
	{
		d[k] = xx;
	}
	void update(int k, double xx)
	{
		x[k] = xx;
	}
	void clear()
	{
		d.clear();
		i.clear(); i2.clear();
		j.clear(); j2.clear();
		x.clear(); x2.clear();
	}
	void print()
	{
		printf("U:\n");
		printf("diag:");
		for(int k=0; k<d.size(); k++)
		{
			printf("%6.4f, ", d[k]);
		}
		printf("\n");
		printf("[%2s] %2s-%-2s: %6s\n","k","i","j","x");
		for(int k=0; k<i.size(); k++)
		{
			printf("[%2i] %2i-%-2i: %6.4f\n",k,i[k],j[k],x[k]);
		}
		printf("------------------------\n");
		for(int k=0; k<i2.size(); k++)
		{
			printf("[%2i] %2i-%-2i: %6.4f\n",k,i2[k],j2[k],x2[k]);
		}
	}
};

void LUDecomposition(const jacMatrix jac, sparseL &L, sparseU &U, bool isFirst)
{
	// clearing output vars if it is first
	if(isFirst)
	{
		L.clear();
		U.clear();
	}

	// calculate upper and left parts
	int n_u = jac.ic[0]; // size of upper part
	int n_jac = jac.len; // size of jac matrix
	// diag of U
	if(isFirst)
	{
		for(int i=0; i<n_u; i++)
		{
			U.insertDiag(jac.d[i]);
		}
		// lower left of L
		for(int i=0; i<jac.ic.size(); i++)
		{
			L.insert(jac.ic[i],jac.jc[i],jac.xc[i]/jac.d[jac.jc[i]]);
		}
		// upper right of U
		for(int i=0; i<jac.ic.size(); i++)
		{
			U.insert(jac.jc[i],jac.ic[i],jac.xc[i]);
		}
	}
	else
	{
		for(int i=0; i<n_u; i++)
		{
			U.updateDiag(i,jac.d[i]);
		}
		// lower left of L
		for(int i=0; i<jac.ic.size(); i++)
		{
			L.update(i,1./jac.d[i]);
		}
	}

	// calculate bottom right corner
	if(isFirst)
	{
		int k_d=0, k_d2=0;
		int k_l=0, k_l2=0;
		int k_l0=0, k_l20=0;
		int k_u=0, k_u2=0;
		for(int i=n_u; i<n_jac; i++)
		{
			// diag of U
			double uii=0.;
			k_l0 = k_d;
			while(i==L.i[k_d])
			{
				uii -= L.x[k_d]*U.x[k_d];
				k_d++;
			}
			if(L.i2.size()>0)
			{	
				while(i>L.i2[k_d2])
				{
					k_d2++;
				}
				k_l2=k_d2;
				while(i==L.i2[k_d2])
				{
					uii -= L.x2[k_d2]*U.x2[k_d2];
					k_d2++;
				}
			}
			uii += jac.d[i];
			U.insertDiag(uii);

			// row of U
			k_u = k_d;
			k_u2 = k_d2;
			for(int j=i+1; j<jac.len; j++)
			{	
				double uij=0.;
				// product from upper part of U
				k_l = k_l0;
				while(j>U.j[k_u])
				{
					k_u++;
				}
				while(i==L.i[k_l] && j==U.j[k_u])
				{
					if(L.j[k_l] == U.i[k_u])
					{
						uij -= L.x[k_l]*U.x[k_u];
						k_u++;
						k_l++;
					}
					else if(L.j[k_l]>U.i[k_u])
					{
						k_u++;
					}
					else
					{
						k_l++;
					}
				}
				// product from lower part of U
				k_l2 = k_l20;
				if(L.i2.size()>0)
				{
					while(j>U.j2[k_u2])
					{
						k_u2++;
					}
					while(i==L.i2[k_l2] && j==U.j2[k_u2])
					{
						if(L.j2[k_l2] == U.i2[k_u2])
						{
							uij -= L.x2[k_l2]*U.x2[k_u2];
							k_u2++;
							k_l2++;
						}
						else if(L.j2[k_l2]>U.i2[k_u2])
						{
							k_u2++;
						}
						else
						{
							k_l2++;
						}
					}
				}
				// todo rest of L*U
				if(uij != 0)
				{
					auto started = std::chrono::high_resolution_clock::now();
					int kk = L.insert2(j,i,uij/uii);
					auto done = std::chrono::high_resolution_clock::now();
					std::cout << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << endl;
					U.insert2(kk,i,j,uij);
				}
			}
		}
	}
	else
	{
		cout << "LOFASZ" << endl;
	}
	printf("FERTICH!!!\n");
}

void loadJac(string filename, jacMatrix &jac)
{
  ifstream ifile(filename);
  string temp, line;
  if(ifile.is_open())
  {
		while(getline(ifile,line))
		{
			int i,j;
			double x;
			int counter=0;
			for (string::iterator s=line.begin(); s!=line.end(); s++)
			{
				if(*s!=',')
				{
					temp +=*s;
				}
				else
				{
					if(counter==0)
					{
						i=stoi(temp);
						counter++;
					}
					else if(counter==1)
					{
						j=stoi(temp);
						counter++;
					}
					temp = "";
				}
			}
			x=stod(temp,0);
			temp="";
			if(i==j)
			{
				jac.insertDiag(x);
			}
			else if(i>j)
			{
				jac.insertNonDiag(i-1,j-1,x);
			}
		}
  }
  ifile.close();
}

int main()
{
	/*int n=9;
	jacMatrix jac(n);
	jac.insertDiag(5.);
	jac.insertDiag(6.);
	jac.insertDiag(7.);
	jac.insertDiag(8.);
	jac.insertDiag(9.);
	jac.insertDiag(14.);
	jac.insertDiag(13.);
	jac.insertDiag(12.);
	jac.insertDiag(11.);
	jac.insertNonDiag(5,1,1.);
	jac.insertNonDiag(5,4,-1.);
	jac.insertNonDiag(6,2,-1.);
	jac.insertNonDiag(6,3,1.);
	jac.insertNonDiag(7,2,1.);
	jac.insertNonDiag(7,4,1.);
	jac.insertNonDiag(8,0,-1.);
	jac.insertNonDiag(8,1,-1.);
	jac.insertNonDiag(8,3,-1.);*/

	/*int n=9;
	jacMatrix jac(n);
	loadJac("jacobi2_grid.txt",jac);
	cout << "LOAD OK " << endl;*/

	/*int n=22;
	jacMatrix jac(n);
	loadJac("jacobi2_ohermes.txt",jac);
	cout << "LOAD OK " << endl;*/

	//jac.print();

	int n=14543;
	jacMatrix jac(n);
	loadJac("jacobi2_ferto.txt",jac);
	cout << "LOAD OK " << endl;

	sparseL L;
	sparseU U;	

 	srand((unsigned int) time(0));
 	clock_t ido = clock();
	LUDecomposition(jac,L,U,true);
	cout << endl << "\nDecomp:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

	cout << "DECOMP OK " << endl;

	cout << "L nonzero: " << L.i.size() << "  " << L.i2.size() << endl;
	cout << "U nonzero: " << U.i.size() << "  " << U.i2.size() << endl;

	//L.print();

	//jac.updateDiag(0,4.);
	//jac.updateDiag(1,3.);
	//jac.updateDiag(2,2.);
	//jac.updateDiag(3,5.);
	//jac.updateDiag(4,6.);
	//LUDecomposition(jac,L,U,false);

	cout << endl << endl;
	return 0;
}
