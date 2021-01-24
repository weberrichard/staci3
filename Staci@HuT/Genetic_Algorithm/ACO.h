#include "Randoms.cpp"

class ACO {
public:
	ACO (int nAnts, int nCities, 
		double alpha, double beta, double q, double ro, double taumax,
		int initCity);
	virtual ~ACO ();
	
	void init ();
	
	void connectNODES (int node_in, int node_out);
	void setCITYPOSITION (int city, double x, double y);
	void Search_Graph_Build_Up();
	double Fitness_Function_AVG_Sens (int node_in, int node_out);
	double Fitness_Function_Peak_Sens (int node_in, int node_out);
	void printPHEROMONES ();
	void printGRAPH ();
	void printRESULTS ();
	
	void optimize (int ITERATIONS);

private:
	double distance (int node_in, int node_out);
	bool exists (int node_in, int node_out);
	bool vizited (int antk, int c);
	double PHI (int node_in, int node_out, int antk);
	
	double fitness_Sum (int antk);
	
	int NODE ();
	void route (int antk);
	int valid (int antk, int iteration);
	
	void updatePHEROMONES ();

	
	int NUMBEROFANTS, NUMBEROFNODES, INITIALNODE;
	double ALPHA, BETA, Q, RO, TAUMAX;
	
	double BESTLENGTH;
	int *BESTROUTE;

	int **GRAPH, **ROUTES;
	double **CITIES, **PHEROMONES, **DELTAPHEROMONES, **PROBS;

	Randoms *randoms;
};

