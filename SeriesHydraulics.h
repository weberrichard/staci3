/*===================================================================*\
                            SeriesHydraulics
                            ----------------
 
  staci3 is using Eigen, see http://eigen.tuxfamily.org

    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef SERIESHYDRAULICS_H
#define SERIESHYDRAULICS_H

#include "Sensitivity.h"
#include "Rule.h"

class SeriesHydraulics : public Sensitivity
{
public:

	SeriesHydraulics(string fileName);
	virtual ~SeriesHydraulics();

  // control parameteres of the series sensitivity
  struct ssc
  {
    string parameter; // roughness, diameter or demand
    vector<double> time; // time appointments where sensitivity calculation must be done
    vector<int> rowIndex; // row index of stored sensitivity matrix, all if empty
    vector<int> colIndex; // col index of stored sensitivity matrix, all if empty
    bool isPressure; // true if pressureSensitivity, false if massFlowRateSensitivity
  };

  // collecting the time stamps where calculations were made
  vector<double> vectorTime;

  // calculating the hydraulics in time
  void seriesSolve(ssc seriesSensitivityControl);
  void seriesSolve();
  // deciding between solveSystem and calculateSensitivity basd on series sensitivity control parameters
  bool isSensitivityCalculation(ssc seriesSensitivityControl);
  // saving the important part of the sensitivity matrix in seriesSensitivity
 // void saveSensitivityMatrix(ssc seriesSensitivityControl);

  // storing the sensitivity matrices in memory, only important part
  vector<MatrixXd> seriesSensitivity;

  // Printing results to console
  void timeTableNode(vector<string> ID, string unit);
  void timeTableNode(vector<int> idx, string unit);
  void timeTableEdge(vector<string> ID, string unit);
  void timeTableEdge(vector<int> idx, string unit);
  // saving heads, consumtions for nodes and flow rate for edges in file
  void saveTimeResult(vector<string> edgeID, vector<string> nodeID, string qUnit, string hUnit);

  void addToPipeCVIndex(int idx);
  void deleteFromPipeCVIndex(int idx);

protected:
  // updating stuffs
  void updateDemand();
  void updatePressurePointPattern();
  void updatePool();
  void updateControl();
  void updateRule();

  // saving heads, consumtions for nodes and flow rate for edges in vector* variables
  void saveOutput();
  // clearing the vectors that are push_back-d
  void clearOutput();

  // Calculating new time step based on controls and rules (tank filling etc.)
  double newHydraulicTimeStep();

  // controlling the series calculations
  double endTime=0., startTime=0.;
  // time stuff
  double time, clockTime, hydraulicTimeStep = 60.;

  // get dt
  double getHydraulicTimeStep();

  // info about the seires calculation
  void seriesInfo();
    
  // setting back to original status, clearing vectors, setting pool to initial level
  void seriesInitialization();
  
private:
  // loading the series settings from INP file
  void loadTimeSettings();

  bool checkCondition(int i, double dt, double dwl);

  // converting string to a boolean operator between val1 and val2, e.g. val1 op val2 -> val1 > val2
  bool booleanWithString(string op, double left, double right);

  // Converting time (double) to a string in for of xy:as
  string secondsToTime(double seconds);

  // Converting everything to seconds like 24:00 or 7.0 DAYS
  double timeToSeconds(string s1, string s2);

  // time stuff
  double patternTimeStep, hydraulicTimeStepOriginal = 60.;

  // demand patterns
  string basePattern;
  vector<string> patternID;
  vector<vector<double> > patternValue;

  // CONTROLS
  // TODO: do this in more general way, handling time controls as well (maybe struct or class)
  vector<string> controlEdgeID; // controlled edge ID
  vector<int> controlEdgeIndex; // controlled edge INDEX in edges list
  vector<bool> controlStatus; // Open (true) or closed (false)
  vector<string> controlNodeID; // checked node (or tank) pressure ID
  vector<int> controlNodeIndex; // checked node (or tank) pressure INDEX
  vector<string> controlType; // Pool or Node
  vector<bool> controlAbove; // true if ABOVE, false if BELOW
  vector<double> controlValue; // head [m]

  // RULES
  vector<Rule*> rules;

  // filling up vectorPatternIndex in Nodes
  void patternIDtoIndex();
};

#endif