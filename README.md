# STACI
Standard hydraulic solver for water distribution networks developed by the Department of Hydrodynamic Systems, Faculty of Mechanical Engineering, Budapest University of Technology and Economics. *STACI* uses a general mathematical solver for nonlinear, algebric equations, namely Newton's iteration, thus there is no restriciton for the nature of the equations. It also has a modular built, any extension/modification can be easily performed. The figure depicts the class hiearchy below.

![Alt text](ClassHiearchy.png?raw=true "Title")

### Functions
- *snapshot simulation:* standard hydraulic simulation for a time instant for a water distribution network
- *extended period simulation:* full day/week/month simulation with demand patterns, controls, rules, active elements and actions
- *leakage modelling:* pressure dependent leakage modelling with arbitrary constants
- *pressure dependent demand:* pressure dependent demand modelling with arbitrary constants
- *shutdown plan:* creating isolation plans for pipe failures, calculating the hydraulics during reconstruction
- *vulnerability/criticality analysis:* determining the exposed segments/valves in the network using hydraulic simulations
- *waterage/chlorine/biofilm:* solving general transport equations to calculate water age/chlorine/biofilm distribution, the RK56 can solve for any source term, i.e. with any chlorine/biofilm model

### Current projects
- *biofilm:* modelling the biofilm distribution in real-life water distribution networks
- *criticality:* calculating the importance of the operation of isolation valves in terms of possible service outage, approxing with complex network theory
- *isolation valve placement:* how to place the isolation valves to minimize the possible service outage using network theory, hydraulics, and NSGA-II
- *leakage reduction:* PRVs can reduce the leakge amount by decreasing the average pressure in the network, the question is where is their optimal place and setting
- *vulnerability:* determining how certain segments are exposed to an accidental pipe burst

### How to use
The code is built upon the base folder and the *Projects* folders. While the former one includes the basic sources of the *STACI*, the latter one contains the projects which are applying the source code. Each project has an individual make file that can compile the whole code. The *Plot* folder contains Matlab scripts for visualisation.

```sh
$ make -f make_*.mk
```

Running *.out* file will run the simulation. The *Networks* folder must contain the *.inp* files of the model with all input data.

### Dependencies
- *C++ compiler:* first_blood uses clang++, but any general C++ compiler should work
- *Eigen:* Eigen solves linear sets of equation ensuring computational efficiency
- *make:* for compyling multiple cpp files at once

### Developement team
Dr Richárd Wéber, assistant professor
Tamás Huzsvár, PhD student

### Publications

Huzsvár, T., Wéber, R., Déllei, Á., & Hős, C. (2021). Increasing the capacity of water distribution networks using fitness function transformation. Water Research, 201. https://doi.org/10.1016/j.watres.2021.117362

Wéber, R., Huzsvár, T., & Hős, C. (2021). Vulnerability of water distribution networks with real-life pipe failure statistics. Water Supply, 00(0), 1–10. https://doi.org/10.2166/ws.2021.447

Huzsvár, T., Weber, R., & Hős, C. J. (2020). Fire and drinking water capacity enhancement in water distribution networks. Water Science and Technology: Water Supply, 20(4), 1207–1214. https://doi.org/10.2166/ws.2020.037

Wéber, R., & Hős, C. (2020). Efficient Technique for Pipe Roughness Calibration and Sensor Placement for Water Distribution Systems. Journal of Water Resources Planning and Management, 146(1). https://doi.org/10.1061/(ASCE)WR.1943-5452.0001150

Wéber, R., Huzsvár, T., & Hős, C. (2020). Vulnerability analysis of water distribution networks to accidental pipe burst. Water Research, 184. https://doi.org/10.1016/j.watres.2020.116178

Huzsvár, T., Wéber, R., & Hős, C. (2019). Analysis of the segment graph of water distribution networks. Periodica Polytechnica Mechanical Engineering, 64(4), 295–300. https://doi.org/10.3311/PPme.13739