# STACI
Standard hydraulic solver for water distribution networks developed by the Department of Hydrodynamic Systems, Faculty of Mechanical Engineering, Budapest University of Technology and Economics. 

### Current projects
- *carotis:* analysing only the carotis stenosis with only 4 branches and nonlinear resistance
- *cerebral:* the effect of the incomplete Willis-circle to the cerebral arteries during surgery of internal carotis
- *heart_modelling:* how different 0D heart models influance the arterial pressure
- *reymond_modell:* early project for building complete arterial system models
- *sensitivity:* sensitivity analysis of physiologically relevant outputs to input parameters
- *virtual_patient_database:* creating virtual patient database (VPD) that mimics the whole population physiologically properly
- *vpd_ref:* finding a reference (average) patient for the VPD with differential evolution


### How to use
The code is built upon the base folder and the *Projects* folder. While the former one includes the basic sources of the *first_blood*, the latter one contains the projects which are applying the source code. Each project has an individual make file that can compile the whole code.

```sh
$ make -f make_*.mk
```

Running *.out* file will run the simulation. The *Networks* folder must contain the *.inp* files of the model with all input data.

### Dependencies
- *C++ compiler:* first_blood uses clang++, but any general C++ compiler should work
- *Eigen:* Eigen solves linear sets of equation ensuring computational efficiency
- *make:* for compyling multiple cpp files at once

### Publications

Huzsvár, T., Wéber, R., Déllei, Á., & Hős, C. (2021). Increasing the capacity of water distribution networks using fitness function transformation. Water Research, 201. https://doi.org/10.1016/j.watres.2021.117362

Wéber, R., Huzsvár, T., & Hős, C. (2021). Vulnerability of water distribution networks with real-life pipe failure statistics. Water Supply, 00(0), 1–10. https://doi.org/10.2166/ws.2021.447

Huzsvár, T., Weber, R., & Hős, C. J. (2020). Fire and drinking water capacity enhancement in water distribution networks. Water Science and Technology: Water Supply, 20(4), 1207–1214. https://doi.org/10.2166/ws.2020.037

Wéber, R., & Hős, C. (2020). Efficient Technique for Pipe Roughness Calibration and Sensor Placement for Water Distribution Systems. Journal of Water Resources Planning and Management, 146(1). https://doi.org/10.1061/(ASCE)WR.1943-5452.0001150

Wéber, R., Huzsvár, T., & Hős, C. (2020). Vulnerability analysis of water distribution networks to accidental pipe burst. Water Research, 184. https://doi.org/10.1016/j.watres.2020.116178

Huzsvár, T., Wéber, R., & Hős, C. (2019). Analysis of the segment graph of water distribution networks. Periodica Polytechnica Mechanical Engineering, 64(4), 295–300. https://doi.org/10.3311/PPme.13739