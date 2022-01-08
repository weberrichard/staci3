STACI_DIR=../bin/
FUNC=TopologyGenerator

CXX=clang++
CXXFLAGS= -std=c++11 -Ofast -static-libgcc -static-libstdc++ -Wall -pedantic -I/usr/include/python2.7 -I/usr/local/include/igraph

OBJS += \
$(STACI_DIR)BasicFileIO.d \
$(STACI_DIR)Edge.d \
$(STACI_DIR)FlowMeter.d \
$(STACI_DIR)HydraulicSolver.d \
$(STACI_DIR)Graph.d \
$(STACI_DIR)IOinp.d \
$(STACI_DIR)IOxml.d \
$(STACI_DIR)Leakage.d \
$(STACI_DIR)Node.d \
$(STACI_DIR)Pipe.d \
$(STACI_DIR)Pool.d \
$(STACI_DIR)PressurePoint.d \
$(STACI_DIR)Pump.d \
$(STACI_DIR)Sensitivity.d \
$(STACI_DIR)Shutdown.d \
$(STACI_DIR)Staci.d \
$(STACI_DIR)Statistic.d \
$(STACI_DIR)Valve.d \
$(STACI_DIR)ValveFCV.d \
$(STACI_DIR)ValveISO.d \
$(STACI_DIR)ValvePRV.d \
$(STACI_DIR)ValvePSV.d \
$(STACI_DIR)ValveTCV.d \
$(STACI_DIR)Vulnerability.d \
$(STACI_DIR)xmlParser.d \
$(FUNC).o

%.d: ../../%.cpp
	@echo '[*] Building file: $<'
	@echo '[*] Invoking: CLANG++ Compiler'
	$(CXX) $(CXXFLAGS) -c -fmessage-length=0 -MMD -MP -MF "$(@:%.d=%.d)" -MT "$(@:%.d=%.d)" -o "$@"  "$<"
	@echo '[*] Finished building: $<'
	@echo ' '

all: $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $(FUNC) $(OBJS)

clean:
	-rm $(STACI_DIR)*.o $(STACI_DIR)*.d $(FUNC) $(FUNC).o