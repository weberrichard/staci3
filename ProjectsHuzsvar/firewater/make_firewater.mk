
STACI_DIR=/mnt/d/Linux/Staci_V3/staci/Staci@HuT/
STACI_SRC=/mnt/d/Linux/Staci_V3/staci/
FUNC=firewater

CXX=clang++
CXXFLAGS=-g -std=c++11 -Wall -pedantic -I/usr/include/python2.7 -I/usr/local/include/igraph -L/usr/local/lib -ligraph

# Add inputs and outputs from these tool invocations to the build variables
O_SRCS += \
$(STACI_DIR)Edge.o \
$(STACI_DIR)AnyOption.o \
$(STACI_DIR)Weir.o \
$(STACI_DIR)Channel.o \
$(STACI_DIR)Pipe.o \
$(STACI_DIR)Node.o \
$(STACI_DIR)Valve.o \
$(STACI_DIR)PressurePoint.o \
$(STACI_DIR)Staci.o \
$(STACI_DIR)Pump.o \
$(STACI_DIR)Pool.o \
$(STACI_DIR)CheckValve.o \
$(STACI_DIR)IOxml.o \
$(STACI_DIR)lubksb.o \
$(STACI_DIR)ludcmp.o \
$(STACI_DIR)xmlParser.o \
$(STACI_DIR)Statistic.o \
$(STACI_DIR)HydraulicSolver.o \
$(STACI_DIR)Sensitivity.o \
$(STACI_DIR)Transport.o \
$(FUNC).o 

CPP_SRCS += \
$(STACI_SRC)/Edge.cpp \
$(STACI_SRC)/AnyOption.cpp \
$(STACI_SRC)/Weir.cpp \
$(STACI_SRC)/Channel.cpp \
$(STACI_SRC)/Pipe.cpp \
$(STACI_SRC)/Node.cpp \
$(STACI_SRC)/Valve.cpp \
$(STACI_SRC)/PressurePoint.cpp \
$(STACI_SRC)/Staci.cpp \
$(STACI_SRC)/Pump.cpp \
$(STACI_SRC)/Pool.cpp \
$(STACI_SRC)/CheckValve.cpp \
$(STACI_SRC)/IOxml.cpp \
$(STACI_SRC)/lubksb.cpp \
$(STACI_SRC)/ludcmp.cpp \
$(STACI_SRC)/xmlParser.cpp \
$(STACI_SRC)/Statistic.cpp \
$(STACI_SRC)/HydraulicSolver.cpp \
$(STACI_SRC)/Sensitivity.cpp \
$(STACI_SRC)/Transport.cpp \
$(FUNC).cpp

OBJS += \
$(STACI_DIR)Edge.o \
$(STACI_DIR)AnyOption.o \
$(STACI_DIR)Weir.o \
$(STACI_DIR)Channel.o \
$(STACI_DIR)Pipe.o \
$(STACI_DIR)Node.o \
$(STACI_DIR)Valve.o \
$(STACI_DIR)PressurePoint.o \
$(STACI_DIR)Staci.o \
$(STACI_DIR)Pump.o \
$(STACI_DIR)Pool.o \
$(STACI_DIR)CheckValve.o \
$(STACI_DIR)IOxml.o \
$(STACI_DIR)lubksb.o \
$(STACI_DIR)ludcmp.o \
$(STACI_DIR)xmlParser.o \
$(STACI_DIR)Statistic.o \
$(STACI_DIR)HydraulicSolver.o \
$(STACI_DIR)Sensitivity.o \
$(STACI_DIR)Transport.o \
$(FUNC).o

CPP_DEPS += \
$(STACI_DIR)Edge.d \
$(STACI_DIR)AnyOption.d \
$(STACI_DIR)Weir.d \
$(STACI_DIR)Channel.d \
$(STACI_DIR)Pipe.d \
$(STACI_DIR)Node.d \
$(STACI_DIR)Valve.d \
$(STACI_DIR)PressurePoint.d \
$(STACI_DIR)Staci.d \
$(STACI_DIR)Pump.d \
$(STACI_DIR)Pool.d \
$(STACI_DIR)CheckValve.d \
$(STACI_DIR)IOxml.d \
$(STACI_DIR)lubksb.d \
$(STACI_DIR)ludcmp.d \
$(STACI_DIR)xmlParser.d \
$(STACI_DIR)/Statistic.d \
$(STACI_DIR)/HydraulicSolver.d \
$(STACI_DIR)/Sensitivity.d \
$(STACI_DIR)/Transport.d \
$(FUNC).d


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I/usr/local/include -I/usr/local/include/eigen3 -O2 -g -Wall -c -fmessage-length=0 -pedantic -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@"  "$<"
	@echo 'Finished building: $<'
	@echo ' '


all: $(OBJS) $(USER_OBJS)
	$(CXX) $(CXXFLAGS) -o $(FUNC) $(OBJS) $(USER_OBJS) -L/usr/local/lib -L/usr/lib -L/usr/local/lib -lm -lumfpack -lpython2.7

clean:
	-rm $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) $(FUNC)
