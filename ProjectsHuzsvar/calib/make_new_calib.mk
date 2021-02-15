
STACI_DIR=/home/rweber/Dropbox/0_PhD/staci/Staci@Weber/
STACI_SRC=/home/rweber/Dropbox/0_PhD/staci/
FUNC=new_calib

CXX=clang++
CXXFLAGS=-g -std=c++11 -Wall -pedantic -I/usr/include/python2.7 -I/usr/local/include/igraph -L/usr/local/lib -ligraph

# Add inputs and outputs from these tool invocations to the build variables
O_SRCS += \
$(STACI_DIR)Agelem.o \
$(STACI_DIR)AnyOption.o \
$(STACI_DIR)BukoMutargy.o \
$(STACI_DIR)Csatorna.o \
$(STACI_DIR)Cso.o \
$(STACI_DIR)Csomopont.o \
$(STACI_DIR)JelleggorbesFojtas.o \
$(STACI_DIR)KonstNyomas.o \
$(STACI_DIR)Staci.o \
$(STACI_DIR)Szivattyu.o \
$(STACI_DIR)Vegakna.o \
$(STACI_DIR)VisszacsapoSzelep.o \
$(STACI_DIR)data_io.o \
$(STACI_DIR)lubksb.o \
$(STACI_DIR)ludcmp.o \
$(STACI_DIR)xmlParser.o \
$(STACI_DIR)Statistic.o \
$(STACI_DIR)HydraulicSolver.o \
$(STACI_DIR)Sensitivity.o \
$(STACI_DIR)Transport.o \
$(STACI_DIR)Calibration.o \
$(STACI_DIR)SVDCalibration.o \
$(FUNC).o 

CPP_SRCS += \
$(STACI_SRC)/Agelem.cpp \
$(STACI_SRC)/AnyOption.cpp \
$(STACI_SRC)/BukoMutargy.cpp \
$(STACI_SRC)/Csatorna.cpp \
$(STACI_SRC)/Cso.cpp \
$(STACI_SRC)/Csomopont.cpp \
$(STACI_SRC)/JelleggorbesFojtas.cpp \
$(STACI_SRC)/KonstNyomas.cpp \
$(STACI_SRC)/Staci.cpp \
$(STACI_SRC)/Szivattyu.cpp \
$(STACI_SRC)/Vegakna.cpp \
$(STACI_SRC)/VisszacsapoSzelep.cpp \
$(STACI_SRC)/data_io.cpp \
$(STACI_SRC)/lubksb.cpp \
$(STACI_SRC)/ludcmp.cpp \
$(STACI_SRC)/xmlParser.cpp \
$(STACI_SRC)/Statistic.cpp \
$(STACI_SRC)/HydraulicSolver.cpp \
$(STACI_SRC)/Sensitivity.cpp \
$(STACI_SRC)/Transport.cpp \
$(STACI_SRC)/Calibration.cpp \
$(STACI_SRC)/SVDCalibration.cpp \
$(FUNC).cpp

OBJS += \
$(STACI_DIR)Agelem.o \
$(STACI_DIR)AnyOption.o \
$(STACI_DIR)BukoMutargy.o \
$(STACI_DIR)Csatorna.o \
$(STACI_DIR)Cso.o \
$(STACI_DIR)Csomopont.o \
$(STACI_DIR)JelleggorbesFojtas.o \
$(STACI_DIR)KonstNyomas.o \
$(STACI_DIR)Staci.o \
$(STACI_DIR)Szivattyu.o \
$(STACI_DIR)Vegakna.o \
$(STACI_DIR)VisszacsapoSzelep.o \
$(STACI_DIR)data_io.o \
$(STACI_DIR)lubksb.o \
$(STACI_DIR)ludcmp.o \
$(STACI_DIR)xmlParser.o \
$(STACI_DIR)Statistic.o \
$(STACI_DIR)HydraulicSolver.o \
$(STACI_DIR)Sensitivity.o \
$(STACI_DIR)Transport.o \
$(STACI_DIR)Calibration.o \
$(STACI_DIR)SVDCalibration.o \
$(FUNC).o

CPP_DEPS += \
$(STACI_DIR)Agelem.d \
$(STACI_DIR)AnyOption.d \
$(STACI_DIR)BukoMutargy.d \
$(STACI_DIR)Csatorna.d \
$(STACI_DIR)Cso.d \
$(STACI_DIR)Csomopont.d \
$(STACI_DIR)JelleggorbesFojtas.d \
$(STACI_DIR)KonstNyomas.d \
$(STACI_DIR)Staci.d \
$(STACI_DIR)Szivattyu.d \
$(STACI_DIR)Vegakna.d \
$(STACI_DIR)VisszacsapoSzelep.d \
$(STACI_DIR)data_io.d \
$(STACI_DIR)lubksb.d \
$(STACI_DIR)ludcmp.d \
$(STACI_DIR)xmlParser.d \
$(STACI_DIR)/Statistic.d \
$(STACI_DIR)/HydraulicSolver.d \
$(STACI_DIR)/Sensitivity.d \
$(STACI_DIR)/Transport.d \
$(STACI_DIR)/Calibration.d \
$(STACI_DIR)/SVDCalibration.d \
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
