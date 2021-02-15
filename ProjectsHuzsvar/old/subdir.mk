################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../Agelem.o \
../AnyOption.o \
../BukoMutargy.o \
../Csatorna.o \
../Cso.o \
../Csomopont.o \
../JelleggorbesFojtas.o \
../KonstNyomas.o \
../Staci.o \
../Staci_Runner.o \
../Szivattyu.o \
../Vegakna.o \
../VisszacsapoSzelep.o \
../data_io.o \
../lubksb.o \
../ludcmp.o \
../staci_main.o \
../xmlParser.o 

CPP_SRCS += \
../Agelem.cpp \
../AnyOption.cpp \
../BukoMutargy.cpp \
../Csatorna.cpp \
../Cso.cpp \
../Csomopont.cpp \
../JelleggorbesFojtas.cpp \
../KonstNyomas.cpp \
../Staci.cpp \
../Szivattyu.cpp \
../Vegakna.cpp \
../VisszacsapoSzelep.cpp \
../data_io.cpp \
../lubksb.cpp \
../ludcmp.cpp \
../staci_main.cpp \
../xmlParser.cpp 

OBJS += \
./Agelem.o \
./AnyOption.o \
./BukoMutargy.o \
./Csatorna.o \
./Cso.o \
./Csomopont.o \
./JelleggorbesFojtas.o \
./KonstNyomas.o \
./Staci.o \
./Szivattyu.o \
./Vegakna.o \
./VisszacsapoSzelep.o \
./data_io.o \
./lubksb.o \
./ludcmp.o \
./staci_main.o \
./xmlParser.o 

CPP_DEPS += \
./Agelem.d \
./AnyOption.d \
./BukoMutargy.d \
./Csatorna.d \
./Cso.d \
./Csomopont.d \
./JelleggorbesFojtas.d \
./KonstNyomas.d \
./Staci.d \
./Szivattyu.d \
./Vegakna.d \
./VisszacsapoSzelep.d \
./data_io.d \
./lubksb.d \
./ludcmp.d \
./staci_main.d \
./xmlParser.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


