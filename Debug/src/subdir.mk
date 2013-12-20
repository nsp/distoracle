################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/qt.cpp 

CU_SRCS += \
../src/sssp.cu 

CU_DEPS += \
./src/sssp.d 

OBJS += \
./src/qt.o \
./src/sssp.o 

CPP_DEPS += \
./src/qt.d 


# Each subdirectory must supply rules for building sources it contributes
src/qt.o: ../src/qt.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"src/qt.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/bin/nvcc -G -g -O0 -gencode arch=compute_20,code=sm_20 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/bin/nvcc --compile -G -O0 -g -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


