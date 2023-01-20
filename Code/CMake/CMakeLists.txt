cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

if (DEFINED NVCCROOT)
    set(CMAKE_CUDA_COMPILER ${NVCCROOT})
    MESSAGE("CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
endif()

set(CUDA_CMP "compute_61")
set(CUDA_SM "sm_61")
if (DEFINED CUDASM)
    set(CUDA_CMP "compute_${CUDASM}")
    set(CUDA_SM "sm_${CUDASM}")
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set (CMAKE_CUDA_ARCHITECTURES ${CUDASM})
    endif(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
endif()

MESSAGE("Note: arch is ${CUDA_CMP} and ${CUDA_SM}.")
MESSAGE("52 for 970, 61 for GTX10, 70 for V100, 75 for RTX20, RTX16, 86 for RTX30")
project(QuantumLabProj LANGUAGES C CXX CUDA)

set(CMAKE_GENERATOR_PLATFORM x64)

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/Ubuntu)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})
# This is our code file dir
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
# Flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-class-memaccess")
add_definitions(-D_UBUNTU)
MESSAGE("CMAKE_CUDA_FLAGS flag = ${CMAKE_CUDA_FLAGS}")
MESSAGE("CMAKE_CXX_FLAGS flag = ${CMAKE_CXX_FLAGS}")




# ==================== 
# QuEST 
# =================

include_directories(${PROJECT_SOURCE_DIR}/QuEST)
add_library(QuEST STATIC
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/include/PlatformFunctions.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/include/QuEST.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/include/QuEST_complex.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/include/QuEST_precision.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/include/QuEST_private.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/CPU/QuEST_cpu_internal.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/mt19937ar.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_debug.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_internal.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_qasm.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_validation.h
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/GPU/QuEST_gpu.cu
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/CPU/QuEST_cpu.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/CPU/QuEST_cpu_distributed.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/CPU/QuEST_cpu_local.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/mt19937ar.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_common.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_private.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_qasm.c
    ${PROJECT_SOURCE_DIR}/QuEST340/QuEST/src/QuEST_validation.c
    )

# Request that QuEST be built with -std=c++14
# As this is a public compile feature anything that links to 
# QuEST will also build with -std=c++14
target_compile_features(QuEST PUBLIC cxx_std_14)
 
# We need to explicitly state that we need all CUDA files in the 
# QuEST library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties(QuEST PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(QuEST PRIVATE ${CMAKE_SOURCE_DIR}/../QuEST340/QuEST/include)


target_include_directories(QuEST PRIVATE ${CMAKE_SOURCE_DIR}/../QuEST340/QuEST/src)



target_link_libraries(QuEST -lcurand)
target_link_libraries(QuEST -lcufft)
target_link_libraries(QuEST -lcusolver)
target_link_libraries(QuEST -lcublas)

# To enable the double, the minimum arch is 6.0
target_compile_options(QuEST PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=${CUDA_CMP},code=${CUDA_SM}>)




# ==================== 
# QuantumLab 
# =================

include_directories(${PROJECT_SOURCE_DIR}/QuantumLab)
add_library(QuantumLab STATIC
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/CudaIncs.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/OtherComplexFunction.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/PlatformDefines.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/PlatformIncs.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/QLDefine.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/QLFloat.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/QLSetup.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Base/STLExport.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/AmplitudeEncode.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/CnUGate.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/CSDDecompose.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/FRGate.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/PauliSimulate.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/QuantumFFT.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/QuantumPhaseEstimate.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/SimpleGates.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/SwapTest.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/QLGate.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Subroutine/HHL.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Subroutine/QKMeans.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/GPUKmeans/GPUKmeans.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/QLMatrix.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/QLRandom.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/QuantumDataViewer.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/ReduceSum.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/QuantumLab.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/QuantumLabPCH.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulator.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulatorMatrix.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulatorMeasure.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulatorVector.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CBitFlag.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CCString.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CLinkedList.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CNMD5.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CSV.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CYAMLParser.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/MemStack.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/STDStringFunctions.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/TArray.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/THashMap.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/Tracer.h
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Subroutine/QKMeans.cu
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/GPUKmeans/GPUKmeans.cu
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/QLMatrix.cu
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/QLRandom.cu
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/ReduceSum.cu
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/AmplitudeEncode.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/CnUGate.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/CSDDecompose.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/FRGate.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/PauliSimulate.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/QuantumFFT.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/QuantumPhaseEstimate.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/SimpleGates.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Gates/SwapTest.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/QLGate.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Circuit/Subroutine/HHL.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/ClassicalTools/QuantumDataViewer.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/QuantumLabPCH.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulator.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulatorMatrix.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulatorMeasure.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/Simulator/QLSimulatorVector.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CCString.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CSV.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/CYAMLParser.cpp
    ${PROJECT_SOURCE_DIR}/QuantumLab/ToolClass/Tracer.cpp
    )

# Request that QuantumLab be built with -std=c++14
# As this is a public compile feature anything that links to 
# QuantumLab will also build with -std=c++14
target_compile_features(QuantumLab PUBLIC cxx_std_14)
 
# We need to explicitly state that we need all CUDA files in the 
# QuantumLab library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties(QuantumLab PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(QuantumLab PRIVATE ${CMAKE_SOURCE_DIR}/../QuEST340/QuEST/include)



target_link_libraries(QuantumLab -lcurand)
target_link_libraries(QuantumLab -lcufft)
target_link_libraries(QuantumLab -lcusolver)
target_link_libraries(QuantumLab -lcublas)
target_link_libraries(QuantumLab QuEST)

# To enable the double, the minimum arch is 6.0
target_compile_options(QuantumLab PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=${CUDA_CMP},code=${CUDA_SM}>)




# ==================== 
# SimpleTest 
# =================

include_directories(${PROJECT_SOURCE_DIR}/Application/SimpleTest)
add_executable(SimpleTest 
    ${PROJECT_SOURCE_DIR}/Application/SimpleTest/SimpleTest.cpp
    )

target_compile_features(SimpleTest PUBLIC cxx_std_14)
target_link_libraries(SimpleTest QuantumLab)



# ==================== 
# QKMeans 
# =================

include_directories(${PROJECT_SOURCE_DIR}/Application/QKMeans)
add_executable(QKMeans 
    ${PROJECT_SOURCE_DIR}/Application/QKMeans/QKMeans.h
    ${PROJECT_SOURCE_DIR}/Application/QKMeans/ClassicalKmeans.cpp
    ${PROJECT_SOURCE_DIR}/Application/QKMeans/QKMeans.cpp
    ${PROJECT_SOURCE_DIR}/Application/QKMeans/QuantumKmeans.cpp
    )

target_compile_features(QKMeans PUBLIC cxx_std_14)
target_link_libraries(QKMeans QuantumLab)
