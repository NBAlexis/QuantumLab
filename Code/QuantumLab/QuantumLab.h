//=============================================================================
// FILENAME : QuantumLab.h
// 
// DESCRIPTION:
// This is the one header file for all
//
// REVISION: [dd/mm/yy]
//  [09/09/2022 nbale]
//=============================================================================

#ifndef _QUANTUMLAB_H_
#define _QUANTUMLAB_H_

#include "Base/QLSetup.h"
#include "Base/QLDefine.h"
#include "Base/PlatformIncs.h"
#include "Base/PlatformDefines.h"

#if defined(_QL_WIN)
#   if !defined(QLAPI)
#        define __LIB_TITLE__    "QuantumLab"
#       ifdef _QL_PRIVATE
#           define QLAPI __DLL_EXPORT
#       else
#           define QLAPI __DLL_IMPORT
#       endif
#       ifndef _QL_PRIVATE
#           ifdef _QL_DEBUG
#                define __LIB_FILE__    __LIB_TITLE__ "_d.lib"
#            else
#                define __LIB_FILE__ __LIB_TITLE__ ".lib"
#            endif
#            pragma __IMPORT_LIB(__LIB_FILE__)
#            pragma message("linking with " __LIB_FILE__ "...")
#            undef __LIB_FILE__
#            undef __LIB_TITLE__
#       endif
#   endif
#else
#    define QLGAPI  
#endif

#include "Base/STLExport.h"
#include "Base/CudaIncs.h"
#include "Base/QLFloat.h"
#include "Base/OtherComplexFunction.h"

//These are my tool classes, I want to use std but don't know how to get rid of warning 4251
//They are tested for years
#include "ToolClass/CLinkedList.h"
#include "ToolClass/STDStringFunctions.h"
#include "ToolClass/TArray.h"
#include "ToolClass/CBitFlag.h"
#include "ToolClass/CCString.h"
#include "ToolClass/CNMD5.h"
#include "ToolClass/MemStack.h"
#include "ToolClass/THashMap.h"
#include "ToolClass/Tracer.h"
#include "ToolClass/CYAMLParser.h"

#include "ClassicalTools/ReduceSum.h"
#include "ClassicalTools/QLRandom.h"
#include "ClassicalTools/QLMatrix.h"
#include "ClassicalTools/QuantumDataViewer.h"
#include "ToolClass/CSV.h"
#include "ClassicalTools/GPUKmeans/GPUKmeans.h"

#include "Circuit/QLGate.h"
#include "Circuit/Gates/SimpleGates.h"
#include "Circuit/Gates/CnUGate.h"
#include "Circuit/Gates/FRGate.h"
#include "Circuit/Gates/CSDDecompose.h"
#include "Circuit/Gates/AmplitudeEncode.h"
#include "Circuit/Gates/QuantumFFT.h"
#include "Circuit/Gates/QuantumPhaseEstimate.h"
#include "Circuit/Gates/PauliSimulate.h"
#include "Circuit/Gates/SwapTest.h"

#include "Circuit/Subroutine/HHL.h"
#include "Circuit/Subroutine/QKMeans.h"

#include "Simulator/QLSimulator.h"
#include "Simulator/QLSimulatorMatrix.h"
#include "Simulator/QLSimulatorVector.h"
#include "Simulator/QLSimulatorMeasure.h"

#ifndef _QL_PRIVATE
__USE_NAMESPACE
#endif

#endif //#ifndef _QUANTUMLAB_H_

//=============================================================================
// END OF FILE
//=============================================================================
