//=============================================================================
// FILENAME : QLSetup.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/09/2022 nbale]
//=============================================================================

#ifndef _QLSETUP_H_
#define _QLSETUP_H_

#ifdef _DEBUG
#define _QL_DEBUG 1
#endif

//This is the tag for windows, msvc specific
#ifdef WIN64
#define _QL_WIN 1
#endif

#ifndef _QL_DOUBLEFLOAT
#define _QL_DOUBLEFLOAT 1
#endif

//_QL_USE_LAUNCH_BOUND = 0 or 1.
//NOTE: If the regcount required is out-numbered, sometimes, there is NO error message!
//So, either be sure to build with _QL_USE_LAUNCH_BOUND = 1, or reduce the thread count
//reduce the thread count is expansive, so _QL_USE_LAUNCH_BOUND = 1 is recommanded
//It's better to complie using the maximum thread per block of the device of the computer.
#if _QL_DEBUG
#define _QL_USE_LAUNCH_BOUND 0
#else
#define _QL_USE_LAUNCH_BOUND 1
#endif

#ifndef _QL_LAUNCH_MAX_THREAD
#if _QL_DEBUG
#define _QL_LAUNCH_MAX_THREAD 512
#define _QL_LAUNCH_MAX_THREADHALF 256
#else
#define _QL_LAUNCH_MAX_THREAD 1024
#define _QL_LAUNCH_MAX_THREADHALF 512
#endif
#endif


#endif //#ifndef _QLSETUP_H_

//=============================================================================
// END OF FILE
//=============================================================================