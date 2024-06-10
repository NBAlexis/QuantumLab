//=============================================================================
// FILENAME : SwapTest.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [14/10/2022 nbale]
//=============================================================================

#ifndef _SWAPTEST_H_
#define _SWAPTEST_H_

__BEGIN_NAMESPACE


extern QLGate QLAPI CreateSwapTest(TArray<QLComplex> v1, TArray<QLComplex> v2);

/**
* when v1 and v2 are real vectors
*/
extern QLGate QLAPI CreateSwapTestReal(TArray<QLComplex> v1, TArray<QLComplex> v2);

extern QLGate QLAPI ZeroTest(TArray<QLComplex> v1, TArray<QLComplex> v2);

/**
* when v1 and v2 are real vectors
*/
extern QLGate QLAPI ZeroTestReal(TArray<QLComplex> v1, TArray<QLComplex> v2);

extern QLGate QLAPI SimpleZeroTest(TArray<QLComplex> v1, TArray<QLComplex> v2, BYTE byEncodeBits);

__END_NAMESPACE


#endif //#ifndef _SWAPTEST_H_

//=============================================================================
// END OF FILE
//=============================================================================