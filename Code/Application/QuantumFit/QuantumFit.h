//=============================================================================
// FILENAME : QuantumFit.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/10/2022 nbale]
//=============================================================================

#ifndef _QUANTUMFIT_H_
#define _QUANTUMFIT_H_

#include "QuantumLab.h"

typedef QLComplex (*fitFunction_t)(const QLComplex& );

inline QLComplex one(const QLComplex& c)
{
    return _mcr(1.0);
}

inline QLComplex x(const QLComplex& c)
{
    return c;
}

inline QLComplex x2(const QLComplex& c)
{
    return _cuCmulf(c, c);
}

inline QLComplex expm(const QLComplex& c)
{
    return __cuCexpf(_make_cuComplex(-c.x, -c.y));
}

/**
* xy is expected as a two-line matrix, for data to be fitted
* 
*/
extern QLGate BuildQuantumFitCircuit(const QLMatrix& xy, const TArray<fitFunction_t>& fitFunctions, UINT trotter, BYTE phaseQubitNum, Real maxAbsoluteEigen, Real fTruncate);

extern void SimulateQuantumFit(const QLGate& gate, UINT size, BYTE phaseQubit);

extern QLMatrix LoadAMatrix(const CCString& sFileName);

extern void BuildIFUsingA(QLMatrix a, QLMatrix& mif, QLMatrix& mifd);

extern QLGate BuildQuantumFitCircuitWithAY(const QLMatrix& ay, UINT trotter, BYTE phaseQubitNum, Real maxAbsoluteEigen, Real fTruncate);

#endif //#ifndef _QUANTUMFIT_H_

//=============================================================================
// END OF FILE
//=============================================================================