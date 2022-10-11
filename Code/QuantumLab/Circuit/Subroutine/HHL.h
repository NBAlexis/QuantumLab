//=============================================================================
// FILENAME : HHL.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [07/10/2022 nbale]
//=============================================================================

#ifndef _HHL_H_
#define _HHL_H_

__BEGIN_NAMESPACE

extern QLGate QLAPI HHLGate(const QLMatrix& h, const TArray<QLComplex>& y, UINT trotterStep, Real maxAbsEigenValue, BYTE phaseQubitNum, Real small = F(0.00001));

extern QLGate QLAPI HermitianMatrixMultiply(const QLMatrix& h, const TArray<QLComplex>& y, UINT trotterStep, Real maxAbsEigenValue, BYTE phaseQubitNum, Real small = F(0.00001));

extern QLGate QLAPI MatrixPowerGate(const QLMatrix& h, INT iPower, UINT trotterStep, Real maxAbsEigenValue, BYTE phaseQubitNum, Real small = F(0.00001));

__END_NAMESPACE


#endif //#ifndef _HHL_H_

//=============================================================================
// END OF FILE
//=============================================================================