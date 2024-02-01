//=============================================================================
// FILENAME : AmplitudeEncode.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#ifndef _AMPLITUDEENCODE_H_
#define _AMPLITUDEENCODE_H_

__BEGIN_NAMESPACE

static inline UINT MostSignificantPowerTwo(UINT n)
{
    UINT k = 0;
    while ((1U << k) < n)
    {
        ++k;
    }
    return k;
}

extern TArray<QLComplex> QLAPI NormalizeV(const TArray<QLComplex>& v, UINT& lenPower);

extern TArray<Real> QLAPI NormalizeVReal(const TArray<Real>& v, UINT& lenPower);

//[[deprecated("use AmplitudeEncodeOneVector instead")]]
extern QLGate QLAPI AmplitudeEncode(const TArray<QLComplex>& v);

//[[deprecated("use AmplitudeEncodeOneVectorReal instead")]]
extern QLGate QLAPI AmplitudeEncodeReal(const TArray<Real>& v);

//=============================================================================
// rewrite amplitude encode
//=============================================================================             


/**
* (factor, r1 + r2 I, r3 + r3 I, ...)
*/
extern void QLAPI NormalizeRealToComplex(const Real* deviceData, QLComplex* deviceNormalized, UINT realVectorLength, UINT complexVectorLength, UINT uiVectorCount, Real factor);

/**
* v is on host
* phaseBuffer and deviceZ can be NULL if v is real
* len(v) = 1 << vectorPower
* 
* absBuffer and vectorCount are middle result
*/
extern void QLAPI CalculateDegrees(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ);
extern void QLAPI CalculateDegreesReal(const QLComplex* deviceV, Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY);

/**
* suppose abs and phase are already calculated
*/
extern void QLAPI CalculateDegrees(Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ);
extern void QLAPI CalculateDegreesForEach(Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ);
extern void QLAPI CalculateDegreesReal(Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY);
extern void QLAPI CalculateDegreesRealForEach(Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY);

/**
*
*/
extern void QLAPI CalculateDegrees(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ, Real* hostY, Real* hostZ);
extern void QLAPI CalculateDegreesReal(const QLComplex* deviceV, Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* hostY);

/**
*
*/
extern QLGate QLAPI ExchangeToYZGate(UINT vectorPower, Real* hostY, Real* hostZ, UBOOL bHasPhase);
extern QLGate QLAPI ExchangeToYGate(UINT vectorPower, Real* hostY);

/**
*
*/
extern QLGate QLAPI ExchangeToYZGate(UINT vectorCountPower, UINT vectorPower, Real* hostY, Real* hostZ, UBOOL bHasPhase);
extern QLGate QLAPI ExchangeToYGate(UINT vectorCountPower, UINT vectorPower, Real* hostY);

/**
* the gate to build |phi>
* Assume hostV has 2-power length
* Assume hostV is already normalized
*/
extern QLGate QLAPI AmplitudeEncodeOneVector(const QLComplex* hostv, UINT vectorPower, UBOOL bHasPhase);
extern QLGate QLAPI AmplitudeEncodeOneVectorReal(const QLComplex* hostv, UINT vectorPower);

/**
* the gate to build |i>|phi_i>
* Assume len(i) and len(phi_i) both have 2-power length
* Assume phi_i are normalized
* |i> use higher qubits, and phi use lower qubits
* For example, 8 phi_i with len(phi) = 4, there are 5 qubits, |i> is 2,3,4 and |phi> is 0,1
*/
extern QLGate QLAPI AmplitudeEncodeVectors(const QLComplex* hostv, UINT vectorCountPower, UINT vectorPower, UBOOL bHasPhase);


__END_NAMESPACE


#endif //#ifndef _AMPLITUDEENCODE_H_

//=============================================================================
// END OF FILE
//=============================================================================