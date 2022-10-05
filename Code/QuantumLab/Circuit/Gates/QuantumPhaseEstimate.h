//=============================================================================
// FILENAME : QuantumPhaseEstimate.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [02/10/2022 nbale]
//=============================================================================

#ifndef _QUANTUMPHASEESTIMATE_H_
#define _QUANTUMPHASEESTIMATE_H_

__BEGIN_NAMESPACE

/**
* assume ugate has bits from 0 to n
* not tested!
*/
extern QLGate QLAPI QuantumPhaseEstimateWithU(const QLGate& ugate, BYTE numberOfPhaseQubit);

extern QLGate QLAPI QuantumPhaseEstimateWithHSimple(const QLMatrix& h, Real t, UINT iTrotterStep, BYTE numberOfPhaseQubit, Real fMinimalKept = F(0.000001));

__END_NAMESPACE


#endif //#ifndef _QUANTUMPHASEESTIMATE_H_

//=============================================================================
// END OF FILE
//=============================================================================