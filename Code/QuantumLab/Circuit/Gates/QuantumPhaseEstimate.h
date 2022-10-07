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

/**
* ugate on 0 to n-1
* n is ancilla for ugate
* from n+1, are phase qubits
*/
extern QLGate QLAPI QuantumPhaseEstimateWithHSimple(const QLMatrix& h, Real t, UINT iTrotterStep, BYTE numberOfPhaseQubit, Real fMinimalKept = F(0.000001));

extern QLGate QLAPI BuildImprovedQPEInitialState(BYTE byQubit);

/**
* |tau> <tau| exp(i H tau t0 / T)
* The shortest interval is 0.5 t0 / T
* 
* ugate on 0 to n-1
* n is ancilla for ugate
* from n+1, are tau qubits
*/
extern QLGate QLAPI ConditionalHamiltonianEvolution(const QLMatrix& h, BYTE phaseQubitNum, Real t0, UINT shortestIntervalTrotter, Real fMinimalKept = F(0.000001));

/**
* using the method discribed in
* 10.1103/PhysRevLett.103.150502
* The eignvalues close to 2pi does not show up
*/
extern QLGate QLAPI QuantumPhaseEstimateWithHImproved(const QLMatrix& h, Real t, UINT iTrotterStep, BYTE numberOfPhaseQubit, Real fMinimalKept = F(0.000001));

__END_NAMESPACE


#endif //#ifndef _QUANTUMPHASEESTIMATE_H_

//=============================================================================
// END OF FILE
//=============================================================================