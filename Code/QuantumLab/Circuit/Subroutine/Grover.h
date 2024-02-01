//=============================================================================
// FILENAME : Grover.h
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [11/28/2023 nbale]
//=============================================================================

#ifndef _GROVER_H_
#define _GROVER_H_

__BEGIN_NAMESPACE

/**
* 1 - 2|n><n| = diag(1,...,1,-1,1,...1) where -1 is at position n
* we need numberOfQubits + 1 qubits, where the last one is ancilla
* we asumme ancilla is at |0>-|1> already
*
* ===========================================
* 1-step apply not to encode n on the qubits
*
* 2-step apply hadard on ancilla
*
* 3-step apply cn-NOT
*
* 4-step apply not to encode n on the qubits again
*
* Assume 2 ancilla
* bit numberOfQubits:  |0>-|1> to do the phase kick back
* bit numberOfQubits+1: Zeroed ancilla for CNOT
*/
extern QLGate QLAPI GroverSXGate(BYTE numberOfQubits, UINT elementToFlip, UBOOL bUseAncilla = TRUE);

/**
* Assume A|x> = \sum _n |phi_n>|n>
* increase the ability to measure |n0> to build the state phi_{n0}
*
* where x is an integer, n is an integer
* number of qubits of x is a+b, number of qubits of phi_n is a, number of qubits of |n> is b
*
* Depend on whether ancilla is used, we need a+b+2 (with ancilla) or a+b+1 qubits.
*
* agate is the gate A|x> = \sum _n |phi_n>|n>
* subspaceQubits is qubits for |n>, with length b
* elementToFlip is the n0 to choose
* bUseAncilla is whether to use ancilla to implement Cn-NOT
* uiInitialState is |x> such that A|x> = \sum _n |phi_n>|n>, usually it should be |0>
* uiRepeat is the power of Q, where Q = Amplitude Amplification
* 
* agate will be daggered, therefore is not constant
*/
QLGate QLAPI AmplitudeAmplification(QLGate agate, TArray<BYTE> subspaceQubits, UINT elementToFlip, UINT uiRepeat, UBOOL bUseAncilla = TRUE, UINT uiInitialState = 0);

__END_NAMESPACE


#endif //#ifndef _GROVER_H_

//=============================================================================
// END OF FILE
//=============================================================================