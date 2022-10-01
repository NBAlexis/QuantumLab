//=============================================================================
// FILENAME : QuantumFFT.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#ifndef _QUANTUMFFT_H_
#define _QUANTUMFFT_H_

__BEGIN_NAMESPACE

/**
* same as Fourier[v] of Mathematica
* same as ifft(v) / |ifft(v)| of numpy and Cuda
*/
extern QLGate QLAPI QuantumFFTGate(BYTE qubit);

__END_NAMESPACE


#endif //#ifndef _QUANTUMFFT_H_

//=============================================================================
// END OF FILE
//=============================================================================