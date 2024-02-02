//=============================================================================
// FILENAME : QuantumDataViewer.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [08/10/2022 nbale]
//=============================================================================

#ifndef _QUNATUMDATAVIEWER_H_
#define _STATEVECTORVIEWER_H_

__BEGIN_NAMESPACE

/**
* if index = 2, it is to be viewed.
* for example, if the index list is 1 2 1 0 2, it means |20121> = |?01?1> (note that 0qubit is on the right)
* it is to view the data with qubit0=1, qubit2=1, qubit3=0
*/
extern QLMatrix QLAPI ShowStateVectorDetail(const QLComplex* statedata, TArray<BYTE> index, UBOOL bNormalize = TRUE);

extern QLMatrix QLAPI ShowStateVectorDetail(const struct Qureg& vec, TArray<BYTE> index, UBOOL bNormalize = TRUE);

extern QLMatrix QLAPI ShowStateVectorDetail(const QLComplex* statedata, BYTE count, BYTE idx0, ...);

extern QLMatrix QLAPI ShowStateVectorDetail(const struct Qureg& vec, BYTE count, BYTE idx0, ...);

extern QLMatrix QLAPI StateToMatrix(const struct Qureg& vec);

extern void QLAPI HostBufferViewer(const Real* buffer, UINT w, UINT h);

extern void QLAPI HostBufferViewer(const QLComplex* buffer, UINT w, UINT h);

extern void QLAPI DeviceBufferViewer(const Real* buffer, UINT w, UINT h);


__END_NAMESPACE


#endif //#ifndef _STATEVECTORVIEWER_H_

//=============================================================================
// END OF FILE
//=============================================================================