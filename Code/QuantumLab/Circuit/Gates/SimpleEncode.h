//=============================================================================
// FILENAME : SimpleEncode.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [05/06/2024 nbale]
//=============================================================================

#ifndef _SIMPLEENCODE_H_
#define _SIMPLEENCODE_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(ELinkStyle,
    Full,
    Linear,
    Circular,
    PairWise,
    DoublePairWise,
    SCA
)

__DEFINE_ENUM(ELinkLayer,
    CX,
    CZ,
    CRX,
)

/**
* ry is put to hostv[n].x
* rz is put to hostv[n].y
* 
* it is ordered as:
* 0 CONT 3
* 1      4
* 2      ...
* 
* 
* 2209.12788
* the Hadamad is removed.
*/
extern QLGate QLAPI SimpleEncodeOneVector(const QLComplex* hostv, BYTE qubits, UINT uiVLength);

/**
*
*/
extern QLGate QLAPI SimpleEncodeVectors(const QLComplex* hostv, BYTE vectorCountPower, BYTE qubits, UINT uiVLength);

extern QLGate QLAPI SimpleEncodeOneVectorWithLinkStype(const QLComplex* hostv, ELinkStyle eStype, ELinkLayer eLinkLayer, BYTE qubits, UINT uiVLength);

extern QLGate QLAPI SimpleEncodeVectorsWithLinkStype(const QLComplex* hostv, ELinkStyle eStype, ELinkLayer eLinkLayer, BYTE vectorCountPower, BYTE qubits, UINT uiVLength);


extern void QLAPI SimpleEncodeAddLinkLayerAll(ELinkStyle eStype, ELinkLayer eLinkLayer, QLGate& gate, UINT level, const TArray<BYTE>& linkQubits);

__END_NAMESPACE


#endif //#ifndef _SIMPLEENCODE_H_

//=============================================================================
// END OF FILE
//=============================================================================