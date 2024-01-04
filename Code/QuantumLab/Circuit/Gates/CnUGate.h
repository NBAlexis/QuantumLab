//=============================================================================
// FILENAME : CnUGate.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [23/09/2022 nbale]
//=============================================================================

#ifndef _CNUGATE_H_
#define _CNUGATE_H_

__BEGIN_NAMESPACE

inline QLMatrix Sqrt2by2(const QLMatrix& m)
{
    const QLComplex u00 = m.Get(0, 0);
    const QLComplex u11 = m.Get(1, 1);
    const QLComplex u10 = m.Get(1, 0);
    const QLComplex u01 = m.Get(0, 1);
    const QLComplex tau = _cuCaddf(u00, u11);
    const QLComplex delta = _cuCsubf(_cuCmulf(u00, u11), _cuCmulf(u10, u01));
    const QLComplex s = __cuCsqrtf_host(delta);
    const QLComplex t = __cuCsqrtf_host(_cuCaddf(tau, cuCmulf_cr(s, F(2.0))));

    QLComplex v[] = { _cuCdivf(_cuCaddf(u00, s), t), _cuCdivf(u01, t),  _cuCdivf(u10, t),  _cuCdivf(_cuCaddf(u11, s), t) };
    return QLMatrix::CopyCreate(2, 2, v);
}

/**
* total bits is number of controller bits + 1
* The target is the last one qubit
*/
extern QLGate QLAPI CreateCnNot(BYTE numOfController);

/**
* 0,1 controller, 2 is target
* it is known that n-bit toffoli needs 2n CNOT to implement
* 10.26421/QIC9.5-6-8
* 
* T = exp(i pi sigmaz / 8)
* 
* ---------.----------.-------.---T---.---
*          |          |       |       |
* ---.-----|-----.----|---T---+---Td--+---
*    |     |     |    |
* -H-+--Td-+--T--+-Td-+---T---H-----------
*/
inline QLGate CreateToffoliGate()
{
    QLGate ret;
    ret.AddQubits(3);
    ret.m_sName = _T("toffoli");

    QLGate h(EBasicOperation::EBO_H);
    QLGate t(EBasicOperation::EBO_RZ, -PI / F(4.0));
    QLGate td(EBasicOperation::EBO_RZ, PI / F(4.0));
    QLGate cx(EBasicOperation::EBO_CX);

    ret.AppendGate(h, 2);
    ret.AppendGate(cx, 1, 2);
    ret.AppendGate(td, 2);
    ret.AppendGate(cx, 0, 2);
    ret.AppendGate(t, 2);
    ret.AppendGate(cx, 1, 2);
    ret.AppendGate(td, 2);
    ret.AppendGate(cx, 0, 2);

    ret.AppendGate(t, 2);
    ret.AppendGate(h, 2);

    ret.AppendGate(t, 1);
    ret.AppendGate(cx, 0, 1);
    ret.AppendGate(t, 0);
    ret.AppendGate(td, 1);
    ret.AppendGate(cx, 0, 1);

    return ret;
}

extern QLGate QLAPI CreateCnU(BYTE numOfController, const QLMatrix& mtr);

extern QLGate QLAPI CreateCnRX(BYTE numOfController, Real fDegree);

extern QLGate QLAPI CreateCnRY(BYTE numOfController, Real fDegree);

extern QLGate QLAPI CreateCnRZ(BYTE numOfController, Real fDegree);

extern QLGate QLAPI CreateCnP(BYTE numOfController, Real fDegree);

extern QLGate QLAPI CreateCnPh(BYTE numOfController, Real fDegree);

/**
* Burnable Bits: 
*   Guaranteed to be OFF initially, but with no restrictions on state afterwards. 
* 
* ----
* 
* Zeroed Bits: 
*   Guaranteed to be OFF initially, and you must ensure they're OFF when you're done.
* 
*  cnnot ?
*    ?   ?
* 
* Garbage Bits: 
*   Could be in any state initially, and you can add more garbage into the state (you don't have to restore the initial value).
* 
* ----
* 
* Borrowed Bits: 
*   Could be in any state beforehand, and must be restored to that same state afterwards.
* 
* cnnot  0
*   0   cnnot
*/
enum class EAncillaType : UINT
{
    EAT_Burnable,
    EAT_Zeroed,
    EAT_Garbage,
    EAT_Borrowed
};


/**
* total bits is number of controller bits + 1
* The ancilla is the last one qubit
* https://algassert.com/circuits/2015/06/05/Constructing-Large-Controlled-Nots.html
*/
extern QLGate QLAPI CreateCnNotWithAncilla(BYTE numOfController, EAncillaType eAT = EAncillaType::EAT_Zeroed);

/**
* Indeed, there is an ancilla free construction of CnU gate:
* https://algassert.com/circuits/2015/06/22/Using-Quantum-Gates-instead-of-Ancilla-Bits.html
* However, it has a large constant factor, and also suffer from exponential priceition problem, and will cost much time to implement
* So we use an ancilla version
* 
* Only EAT_Zeroed & EAT_Borrowed cases are tested
* 
* cui2 is control-U^-2
*/
extern QLGate QLAPI CreateCnUWithAncillaZeroed(BYTE numOfController, const QLGate& ccu);

extern QLGate QLAPI CreateCnUWithAncilla(BYTE numOfController, const QLGate& cu, const QLGate& cui2, EAncillaType eAT = EAncillaType::EAT_Zeroed);

__END_NAMESPACE


#endif //#ifndef _CNUGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================