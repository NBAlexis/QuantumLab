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
    const QLComplex s = __cuCsqrtf(delta);
    const QLComplex t = __cuCsqrtf(_cuCaddf(tau, cuCmulf_cr(s, F(2.0))));

    QLComplex v[] = { _cuCdivf(_cuCaddf(u00, s), t), _cuCdivf(u01, t),  _cuCdivf(u10, t),  _cuCdivf(_cuCaddf(u11, s), t) };
    return QLMatrix::CopyCreate(2, 2, v);
}

/**
* total bits is number of controller bits + 1
* The target is the last one qubit
*/
extern QLGate QLAPI CreateCnNot(BYTE controllerBits);

inline QLGate CreateToffoliGate()
{
    QLGate ret = CreateCnNot(2);
    ret.m_sName = _T("CCN");
    return ret;
}

extern QLGate QLAPI CreateCnU(BYTE controllerBits, const QLMatrix& mtr);

__END_NAMESPACE


#endif //#ifndef _CNUGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================