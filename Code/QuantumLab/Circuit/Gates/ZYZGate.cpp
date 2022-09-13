//=============================================================================
// FILENAME : ZYZGate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

std::vector<Real> GetZYZDecompose(const QLMatrix& u, UBOOL bNormalize)
{
    std::vector<Real> ret;

    Real theta = F(0.0);
    Real alpha = F(0.0);
    Real delta = F(0.0);
    Real beta = F(0.0);

    QLComplex u11 = u.Get(0, 0);
    QLComplex u12 = u.Get(0, 1);
    QLComplex u21 = u.Get(1, 0);
    QLComplex u22 = u.Get(1, 1);

    QLComplex m12 = _make_cuComplex(-u12.x, -u12.y);

    if (_cuCabsf(u11) < REAL_EPS)
    {
        theta = PI;
        delta = (__cuCargf(m12) + __cuCargf(u21)) / F(2.0);
        beta = __cuCargf(m12) - __cuCargf(u21);
        ret.push_back(theta);
        ret.push_back(beta);
        ret.push_back(alpha);
        ret.push_back(delta);
        return ret;
    }

    if (_cuCabsf(u12) < REAL_EPS)
    {
        delta = (__cuCargf(u22) + __cuCargf(u11)) / F(2.0);
        beta = __cuCargf(u22) - __cuCargf(u11);
        ret.push_back(theta);
        ret.push_back(beta);
        ret.push_back(alpha);
        ret.push_back(delta);
        return ret;
    }

    if (bNormalize)
    {
        Real fNormA = __cuCabsSqf(u11) + __cuCabsSqf(u12);
        QLComplex dotAB = _cuCaddf(_cuCmulf(_cuConjf(u11), u21), _cuCmulf(_cuConjf(u12), u22));
        u21 = _cuCsubf(u21, cuCdivf_cr(_cuCmulf(u11, dotAB), fNormA));
        u22 = _cuCsubf(u22, cuCdivf_cr(_cuCmulf(u12, dotAB), fNormA));
        Real fNormC = _sqrt(__cuCabsSqf(u21) + __cuCabsSqf(u22));
        fNormA = _sqrt(fNormA);

        u11 = cuCdivf_cr(u11, fNormA);
        u12 = cuCdivf_cr(u12, fNormA);
        u21 = cuCdivf_cr(u21, fNormA);
        u22 = cuCdivf_cr(u22, fNormA);
    }

    theta = F(2.0) * _acos(_cuCabsf(u11));
    Real cs = _cos(theta / 2);
    Real sn = _sin(theta / 2);

    //this can be simplified since the arg only changes when cs, sn < 0
    u11 = cuCdivf_cr(u11, cs);
    u12 = cuCdivf_cr(u12, sn);
    u22 = cuCdivf_cr(u22, cs);
    beta = __cuCargf(u11) - __cuCargf(u12);
    alpha = __cuCargf(u12) - __cuCargf(u22);
    delta = (__cuCargf(u11) - alpha / 2 - beta / 2);

    ret.push_back(theta);
    ret.push_back(beta);
    ret.push_back(alpha);
    ret.push_back(delta);
    return ret;
}

extern QLGate QLAPI CreateZYZGate(const QLMatrix& u, UBOOL bNormalize)
{
    QLGate retGate;
    retGate.m_bBasicOperation = FALSE;
    retGate.m_lstQubits.AddItem(0);
    
    std::vector<Real> degrees = GetZYZDecompose(u, bNormalize);

    QLGate rz1(EBasicOperation::EBO_RZ, -degrees[1]);
    QLGate ry(EBasicOperation::EBO_RY, degrees[0]);
    QLGate rz2(EBasicOperation::EBO_RZ, -degrees[2]);
    QLGate ph(EBasicOperation::EBO_Phase, degrees[3]);

    retGate.AppendGate(rz1, retGate.m_lstQubits);
    retGate.AppendGate(ry, retGate.m_lstQubits);
    retGate.AppendGate(rz2, retGate.m_lstQubits);
    retGate.AppendGate(ph, retGate.m_lstQubits);

    return retGate;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================