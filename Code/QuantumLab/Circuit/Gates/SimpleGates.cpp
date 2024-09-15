//=============================================================================
// FILENAME : SimpleGates.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

/**
* U = [u11 u12]
*     [u21 u22]
*
* U = exp(i delta) [ exp(i (alpha+beta)/2) cos(theta/2)  exp(i (alpha-beta)/2) sin(theta/2)] = exp(i delta) Rz(-alpha) Ry(-theta) Rz(-beta)
*                  [-exp(-i(alpha-beta)/2) sin(theta/2)  exp(-i(alpha+beta)/2) cos(theta/2)]
*
* for OpenQASM
*
*
*
* U(a,b,c)=Rz(b)Ry(a)Rz(c)|psi>, where a=-theta, b=-alpha, c=-beta
* a,b,c = degree(0), degree(2), degree(1)
*/
TArray<Real> GetZYZDecompose(const QLMatrix& u)
{
    TArray<Real> ret;

    const QLComplex u11 = u.Get(0, 0);
    const Real amp11 = _cuCabsf(u11);
    const QLComplex u12 = u.Get(1, 0);
    const Real amp12 = _cuCabsf(u12);
    

    Real beta = F(0.0);
    Real alpha = F(0.0);
    Real delta = F(0.0);

    if (amp11 < _QL_FLT_MIN_)
    {
        const QLComplex u21 = u.Get(0, 1);
        const Real a12 = __cuCargf(u12);
        const Real a21 = __cuCargf(u21);
        alpha = (a12 - a21 + PI) / 2;
        beta = -alpha;
        delta = (a12 + a21 - PI) / 2;
    }
    else
    {
        const QLComplex u22 = u.Get(1, 1);
        const Real a12 = __cuCargf(u12);
        const Real a11 = __cuCargf(u11);
        const Real a22 = __cuCargf(u22);
        beta = a11 - a12;
        alpha = a12 - a22;
        delta = a11 - alpha / 2 - beta / 2;
    }

    ret.AddItem(F(2.0) * _atan2(amp12, amp11)); //theta
    ret.AddItem(beta);
    ret.AddItem(alpha);
    ret.AddItem(delta); 

    //appGeneral(_T("==============\n"));

    //u.Print("u");

    //QLMatrix r = QLGate::CreateSingleQubitMatrix(EBasicOperation::EBO_RZ, -ret[1]);
    //r = QLGate::CreateSingleQubitMatrix(EBasicOperation::EBO_RY, -ret[0]) * r;
    //r = QLGate::CreateSingleQubitMatrix(EBasicOperation::EBO_RZ, -ret[2]) * r;
    //r = r * _make_cuComplex(cos(ret[3]), sin(ret[3]));

    //appGeneral(_T("delta = %f\n"), ret[3]);

    //r.Print("r");

    //appGeneral(_T("==============\n"));

    return ret;
}

QLGate QLAPI CreateZYZGate(const QLMatrix& u)
{
    QLGate retGate;
    retGate.m_lstQubits.AddItem(0);
    retGate.m_sName = _T("ZYZ");
    
    TArray<Real> degrees = GetZYZDecompose(u);

    QLGate rz1(EBasicOperation::EBO_RZ, -degrees[1]);
    QLGate ry(EBasicOperation::EBO_RY, -degrees[0]);
    QLGate rz2(EBasicOperation::EBO_RZ, -degrees[2]);
    QLGate ph(EBasicOperation::EBO_Phase, degrees[3]);

    retGate.AppendGate(rz1, retGate.m_lstQubits);
    retGate.AppendGate(ry, retGate.m_lstQubits);
    retGate.AppendGate(rz2, retGate.m_lstQubits);
    retGate.AppendGate(ph, retGate.m_lstQubits);

    return retGate;
}

QLGate QLAPI CreateControlledZYZGate(const QLMatrix& u)
{
    QLGate retGate;
    retGate.m_lstQubits.AddItem(0);
    retGate.m_lstQubits.AddItem(1);
    retGate.m_sName = _T("CZYZ");

    TArray<BYTE> controller;
    TArray<BYTE> target;
    controller.AddItem(0);
    target.AddItem(1);

    TArray<Real> degrees = GetZYZDecompose(u);

    QLGate rz1(EBasicOperation::EBO_RZ, (degrees[2] - degrees[1]) / 2);
    QLGate cnot1(EBasicOperation::EBO_CX);
    QLGate rz2(EBasicOperation::EBO_RZ, (degrees[2] + degrees[1]) / 2);
    QLGate ry1(EBasicOperation::EBO_RY, degrees[0] / 2);
    QLGate cnot2(EBasicOperation::EBO_CX);
    QLGate ry2(EBasicOperation::EBO_RY, -degrees[0] / 2);
    QLGate rz3(EBasicOperation::EBO_RZ, -degrees[2]);
    QLGate p(EBasicOperation::EBO_P, degrees[3]);

    retGate.AppendGate(rz1, target);
    retGate.AppendGate(cnot1, retGate.m_lstQubits);
    retGate.AppendGate(rz2, target);
    retGate.AppendGate(ry1, target);
    retGate.AppendGate(cnot2, retGate.m_lstQubits);
    retGate.AppendGate(ry2, target);
    retGate.AppendGate(rz3, target);
    retGate.AppendGate(p, controller);

    return retGate;
}

QLGate QLAPI CreateSwapGate()
{
    QLGate retGate;
    retGate.m_lstQubits.AddItem(0);
    retGate.m_lstQubits.AddItem(1);
    retGate.m_sName = _T("Swap");

    TArray<BYTE> inverseQubit;
    inverseQubit.AddItem(1);
    inverseQubit.AddItem(0);

    QLGate cnot(EBasicOperation::EBO_CX);
    retGate.AppendGate(cnot, retGate.m_lstQubits);
    retGate.AppendGate(cnot, inverseQubit);
    retGate.AppendGate(cnot, retGate.m_lstQubits);

    return retGate;
}

QLGate QLAPI CreateControlledSwap(BYTE controller)
{
    QLGate cnot(EBasicOperation::EBO_CX);
    QLGate cnnot = CreateCnNot(controller + 1);

    TArray<BYTE> qubits;
    qubits.Append(ByteSequnce, controller + 2);

    QLGate ret;
    ret.AddQubits(controller + 2);
    ret.m_sName = _T("cswap");
    ret.AppendGate(cnot, controller + 1, controller);
    ret.AppendGate(cnnot, qubits);
    ret.AppendGate(cnot, controller + 1, controller);
    return ret;
}

QLGate QLAPI CreateControlledHadamardGate()
{
    QLGate retGate;
    retGate.m_lstQubits.AddItem(0);
    retGate.m_lstQubits.AddItem(1);
    retGate.m_sName = _T("CH");

    TArray<BYTE> target;
    target.AddItem(1);

    QLGate cy1(EBasicOperation::EBO_RY, PI / 4);
    QLGate cnot(EBasicOperation::EBO_CX);
    QLGate cy2(EBasicOperation::EBO_RY, -PI / 4);

    retGate.AppendGate(cy1, target);
    retGate.AppendGate(cnot, retGate.m_lstQubits);
    retGate.AppendGate(cy2, target);

    return retGate;
}
__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================