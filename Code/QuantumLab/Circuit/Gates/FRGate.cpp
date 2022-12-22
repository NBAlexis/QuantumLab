//=============================================================================
// FILENAME : FRGate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [30/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

TArray<Real> QLAPI SpliteAngles(const TArray<Real>& angles, UINT length)
{
    TArray<Real> retList;
    for (UINT i = 0; i < length; ++i)
    {
        Real theta = F(0.0);
        for (UINT j = 0; j < length; ++j)
        {
            UINT uiSign = BitWiseInnerProduct(j, GrayCode(i));
            theta = theta + ((uiSign & 1) ? (-angles[j]) : (angles[j]));
        }
        retList.AddItem(theta / length);
    }
    return retList;
}

TArray<Real> QLAPI SpliteAngles(const Real* angles, UINT length)
{
    TArray<Real> retList;
    for (UINT i = 0; i < length; ++i)
    {
        Real theta = F(0.0);
        for (UINT j = 0; j < length; ++j)
        {
            UINT uiSign = BitWiseInnerProduct(j, GrayCode(i));
            theta = theta + ((uiSign & 1) ? (-angles[j]) : (angles[j]));
        }
        retList.AddItem(theta / length);
    }
    return retList;
}

QLGate QLAPI FRy(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_RY, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(numberOfQubits);
    ret.m_sName = _T("FRy");

    TArray <Real> theta = SpliteAngles(angles, degreeNumber);
    TArray <BYTE> target;
    target.AddItem(numberOfQubits - 1);
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, theta[i]);
        ret.AppendGate(ry, target);

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        ctrIdx = numberOfQubits - 2 - ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(ctrIdx);
        cnotbits.AddItem(numberOfQubits - 1);
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }

    return ret;
}

QLGate QLAPI FRz(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_RZ, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(numberOfQubits);
    ret.m_sName = _T("FRz");

    TArray <Real> theta = SpliteAngles(angles, degreeNumber);
    TArray <BYTE> target;
    target.AddItem(numberOfQubits - 1);
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate rz(EBasicOperation::EBO_RZ, theta[i]);
        ret.AppendGate(rz, target);

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        ctrIdx = numberOfQubits - 2 - ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(ctrIdx);
        cnotbits.AddItem(numberOfQubits - 1);
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }

    return ret;
}

QLGate QLAPI FRyz(const TArray<Real>& anglesY, const TArray<Real>& anglesZ, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(anglesY.Num()) < degreeNumber || static_cast<UINT>(anglesZ.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        QLGate rzdagger(EBasicOperation::EBO_RZ, -anglesZ[0]);
        //rzdagger.Dagger();
        QLGate ry(EBasicOperation::EBO_RY, anglesY[0]);

        QLGate ret;
        ret.AddQubits(1);
        ret.m_sName = _T("FRyz");
        TArray<BYTE> lstQubits;
        lstQubits.AddItem(0);

        ret.AppendGate(ry, lstQubits);
        ret.AppendGate(rzdagger, lstQubits);
        return ret;
    }

    QLGate ret;
    ret.AddQubits(numberOfQubits);
    ret.m_sName = _T("FRyz");

    TArray <Real> thetaY = SpliteAngles(anglesY, degreeNumber);
    TArray <Real> thetaZ = SpliteAngles(anglesZ, degreeNumber);

    TArray <BYTE> target;
    target.AddItem(numberOfQubits - 1);
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, thetaY[i]);
        ret.AppendGate(ry, target);

        if (i != degreeNumber - 1)
        {
            UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;
            TArray <BYTE> cnotbits;
            cnotbits.AddItem(ctrIdx);
            cnotbits.AddItem(numberOfQubits - 1);
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }
    }

    for (UINT i = 0; i < degreeNumber; ++i)
    {
        UINT j = degreeNumber - 1 - i;
        if (0 != i)
        {
            UINT ctrIdx = GrayCodeDifferent(j, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;

            TArray <BYTE> cnotbits;
            cnotbits.AddItem(ctrIdx);
            cnotbits.AddItem(numberOfQubits - 1);
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }

        QLGate rz(EBasicOperation::EBO_RZ, -thetaZ[j]);
        ret.AppendGate(rz, target);
    }

    return ret;
}

QLGate QLAPI FRyz(const Real* anglesY, const Real* anglesZ, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);

    if (1 == numberOfQubits)
    {
        QLGate rzdagger(EBasicOperation::EBO_RZ, -anglesZ[0]);
        //rzdagger.Dagger();
        QLGate ry(EBasicOperation::EBO_RY, anglesY[0]);

        QLGate ret;
        ret.AddQubits(1);
        ret.m_sName = _T("FRyz");
        TArray<BYTE> lstQubits;
        lstQubits.AddItem(0);

        ret.AppendGate(ry, lstQubits);
        ret.AppendGate(rzdagger, lstQubits);
        return ret;
    }

    QLGate ret;
    ret.AddQubits(numberOfQubits);
    ret.m_sName = _T("FRyz");

    TArray <Real> thetaY = SpliteAngles(anglesY, degreeNumber);
    TArray <Real> thetaZ = SpliteAngles(anglesZ, degreeNumber);

    TArray <BYTE> target;
    target.AddItem(numberOfQubits - 1);
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, thetaY[i]);
        ret.AppendGate(ry, target);

        if (i != degreeNumber - 1)
        {
            UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;
            TArray <BYTE> cnotbits;
            cnotbits.AddItem(ctrIdx);
            cnotbits.AddItem(numberOfQubits - 1);
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }
    }

    for (UINT i = 0; i < degreeNumber; ++i)
    {
        UINT j = degreeNumber - 1 - i;
        if (0 != i)
        {
            UINT ctrIdx = GrayCodeDifferent(j, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;

            TArray <BYTE> cnotbits;
            cnotbits.AddItem(ctrIdx);
            cnotbits.AddItem(numberOfQubits - 1);
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }

        QLGate rz(EBasicOperation::EBO_RZ, -thetaZ[j]);
        ret.AppendGate(rz, target);
    }

    return ret;
}

QLGate QLAPI FRp(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_P, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(numberOfQubits);
    ret.m_sName = _T("FRp");

    QLGate x(EBasicOperation::EBO_X);
    TArray<BYTE> allBits;
    for (UINT i = 0; i < numberOfQubits - 1; ++i)
    {
        TArray<BYTE> targetbit;
        targetbit.AddItem(i);
        ret.AppendGate(x, targetbit);

        allBits.AddItem(i);
    }
    allBits.AddItem(numberOfQubits - 1);

    QLGate cnp = CreateCnP(numberOfQubits - 1, angles[0]);
    ret.AppendGate(cnp, allBits);

    for (UINT i = 1; i < degreeNumber; ++i)
    {
        UINT diff = i ^ (i - 1);
        for (UINT j = 0; j < numberOfQubits - 1; ++j)
        {
            if (diff & (1 << j))
            {
                TArray<BYTE> cnotbit;
                cnotbit.AddItem(j);
                ret.AppendGate(x, cnotbit);
            }
        }

        QLGate cnp2 = CreateCnP(numberOfQubits - 1, angles[i]);
        ret.AppendGate(cnp2, allBits);
    }
    return ret;
}

QLGate QLAPI FRPh(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_Phase, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(numberOfQubits);
    ret.m_sName = _T("FRPh");

    QLGate x(EBasicOperation::EBO_X);
    TArray<BYTE> allBits;
    for (UINT i = 0; i < numberOfQubits - 1; ++i)
    {
        TArray<BYTE> targetbit;
        targetbit.AddItem(i);
        ret.AppendGate(x, targetbit);

        allBits.AddItem(i);
    }
    allBits.AddItem(numberOfQubits - 1);

    QLGate cnp = CreateCnPh(numberOfQubits - 1, angles[0]);
    ret.AppendGate(cnp, allBits);

    for (UINT i = 1; i < degreeNumber; ++i)
    {
        UINT diff = i ^ (i - 1);
        for (UINT j = 0; j < numberOfQubits - 1; ++j)
        {
            if (diff & (1 << j))
            {
                TArray<BYTE> cnotbit;
                cnotbit.AddItem(j);
                ret.AppendGate(x, cnotbit);
            }
        }

        QLGate cnp2 = CreateCnPh(numberOfQubits - 1, angles[i]);
        ret.AppendGate(cnp2, allBits);
    }
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================