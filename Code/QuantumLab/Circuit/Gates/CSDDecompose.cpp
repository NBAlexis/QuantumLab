//=============================================================================
// FILENAME : CSDDecompose.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [30/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

void CSDMatrix::SetAsUMatrix(INT level, INT totLevel, const TArray<QLMatrix>& u)
{
    UINT halfSize = 1U << (level - 1);
    m_iLevel = level;
    m_iTotLevel = totLevel;

    if (level == totLevel)
    {
        m_eType = ECSDMatrix::DiagonalUMatrix;
        m_lstDegreeContent.RemoveAll();
        for (UINT i = 0; i < halfSize; ++i)
        {
            //m_lstDegreeContent.AddItem(__cuCargf(u[i].Get(0, 0)));
            //m_lstDegreeContent.AddItem(__cuCargf(u[i].Get(1, 1)));

            m_lstDegreeContent.AddItem(__cuCargf(u[2 * i].Get(0, 0)));
            m_lstDegreeContent.AddItem(__cuCargf(u[2 * i + 1].Get(0, 0)));
        }
        return;
    }

    m_eType = ECSDMatrix::BlockedUMatrix;
    UINT halfMatrixSize = 1U << (totLevel - level);
    UINT totalMatrixSize = halfMatrixSize << 1;
    m_lstSubMatrix.RemoveAll();
    for (UINT i = 0; i < halfSize; ++i)
    {
        //m_lstSubMatrix.AddItem(u[i].GetBlock(0, halfMatrixSize, 0 , halfMatrixSize));
        //m_lstSubMatrix.AddItem(u[i].GetBlock(halfMatrixSize, halfMatrixSize, halfMatrixSize, halfMatrixSize));
        m_lstSubMatrix.AddItem(u[2 * i]);
        m_lstSubMatrix.AddItem(u[2 * i + 1]);
    }
}

void CSDMatrix::SetAsAMatrix(INT level, INT totLevel, const TArray<QLMatrix>& a)
{
    UINT matrixCount = 1U << (level - 1);
    m_iLevel = level;
    m_iTotLevel = totLevel;
    m_eType = ECSDMatrix::AMatrix;
    UINT halfMatrixSize = 1 << (totLevel - level);
    m_lstDegreeContent.RemoveAll();
    for (UINT i = 0; i < matrixCount; ++i)
    {
        for (UINT j = 0; j < halfMatrixSize; ++j)
        {
            Real fReal = a[i].Get(j, j).x;
            if (fReal > F(1.0))
            {
                fReal = F(1.0);
            }
            else if (fReal < F(-1.0))
            {
                fReal = F(-1.0);
            }
            m_lstDegreeContent.AddItem(_acos(fReal));
        }
    }
}

void CSDMatrix::SetAsRzMatrix(INT level, INT totLevel, const TArray<Real>& degrees)
{
    appCrucial(_T("should never call me!\n"));

    UINT degreeCount = 1U << (totLevel - 1);
    if (degreeCount != static_cast<UINT>(degrees.Num()))
    {
        appCrucial(_T("SetAsRzMatrix: something bad happens!\n"));
        return;
    }

    m_iLevel = level;
    m_iTotLevel = totLevel;
    m_eType = ECSDMatrix::FRZMatrix;
    m_lstDegreeContent = degrees;
}

void CSDMatrix::SetAsGlobalPhase(Real fPhase, INT totLevel)
{
    m_eType = ECSDMatrix::GlobalPhase;
    m_iTotLevel = totLevel;
    m_lstDegreeContent.RemoveAll();
    m_lstDegreeContent.AddItem(fPhase);
}

void CSDMatrix::SetAsPartialRZ(INT level, INT totLevel, const TArray<Real>& degrees)
{
    m_eType = ECSDMatrix::PartialFRZMatrix;
    m_iLevel = level;
    m_iTotLevel = totLevel;
    m_lstDegreeContent = degrees;
}

void CSDMatrix::SetAsDiagonalUMatrix(INT level, INT totLevel, const TArray<Real>& degrees)
{
    UINT degreeCount = 1U << totLevel;
    if (degreeCount != static_cast<UINT>(degrees.Num()))
    {
        appCrucial("something bad happens! SetAsDiagonalUMatrix\n");
        return;
    }
    m_iLevel = level;
    m_iTotLevel = totLevel;
    m_eType = ECSDMatrix::DiagonalUMatrix;
    m_lstDegreeContent = degrees;
}

QLMatrix CSDMatrix::PrintAsAMatrix() const
{
    UINT totalSize = 1U << m_iTotLevel;
    QLMatrix toPrint(totalSize, totalSize);
    UINT matrixCount = 1U << (m_iLevel - 1);
    UINT halfMatrixSize = 1U << (m_iTotLevel - m_iLevel);
    UINT matrixSize = halfMatrixSize << 1;
    QLComplex* oneMatrix = (QLComplex*)malloc(sizeof(QLComplex) * matrixSize * matrixSize);
    if (NULL == oneMatrix)
    {
        return QLMatrix();
    }
    for (UINT i = 0; i < matrixCount; ++i)
    {
        memset(oneMatrix, 0, sizeof(QLComplex) * matrixSize * matrixSize);
        for (UINT j = 0; j < halfMatrixSize; ++j)
        {
            Real degree = m_lstDegreeContent[i * halfMatrixSize + j];
            //matrix x * ylength + y
            oneMatrix[j * matrixSize + j] = _make_cuComplex(_cos(degree), F(0.0));
            oneMatrix[j * matrixSize + j + halfMatrixSize] = _make_cuComplex(_sin(degree), F(0.0));
            oneMatrix[(j + halfMatrixSize) * matrixSize + j] = _make_cuComplex(-_sin(degree), F(0.0));
            oneMatrix[(j + halfMatrixSize) * matrixSize + j + halfMatrixSize] = _make_cuComplex(_cos(degree), F(0.0));
        }
        toPrint.SetBlock(i * matrixSize, matrixSize, i * matrixSize, matrixSize, oneMatrix);
    }
    free(oneMatrix);
    return toPrint;
}

QLMatrix CSDMatrix::PrintAsBlockedUMatrix() const
{
    UINT totalSize = 1U << m_iTotLevel;
    QLMatrix toPrint(totalSize, totalSize);
    UINT matrixCount = 1U << m_iLevel;
    UINT matrixSize = 1U << (m_iTotLevel - m_iLevel);
    for (UINT i = 0; i < matrixCount; ++i)
    {
        toPrint.SetBlock(i * matrixSize, matrixSize, i * matrixSize, matrixSize, m_lstSubMatrix[i].HostBuffer());
    }
    return toPrint;
}

QLMatrix CSDMatrix::PrintAsDiagonalUMatrix() const
{
    UINT totalSize = 1U << m_iTotLevel;
    QLComplex* oneMatrix = (QLComplex*)malloc(sizeof(QLComplex) * totalSize * totalSize);
    if (NULL == oneMatrix)
    {
        return QLMatrix();
    }
    memset(oneMatrix, 0, sizeof(QLComplex) * totalSize * totalSize);
    for (UINT i = 0; i < totalSize; ++i)
    {
        oneMatrix[i * totalSize + i] = _make_cuComplex(_cos(m_lstDegreeContent[i]), _sin(m_lstDegreeContent[i]));
    }
    return QLMatrix(totalSize, totalSize, oneMatrix);
}

QLMatrix CSDMatrix::PrintAsFRZMatrix() const
{
    UINT totalSize = 1U << m_iTotLevel;
    QLComplex* oneMatrix = (QLComplex*)malloc(sizeof(QLComplex) * totalSize * totalSize);
    if (NULL == oneMatrix)
    {
        return QLMatrix();
    }
    memset(oneMatrix, 0, sizeof(QLComplex) * totalSize * totalSize);
    UINT matrixCount = 1U << (m_iLevel - 1);
    UINT halfMatrixSize = 1U << (m_iTotLevel - m_iLevel);
    UINT matrixSize = halfMatrixSize << 1;
    for (UINT i = 0; i < matrixCount; ++i)
    {
        for (UINT j = 0; j < halfMatrixSize; ++j)
        {
            Real degree = m_lstDegreeContent[i * halfMatrixSize + j];
            UINT position = i * matrixSize + j;
            oneMatrix[position * totalSize + position] = _make_cuComplex(_cos(degree), _sin(degree));
            oneMatrix[(position + halfMatrixSize) * totalSize + position + halfMatrixSize] = _make_cuComplex(_cos(degree),- _sin(degree));
        }
    }
    return QLMatrix(totalSize, totalSize, oneMatrix);
}

QLMatrix CSDMatrix::PrintAsPartialFRZMatrix() const
{
    UINT totalSize = 1U << m_iTotLevel;
    QLComplex* oneMatrix = (QLComplex*)malloc(sizeof(QLComplex) * totalSize * totalSize);
    if (NULL == oneMatrix)
    {
        return QLMatrix();
    }
    memset(oneMatrix, 0, sizeof(QLComplex) * totalSize * totalSize);
    UINT degreeCount = 1U << m_iLevel;
    UINT repeatLength = 1U << (m_iTotLevel - m_iLevel - 1);
    for (UINT i = 0; i < degreeCount; ++i)
    {
        for (UINT j = 0; j < repeatLength; ++j)
        {
            oneMatrix[(2 * i * repeatLength + j) * totalSize + 2 * i * repeatLength + j] = _make_cuComplex(_cos(m_lstDegreeContent[i]), _sin(m_lstDegreeContent[i]));
            oneMatrix[((2 * i + 1) * repeatLength + j) * totalSize + (2 * i + 1) * repeatLength + j] = _make_cuComplex(_cos(m_lstDegreeContent[i]), -_sin(m_lstDegreeContent[i]));
        }
    }
    return QLMatrix(totalSize, totalSize, oneMatrix);
}

QLMatrix CSDMatrix::PrintAsGlobalPhase() const
{
    UINT totalSize = 1U << m_iTotLevel;
    QLComplex* oneMatrix = (QLComplex*)malloc(sizeof(QLComplex) * totalSize * totalSize);
    if (NULL == oneMatrix)
    {
        return QLMatrix();
    }
    memset(oneMatrix, 0, sizeof(QLComplex) * totalSize * totalSize);
    for (UINT i = 0; i < totalSize; ++i)
    {
        oneMatrix[i * totalSize + i] = _make_cuComplex(_cos(m_lstDegreeContent[0]), _sin(m_lstDegreeContent[0]));
    }
    return QLMatrix(totalSize, totalSize, oneMatrix);
}

void CSDMatrix::GateAsAMatrix(QLGate& gate) const
{
    TArray <BYTE> controller;
    BYTE targetBit = 0;
    for (INT i = 0; i < m_iTotLevel; ++i)
    {
        BYTE consideringBit = static_cast<BYTE>(m_iTotLevel - i - 1);
        if (consideringBit == static_cast<BYTE>(m_iTotLevel - m_iLevel))
        {
            targetBit = consideringBit;
        }
        else
        {
            controller.AddItem(consideringBit);
        }
    }
    controller.AddItem(targetBit);

    TArray <Real> degrees;
    for (INT i = 0; i < m_lstDegreeContent.Num(); ++i)
    {
        degrees.AddItem(F(2.0) * m_lstDegreeContent[i]);
    }

    QLGate fry = FRy(degrees, static_cast<BYTE>(controller.Num()));
    fry.m_sName = _T("cs");
    gate.AppendGate(fry, controller);
}

void CSDMatrix::GateAsFRZMatrix(QLGate& gate) const
{
    TArray <BYTE> controller;
    BYTE targetBit = 0;
    for (INT i = 0; i < m_iTotLevel; ++i)
    {
        BYTE consideringBit = static_cast<BYTE>(m_iTotLevel - i - 1);
        if (consideringBit == static_cast<BYTE>(m_iTotLevel - m_iLevel))
        {
            targetBit = consideringBit;
        }
        else
        {
            controller.AddItem(consideringBit);
        }
    }
    controller.AddItem(targetBit);

    TArray <Real> degrees;
    for (INT i = 0; i < m_lstDegreeContent.Num(); ++i)
    {
        degrees.AddItem(-F(2.0) * m_lstDegreeContent[i]);
    }

    QLGate frz = FRz(degrees, static_cast<BYTE>(controller.Num()));
    frz.m_sName = _T("u");
    gate.AppendGate(frz, controller);
}

void CSDMatrix::GateAsPartialFRZMatrix(QLGate& gate) const
{
    if (0 == m_iLevel)
    {
        QLGate rz(EBasicOperation::EBO_RZ, -m_lstDegreeContent[0] * F(2.0));
        TArray <BYTE> bit;
        bit.AddItem(m_iTotLevel - 1);
        gate.AppendGate(rz, bit);
    }
    else
    {
        TArray <BYTE> allBits;
        for (INT i = 0; i < m_iLevel; ++i)
        {
            allBits.AddItem(static_cast<BYTE>(m_iTotLevel - 1 - i));
        }
        allBits.AddItem(static_cast<BYTE>(m_iTotLevel - 1 - m_iLevel));

        TArray <Real> degree;
        for (INT i = 0; i < m_lstDegreeContent.Num(); ++i)
        {
            degree.AddItem(-F(2.0) * m_lstDegreeContent[i]);
        }

        QLGate frz = FRz(degree, static_cast<BYTE>(allBits.Num()));
        gate.AppendGate(frz, allBits);
    }
}

void CSDMatrix::GateAsGlobalPhase(QLGate& gate) const
{
    QLGate ph(EBasicOperation::EBO_Phase, m_lstDegreeContent[0]);
    TArray <BYTE> bit0;
    bit0.AddItem(0);
    gate.AppendGate(ph, bit0);
}

void CSDMatrix::CombinedYZGate(QLGate& gate, const TArray<Real>& zDegrees) const
{
    if (ECSDMatrix::AMatrix != m_eType)
    {
        appCrucial(_T("something very bad happens CombinedYZGate\n"));
        return;
    }

    TArray<BYTE> allBits;
    BYTE target = 0;
    for (INT i = 0; i < m_iTotLevel; ++i)
    {
        BYTE consideringBit = static_cast<BYTE>(m_iTotLevel - i - 1);
        if (consideringBit == static_cast<BYTE>(m_iTotLevel - m_iLevel))
        {
            target = consideringBit;
        }
        else
        {
            allBits.AddItem(consideringBit);
        }
    }
    allBits.AddItem(target);
    
    TArray <Real> degreeY;
    TArray <Real> degreeZ;
    for (INT i = 0; i < m_lstDegreeContent.Num(); ++i)
    {
        degreeY.AddItem(F(2.0) * m_lstDegreeContent[i]);
        degreeZ.AddItem(zDegrees[i] * F(2.0));
    }
    QLGate fryz = FRyz(degreeY, degreeZ, static_cast<BYTE>(allBits.Num()));

    gate.AppendGate(fryz, allBits);
}

TArray<CSDMatrix> CSDMatrix::DecomposeU() const
{
    TArray<CSDMatrix> ret;
    if (ECSDMatrix::BlockedUMatrix != m_eType)
    {
        appCrucial("something bad happens DecomposeU!\n");
        return ret;
    }
    
    INT matrixCount = m_lstSubMatrix.Num();
    TArray <QLMatrix> l;
    TArray <QLMatrix> a;
    TArray <QLMatrix> r;
    UINT halfSize = 1U << (m_iTotLevel - m_iLevel - 1);
    for (INT i = 0; i < matrixCount; ++i)
    {
        QLMatrix q1, q2, c, s, v1, v2;
        m_lstSubMatrix[i].CSD(q1, q2, c, s, v1, v2, halfSize, halfSize);

        //QLMatrix q, cs, v;
        //QLMatrix::CombinedCSD(q, cs, v, q1, q2, c, s, v1, v2, halfSize, halfSize);

        //l.AddItem(q);
        //a.AddItem(cs);
        //v.Dagger();
        //r.AddItem(v);

        v1.Dagger();
        v2.Dagger();

        l.AddItem(q1);
        l.AddItem(q2);
        a.AddItem(c);
        r.AddItem(v1);
        r.AddItem(v2);
    }

    CSDMatrix lmatrix;
    lmatrix.SetAsUMatrix(m_iLevel + 1, m_iTotLevel, l);
    CSDMatrix amatrix;
    amatrix.SetAsAMatrix(m_iLevel + 1, m_iTotLevel, a);
    CSDMatrix rmatrix;
    rmatrix.SetAsUMatrix(m_iLevel + 1, m_iTotLevel, r);

    ret.AddItem(lmatrix);
    ret.AddItem(amatrix);
    ret.AddItem(rmatrix);
    return ret;
}

TArray<Real> CSDMatrix::CalculateUDegrees(INT iLevel)
{
    TArray<Real> degreesToShift;
    if (ECSDMatrix::DiagonalUMatrix != m_eType)
    {
        appCrucial(_T("something bad happens CalculateUDegrees!\n"));
        return degreesToShift;
    }
    TArray<Real> degreesResult;
    UINT matrixCount = 1U << (iLevel - 1);
    UINT halfMatrixSize = 1U << (m_iTotLevel - iLevel);
    UINT matrixSize = halfMatrixSize << 1;

    for (UINT i = 0; i < matrixCount; ++i)
    {
        for (UINT j = 0; j < halfMatrixSize; ++j)
        {
            UINT position = i * matrixSize + j;
            Real degree1 = m_lstDegreeContent[position];
            Real degree2 = m_lstDegreeContent[position + halfMatrixSize];
            Real toShift = (degree1 + degree2) / 2;
            degreesToShift.AddItem(toShift);
            degreesResult.AddItem(degree1 - toShift);
        }
    }

    m_eType = ECSDMatrix::FRZMatrix;
    m_iLevel = iLevel;
    m_lstDegreeContent = degreesResult;
    return degreesToShift;
}

void CSDMatrix::PhaseShiftU(const TArray<Real>& toShift, INT iLevel)
{
    if (ECSDMatrix::DiagonalUMatrix != m_eType)
    {
        appCrucial(_T("something bad happens PhaseShiftU!\n"));
        return;
    }
    UINT matrixCount = 1U << (iLevel - 1);
    UINT halfMatrixSize = 1U << (m_iTotLevel - iLevel);
    UINT matrixSize = halfMatrixSize << 1;

    for (UINT i = 0; i < matrixCount; ++i)
    {
        for (UINT j = 0; j < halfMatrixSize; ++j)
        {
            UINT position = i * matrixSize + j;
            Real degree = toShift[i * halfMatrixSize + j];
            m_lstDegreeContent[position] = m_lstDegreeContent[position] + degree;
            m_lstDegreeContent[position + halfMatrixSize] = m_lstDegreeContent[position + halfMatrixSize] + degree;
        }
    }
}

TArray<CSDMatrix> CSDMatrix::DecomposeUIntoRZ() const
{
    TArray<CSDMatrix> allDegrees;
    if (ECSDMatrix::DiagonalUMatrix != m_eType)
    {
        appCrucial("something bad happens DecomposeU!\n");
        return allDegrees;
    }

    Real a0 = F(0.0);
    for (INT i = 0; i < m_lstDegreeContent.Num(); ++i)
    {
        a0 += m_lstDegreeContent[i];
    }
    a0 = a0 / (1U << m_iTotLevel);
    CSDMatrix a0Matrix;
    a0Matrix.SetAsGlobalPhase(a0, m_iTotLevel);
    allDegrees.AddItem(a0Matrix);

    for (INT i = 0; i < m_iTotLevel; ++i)
    {
        UINT degreeCount = 1U << i;
        Real factor = F(1.0) / (1U << (m_iTotLevel - i));
        UINT sumLength = 1U << (m_iTotLevel - i - 1);
        TArray<Real> degreeList;
        for (UINT j = 0; j < degreeCount; ++j)
        {
            Real plus = F(0.0);
            for (UINT k = 2 * j * sumLength; k < (2 * j + 1) * sumLength; ++k)
            {
                plus = plus + m_lstDegreeContent[k];
            }
            Real minus = F(0.0);
            for (UINT k = (2 * j + 1) * sumLength; k < (2 * j + 2) * sumLength; ++k)
            {
                minus = minus + m_lstDegreeContent[k];
            }
            degreeList.AddItem(factor * (plus - minus));
        }

        CSDMatrix partialZMatrix;
        partialZMatrix.SetAsPartialRZ(i, m_iTotLevel, degreeList);
        allDegrees.AddItem(partialZMatrix);
    }
    return allDegrees;
}

TArray<CSDMatrix> CSDMatrix::RecursiveUDecompose(const TArray<CSDMatrix>& mtrList)
{
    UBOOL toContinue = TRUE;
    TArray<CSDMatrix> retlist = mtrList;
    while (toContinue)
    {
        TArray<CSDMatrix> newList;
        toContinue = FALSE;
        for (INT i = 0; i < retlist.Num(); ++i)
        {
            if (ECSDMatrix::BlockedUMatrix == retlist[i].m_eType)
            {
                TArray<CSDMatrix> csd = retlist[i].DecomposeU();
                newList.Append(csd);
                toContinue = TRUE;
            }
            else
            {
                newList.AddItem(retlist[i]);
            }
        }
        retlist = newList;
    }
    return retlist;
}

void CSDMatrix::PhaseShiftUAll(TArray<CSDMatrix>& csdMatrix)
{
    for (INT i = 0; i < (csdMatrix.Num() - 1) / 2; ++i)
    {
        if (   ECSDMatrix::DiagonalUMatrix != csdMatrix[2 * i].m_eType
            || ECSDMatrix::AMatrix         != csdMatrix[2 * i + 1].m_eType
            || ECSDMatrix::DiagonalUMatrix != csdMatrix[2 * i + 2].m_eType)
        {
            appCrucial("something bad happens PhaseShiftUAll! %d %d %d %d\n", 2 * i, 
                csdMatrix[2 * i].m_eType,
                csdMatrix[2 * i + 1].m_eType,
                csdMatrix[2 * i + 2].m_eType);
            return;
        }

        TArray<Real> toShift = csdMatrix[2 * i].CalculateUDegrees(csdMatrix[2 * i + 1].m_iLevel);
        csdMatrix[2 * i + 2].PhaseShiftU(toShift, csdMatrix[2 * i + 1].m_iLevel);
    }
}

TArray<CSDMatrix> CSDMatrix::CSDecomposeMatrix(const QLMatrix& u, INT iLevel)
{
    assert(iLevel > 1);
    UINT halfSize = 1U << (iLevel - 1);

    QLMatrix q1, q2, c, s, v1, v2;
    u.CSD(q1, q2, c, s, v1, v2, halfSize, halfSize);

    v1.Dagger();
    v2.Dagger();
    //QLMatrix q, cs, v;
    //QLMatrix::CombinedCSD(q, cs, v, q1, q2, c, s, v1, v2, halfSize, halfSize);
    //v.Dagger();
    TArray<QLMatrix> l;
    TArray<QLMatrix> a;
    TArray<QLMatrix> r;
    //l.AddItem(q);
    //a.AddItem(cs);
    //r.AddItem(v);
    l.AddItem(q1);
    l.AddItem(q2);
    a.AddItem(c);
    r.AddItem(v1);
    r.AddItem(v2);

    CSDMatrix lMatrix;
    lMatrix.SetAsUMatrix(1, iLevel, l);
    CSDMatrix aMatrix;
    aMatrix.SetAsAMatrix(1, iLevel, a);
    CSDMatrix rMatrix;
    rMatrix.SetAsUMatrix(1, iLevel, r);
    TArray<CSDMatrix> result;
    result.AddItem(lMatrix);
    result.AddItem(aMatrix);
    result.AddItem(rMatrix);

    result = RecursiveUDecompose(result);
    PhaseShiftUAll(result);

    CSDMatrix lastCSDM = result.Pop();
    TArray<CSDMatrix> lastDecompose = lastCSDM.DecomposeUIntoRZ();
    INT lastDecomposeNum = lastDecompose.Num();
    for (INT i = 0; i < lastDecomposeNum; ++i)
    {
        result.AddItem(lastDecompose[lastDecomposeNum - 1 - i]);
    }

    //for (INT i = 0; i < result.Num(); ++i)
    //{
    //    CCString sName;
    //    sName.Format(_T("u%d"), i);
    //    result[i].PrintMe().Print(sName);
    //}
    return result;
}

QLGate QLAPI CSDDecompose(const QLMatrix& u, INT iLevel)
{
    assert(iLevel > 0);
    if (iLevel <= 0)
    {
        return QLGate();
    }
    if (1 == iLevel)
    {
        return CreateZYZGate(u);
    }
    assert(u.X() >= (1U << iLevel));
    assert(u.Y() >= (1U << iLevel));

    QLGate ret;
    ret.AddQubits(iLevel);
    ret.m_sName = _T("CSD");
    TArray<CSDMatrix> decompose = CSDMatrix::CSDecomposeMatrix(u, iLevel);
    for (INT i = 0; i < decompose.Num(); ++i)
    {
        INT idx = decompose.Num() - 1 - i;
        if (ECSDMatrix::AMatrix == decompose[idx].m_eType)
        {
            decompose[idx].CombinedYZGate(ret, decompose[idx - 1].m_lstDegreeContent);
            //DebugGate(decompose[idx], decompose[idx - 1], iLevel);

            //DebugGate(decompose[idx - 1], iLevel);
            //DebugGate(decompose[idx], iLevel);
        }
        else if (ECSDMatrix::FRZMatrix != decompose[idx].m_eType)
        {
            decompose[idx].Gate(ret);
            //DebugGate(decompose[idx], iLevel);
        }
        //DebugGate(decompose[idx], iLevel);
    }

    return ret;
}

void QLAPI DebugGate(const CSDMatrix& m, INT iLevel)
{
    QLGate ch;
    ch.AddQubits(iLevel);
    ch.m_sName = _T("CSD");

    m.Gate(ch);

    appGeneral(_T("matrix type: %d\n"), m.m_eType);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = iLevel;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);

    m.PrintMe().Print("m");
}

void QLAPI DebugGate(const CSDMatrix& m1, const CSDMatrix& m2, INT iLevel)
{
    QLGate ch;
    ch.AddQubits(iLevel);
    ch.m_sName = _T("CSD");

    m1.CombinedYZGate(ch, m2.m_lstDegreeContent);

    appGeneral(_T("matrix type: %d\n"), m1.m_eType);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = iLevel;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);

    (m1.PrintMe() * m2.PrintMe()).Print("m");
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================