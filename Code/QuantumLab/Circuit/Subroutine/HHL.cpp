//=============================================================================
// FILENAME : HHL.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [07/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLGate QLAPI HHLGate(const QLMatrix& h, const TArray<QLComplex>& y, UINT trotterStep, Real maxAbsEigenValue, BYTE phaseQubitNum, Real small)
{
    UINT len;
    TArray<QLComplex> normalizedY = NormalizeV(y, len);
    Real t = PI / maxAbsEigenValue;
    Real rotationC = F(-0.000001) + PI2 / (1U << phaseQubitNum) / t;

    QLGate hhlGate;
    hhlGate.AddQubits(static_cast<BYTE>(len + 2 + phaseQubitNum));
    hhlGate.m_sName = _T("HHL");

    //========= step1: put y ==============
    QLGate yGate = AmplitudeEncode(normalizedY);
    TArray<BYTE> yBits;
    for (UINT i = 0; i < len; ++i)
    {
        yBits.AddItem(static_cast<BYTE>(i));
    }
    hhlGate.AppendGate(yGate, yBits);

    
    //========= step2: QPE on y ==============
    TArray<BYTE> qpeBits;
    QLGate qpeGate = QuantumPhaseEstimateWithHImproved(h, t, trotterStep, phaseQubitNum);
    for (UINT i = 0; i < len + 1 + phaseQubitNum; ++i)
    {
        qpeBits.AddItem(static_cast<BYTE>(i));
    }
    hhlGate.AppendGate(qpeGate, qpeBits);

    //========= step3: rotate on ancilla =========
    // phasebit start from len + 1
    // len + 1 is the lowest one

    TArray<Real> degreesToRot;
    UINT phaseSeparation = 1U << phaseQubitNum;
    UINT halfPhase = phaseSeparation >> 1U;
    for (UINT i = 0; i < phaseSeparation; ++i)
    {
        if (0 == i)
        {
            degreesToRot.AddItem(PI);
            continue;
        }
        Real lambdaToInverse = i / static_cast<Real>(phaseSeparation);
        if (i > halfPhase)
        {
            lambdaToInverse = lambdaToInverse - F(1.0);

        }
        lambdaToInverse = lambdaToInverse * PI2 / t;
        Real rotDegree = _atan2(rotationC / lambdaToInverse,
            _sqrt(1 - rotationC * rotationC / (lambdaToInverse * lambdaToInverse)));
        degreesToRot.AddItem(rotDegree * 2);
    }

    TArray<BYTE> rotationBits;
    for (BYTE i = 0; i < phaseQubitNum; ++i)
    {
        rotationBits.AddItem(len + phaseQubitNum - i);
    }
    rotationBits.AddItem(len + phaseQubitNum + 1);
    QLGate fry = FRy(degreesToRot, phaseQubitNum + 1);
    hhlGate.AppendGate(fry, rotationBits);

    //========= step4: collapse ancilla =========
    QLGate collaps(EBasicOperation::EBO_CC, 1);
    TArray<BYTE> collapsbit;
    collapsbit.AddItem(static_cast<BYTE>(len + 1 + phaseQubitNum));
    hhlGate.AppendGate(collaps, collapsbit);

    //========= step5: undo the qpe =========
    qpeGate.Dagger();
    hhlGate.AppendGate(qpeGate, qpeBits);
    
    return hhlGate;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================