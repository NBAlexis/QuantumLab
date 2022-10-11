//=============================================================================
// FILENAME : QuantumPhaseEstimate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [02/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLGate QLAPI QuantumPhaseEstimateWithU(const QLGate& ugate, BYTE numberOfPhaseQubit)
{
    QLGate ret;
    ret.m_sName = _T("QPE-U");
    BYTE uBitNumber = static_cast<BYTE>(ugate.m_lstQubits.Num());
    ret.AddQubits(uBitNumber + numberOfPhaseQubit);
    QLGate cu_gate = ugate.CreateControlled();
    QLGate fft = QuantumFFTGate(numberOfPhaseQubit);
    fft.Dagger();

    QLGate x(EBasicOperation::EBO_X);
    QLGate h(EBasicOperation::EBO_H);
    TArray<BYTE> ubits;
    for (BYTE i = 0; i < uBitNumber; ++i)
    {
        TArray<BYTE> onebit;
        onebit.AddItem(i);
        ret.AppendGate(x, onebit);
        ubits.AddItem(i);
    }
    
    TArray<BYTE> phasebits;
    for (BYTE i = 0; i < numberOfPhaseQubit; ++i)
    {
        TArray<BYTE> onebit;
        onebit.AddItem(i + uBitNumber);
        ret.AppendGate(h, onebit);
        phasebits.AddItem(i + uBitNumber);
    }

    for (BYTE i = 0; i < numberOfPhaseQubit; ++i)
    {
        for (BYTE j = 0; j < static_cast<BYTE>(1U << i); ++j)
        {
            TArray<BYTE> cubits;
            cubits.AddItem(i);
            cubits.Append(ubits);
            ret.AppendGate(cu_gate, cubits);
        }
    }

    ret.AppendGate(fft, phasebits);

    return ret;
}


QLGate QLAPI QuantumPhaseEstimateWithHSimple(const QLMatrix& h, Real t, UINT iTrotterStep, BYTE numberOfPhaseQubit, Real fMinimalKept)
{
    UINT lenHx = h.X();
    BYTE uQubitLength = static_cast<BYTE>(Log2(lenHx) + 1);

    QLGate qpe;
    qpe.AddQubits(uQubitLength + numberOfPhaseQubit);
    qpe.m_sName = _T("QPE");

    TArray<PauliProduct> decomposed = DecomposePauliProducts(h, fMinimalKept);
    for (INT i = 1; i < decomposed.Num(); ++i)
    {
        for (INT j = i + 1; j < decomposed.Num(); ++j)
        {
            if (abs(decomposed[i].m_fCoefficient) < abs(decomposed[j].m_fCoefficient))
            {
                PauliProduct temp = decomposed[i];
                decomposed[i] = decomposed[j];
                decomposed[j] = temp;
            }
        }
    }

    QLGate u= PauliSimulateGateLeapfrog(static_cast<BYTE>(uQubitLength - 1), decomposed, t, iTrotterStep);
    QLGate cu = u.CreateControlled();

    QLGate qft_circuit = QuantumFFTGate(numberOfPhaseQubit);
    qft_circuit.Dagger();

    QLGate xgate(EBasicOperation::EBO_X);
    QLGate hgate(EBasicOperation::EBO_H);

    TArray<BYTE> ubits;
    for (BYTE ubit = 0; ubit < uQubitLength; ++ubit)
    {
        ubits.AddItem(ubit);
        if (ubit < uQubitLength - 1)
        {
            TArray<BYTE> xbit;
            xbit.AddItem(ubit);
            qpe.AppendGate(xgate, xbit);
        }
    }
    TArray<BYTE> pbits;
    for (BYTE pbit = 0; pbit < numberOfPhaseQubit; ++pbit)
    {
        pbits.AddItem(uQubitLength + pbit);

        TArray<BYTE> hardbit;
        hardbit.AddItem(uQubitLength + pbit);
        qpe.AppendGate(hgate, hardbit);

        hardbit.Append(ubits);
        
        for (UINT loop = 0; loop < (1U << pbit); ++loop)
        {
            qpe.AppendGate(cu, hardbit);
        }
    }

    qpe.AppendGate(qft_circuit, pbits);

    return qpe;
}


QLGate QLAPI BuildImprovedQPEInitialState(BYTE byQubit)
{
    TArray<Real> coeffs;
    Real T = static_cast<Real>(1U << byQubit);
    for (UINT i = 0; i < (1U << byQubit); ++i)
    {
        coeffs.AddItem(_sqrt(2 / T) * _sin(PI * (i + F(0.5)) / T));
    }

    QLGate ret = AmplitudeEncodeReal(coeffs);
    ret.m_sName = _T("IQPEInit");
    return ret;
}

QLGate QLAPI ConditionalHamiltonianEvolution(const QLMatrix& h, BYTE phaseQubitNum, Real t, UINT shortestIntervalTrotter, Real fMinimalKept)
{
    BYTE uQubitLength = static_cast<BYTE>(Log2(h.X()) + 1);
    TArray<Real> evaluationTimes;
    TArray<BYTE> uqubits;
    for (BYTE i = 0; i < uQubitLength; ++i)
    {
        uqubits.AddItem(i);
    }

    for (UINT i = 0; i < (1U << phaseQubitNum); ++i)
    {
        evaluationTimes.AddItem(static_cast<Real>(i));
    }
    Real fevaluationFactor = t / static_cast<Real>(1U << phaseQubitNum);
    TArray<Real> theta = SpliteAngles(evaluationTimes, 1U << phaseQubitNum);

    TArray<PauliProduct> decomposed = DecomposePauliProducts(h, fMinimalKept);
    for (INT i = 1; i < decomposed.Num(); ++i)
    {
        for (INT j = i + 1; j < decomposed.Num(); ++j)
        {
            if (abs(decomposed[i].m_fCoefficient) < abs(decomposed[j].m_fCoefficient))
            {
                PauliProduct temp = decomposed[i];
                decomposed[i] = decomposed[j];
                decomposed[j] = temp;
            }
        }
    }
    appGeneral(_T("Pauli decompose %d matrix \n"), decomposed.Num());

    QLGate ret;
    ret.AddQubits(uQubitLength + phaseQubitNum);
    ret.m_sName = _T("CHE");

    UINT degreeNumber = (1U << phaseQubitNum);
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        INT trotterStep = static_cast<INT>(round(abs(theta[i] / F(0.5))));
        if (trotterStep > 0)
        {
            QLGate u = PauliSimulateGateLeapfrog(static_cast<BYTE>(uQubitLength - 1), decomposed, theta[i] * fevaluationFactor, trotterStep * shortestIntervalTrotter);
            ret.AppendGate(u, uqubits);
        }

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        //ctrIdx = uQubitLength + (phaseQubitNum - 1 - ctrIdx);
        ctrIdx = uQubitLength + ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(ctrIdx);
        cnotbits.AddItem(uQubitLength - 1);
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }
    return ret;
}

QLGate QLAPI QuantumPhaseEstimateWithHImproved(const QLMatrix& h, Real t, UINT iTrotterStep, BYTE numberOfPhaseQubit, Real fMinimalKept)
{
    t = t * (1U << numberOfPhaseQubit);
    BYTE uQubitLength = static_cast<BYTE>(Log2(h.X()) + 1);

    QLGate qpe;
    qpe.AddQubits(uQubitLength + numberOfPhaseQubit);
    qpe.m_sName = _T("QPE");

    TArray<BYTE> phaseQubits;
    for (BYTE i = 0; i < numberOfPhaseQubit; ++i)
    {
        phaseQubits.AddItem(i + uQubitLength);
    }
    QLGate preparestate = BuildImprovedQPEInitialState(numberOfPhaseQubit);
    qpe.AppendGate(preparestate, phaseQubits);

    QLGate che = ConditionalHamiltonianEvolution(h, numberOfPhaseQubit, t, iTrotterStep, fMinimalKept);
    qpe.AppendGate(che, qpe.m_lstQubits);

    QLGate qft_circuit = QuantumFFTGate(numberOfPhaseQubit);
    qft_circuit.Dagger();
    qpe.AppendGate(qft_circuit, phaseQubits);
    
    return qpe;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================