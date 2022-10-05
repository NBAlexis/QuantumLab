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

    QLGate u= PauliSimulateGateLeapfrog(h, t, iTrotterStep, fMinimalKept);
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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================