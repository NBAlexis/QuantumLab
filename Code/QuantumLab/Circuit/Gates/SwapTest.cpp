//=============================================================================
// FILENAME : SwapTest.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [14/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLGate QLAPI CreateSwapTest(TArray<QLComplex> v1, TArray<QLComplex> v2)
{
    UINT l = static_cast<UINT>(v1.Num() > v2.Num() ? v1.Num() : v2.Num());
    BYTE qubitForEachVector = MostSignificantPowerTwo(l);
    for (UINT i = 0; i < (1U << qubitForEachVector); ++i)
    {
        if (v1.Num() <= static_cast<INT>(i))
        {
            v1.AddItem(_zeroc);
        }

        if (v2.Num() <= static_cast<INT>(i))
        {
            v2.AddItem(_zeroc);
        }
    }

    QLGate ret;
    ret.AddQubits(2 * qubitForEachVector + 1);

    QLGate swap = CreateSwapGate();
    QLGate cswap = swap.CreateControlled();

    QLGate aev1 = AmplitudeEncodeOneVector(v1.GetData(), qubitForEachVector, FALSE);
    QLGate aev2 = AmplitudeEncodeOneVector(v2.GetData(), qubitForEachVector, FALSE);

    QLGate h(EBasicOperation::EBO_H);

    TArray<BYTE> aequbit1;
    TArray<BYTE> aequbit2;

    for (BYTE byi = 0; byi < qubitForEachVector; ++byi)
    {
        aequbit1.AddItem(byi + 1);
        aequbit2.AddItem(byi + 1 + qubitForEachVector);
    }

    ret.AppendGate(h, 0);
    ret.AppendGate(aev1, aequbit1);
    ret.AppendGate(aev2, aequbit2);

    for (BYTE i = 0; i < qubitForEachVector; ++i)
    {
        TArray<BYTE> qubits;
        qubits.AddItem(0);
        qubits.AddItem(i + 1);
        qubits.AddItem(i + qubitForEachVector + 1);

        ret.AppendGate(cswap, qubits);
    }

    ret.AppendGate(h, 0);

    return ret;
}

QLGate QLAPI CreateSwapTestReal(TArray<QLComplex> v1, TArray<QLComplex> v2)
{
    UINT l = static_cast<UINT>(v1.Num() > v2.Num() ? v1.Num() : v2.Num());
    BYTE qubitForEachVector = MostSignificantPowerTwo(l);
    for (UINT i = 0; i < (1U << qubitForEachVector); ++i)
    {
        if (v1.Num() <= static_cast<INT>(i))
        {
            v1.AddItem(_zeroc);
        }

        if (v2.Num() <= static_cast<INT>(i))
        {
            v2.AddItem(_zeroc);
        }
    }

    QLGate ret;
    ret.AddQubits(2 * qubitForEachVector + 1);

    QLGate swap = CreateSwapGate();
    QLGate cswap = swap.CreateControlled();

    QLGate aev1 = AmplitudeEncodeOneVectorReal(v1.GetData(), qubitForEachVector);
    QLGate aev2 = AmplitudeEncodeOneVectorReal(v2.GetData(), qubitForEachVector);

    QLGate h(EBasicOperation::EBO_H);

    TArray<BYTE> aequbit1;
    TArray<BYTE> aequbit2;

    for (BYTE byi = 0; byi < qubitForEachVector; ++byi)
    {
        aequbit1.AddItem(byi + 1);
        aequbit2.AddItem(byi + 1 + qubitForEachVector);
    }

    ret.AppendGate(h, 0);
    ret.AppendGate(aev1, aequbit1);
    ret.AppendGate(aev2, aequbit2);

    for (BYTE i = 0; i < qubitForEachVector; ++i)
    {
        TArray<BYTE> qubits;
        qubits.AddItem(0);
        qubits.AddItem(i + 1);
        qubits.AddItem(i + qubitForEachVector + 1);

        ret.AppendGate(cswap, qubits);
    }

    ret.AppendGate(h, 0);

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================