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

QLGate QLAPI CreateSwapTest(BYTE vectorQubitNum)
{
    QLGate ret;
    ret.AddQubits(2 * vectorQubitNum);

    QLGate swap = CreateSwapGate();
    QLGate cswap = swap.CreateControlled();
    QLGate h(EBasicOperation::EBO_H);
    TArray<BYTE> hqubit;
    hqubit.AddItem(0);

    ret.AppendGate(h, hqubit);

    for (BYTE i = 0; i < vectorQubitNum; ++i)
    {
        TArray<BYTE> qubits;
        qubits.AddItem(0);
        qubits.AddItem(i + 1);
        qubits.AddItem(i + vectorQubitNum + 1);

        ret.AppendGate(cswap, qubits);
    }

    ret.AppendGate(h, hqubit);

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================