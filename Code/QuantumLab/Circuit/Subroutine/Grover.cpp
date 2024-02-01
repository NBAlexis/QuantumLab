//=============================================================================
// FILENAME : Grover.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/28/2023 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE


/**
*
*/
QLGate QLAPI GroverSXGate(BYTE numberOfQubits, UINT elementToFlip, UBOOL bUseAncilla)
{
    QLGate ret;
    if (bUseAncilla)
    {
        ret.AddQubits(numberOfQubits + 2);
    }
    else
    {
        ret.AddQubits(numberOfQubits + 1);
    }

    QLGate x(EBasicOperation::EBO_X);

    for (BYTE i = 0; i < numberOfQubits; ++i)
    {
        if (0 == (elementToFlip & (1 << i)))
        {
            TArray<BYTE> addX;
            addX.AddItem(i);
            ret.AppendGate(x, addX);
        }
    }

    if (bUseAncilla)
    {
        QLGate cnotWithAncilla = CreateCnNotWithAncilla(numberOfQubits);
        TArray<BYTE> addCnot;
        if (numberOfQubits > 2)
        {
            //Only when number of qubits > 2, one ancilla is used
            addCnot.Append(ByteSequnce, numberOfQubits + 2);
        }
        else
        {
            addCnot.Append(ByteSequnce, numberOfQubits + 1);
        }
        ret.AppendGate(cnotWithAncilla, addCnot);
    }
    else
    {
        QLGate cnotWithAncilla = CreateCnNot(numberOfQubits);
        TArray<BYTE> addCnot;
        addCnot.Append(ByteSequnce, numberOfQubits + 1);
        ret.AppendGate(cnotWithAncilla, addCnot);
    }

    for (BYTE i = 0; i < numberOfQubits; ++i)
    {
        if (0 == (elementToFlip & (1 << i)))
        {
            TArray<BYTE> addX;
            addX.AddItem(i);
            ret.AppendGate(x, addX);
        }
    }

    return ret;
}

/**
*
*/
QLGate QLAPI AmplitudeAmplification(QLGate agate, TArray<BYTE> subspaceQubits, UINT elementToFlip, UINT uiRepeat, UBOOL bUseAncilla, UINT uiInitialState)
{
    QLGate totalWithAA;
    UINT totalQubits = agate.GetQubitCount() + (bUseAncilla ? 2 : 1);
    totalWithAA.AddQubits(static_cast<BYTE>(totalQubits));
    TArray<BYTE> addAQubits;
    addAQubits.Append(ByteSequnce, agate.GetQubitCount());
    TArray<BYTE> sxQubits;
    sxQubits.Append(ByteSequnce, agate.GetQubitCount() + (bUseAncilla ? 2 : 1));
    QLGate s0Gate = GroverSXGate(agate.GetQubitCount(), uiInitialState, bUseAncilla);

    QLGate s0GateSubspace = GroverSXGate(static_cast<BYTE>(subspaceQubits.Num()), 0, bUseAncilla);
    subspaceQubits.AddItem(agate.GetQubitCount());
    if (bUseAncilla)
    {
        subspaceQubits.AddItem(agate.GetQubitCount() + 1);
    }

    QLGate x(EBasicOperation::EBO_X);
    QLGate h(EBasicOperation::EBO_H);
    TArray<BYTE> ancillaQubit;
    ancillaQubit.AddItem(agate.GetQubitCount());
    totalWithAA.AppendGate(agate, addAQubits);
    totalWithAA.AppendGate(x, ancillaQubit);
    totalWithAA.AppendGate(h, ancillaQubit);

    for (UINT i = 0; i < uiRepeat; ++i)
    {
        totalWithAA.AppendGate(s0GateSubspace, subspaceQubits);
        agate.Dagger();
        totalWithAA.AppendGate(agate, addAQubits);
        totalWithAA.AppendGate(s0Gate, sxQubits);
        agate.Dagger();
        totalWithAA.AppendGate(agate, addAQubits);
    }

    totalWithAA.AppendGate(h, ancillaQubit);
    totalWithAA.AppendGate(x, ancillaQubit);

    return totalWithAA;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================