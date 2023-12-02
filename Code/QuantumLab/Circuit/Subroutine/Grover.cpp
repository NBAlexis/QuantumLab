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
* 
*/
QLGate QLAPI GroverS0Gate(BYTE numberOfQubits)
{
    QLGate ret;
    return ret;
}

/**
* 1 - 2|n><n| = diag(1,...,1,-1,1,...1) where -1 is at position n
* we need numberOfQubits + 1 qubits, where the last one is ancilla
* we asumme ancilla is at |0>-|1> already
* 
* ===========================================
* 1-step apply not to encode n on the qubits
* 
* 2-step apply hadard on ancilla
* 
* 3-step apply cn-NOT
* 
* 4-step apply not to encode n on the qubits again
* 
* Assume 2 ancilla
* bit numberOfQubits:  |0>-|1> to do the phase kick back
* bit numberOfQubits+1: Zeroed ancilla for CNOT
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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================