//=============================================================================
// FILENAME : QuantumFFT.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLGate QLAPI QuantumFFTGate(BYTE qubitCount)
{
    QLGate ret;
    ret.AddQubits(qubitCount);
    ret.m_sName = _T("FFT");
    for (BYTE i = 0; i < qubitCount; ++i)
    {
        TArray<BYTE> had;
        had.AddItem(qubitCount - 1 - i);
        QLGate h(EBasicOperation::EBO_H);
        ret.AppendGate(h, had);
        for (BYTE j = i + 1; j < qubitCount; ++j)
        {
            TArray<BYTE> cpq;
            cpq.AddItem(qubitCount - 1 - j);
            cpq.AddItem(qubitCount - 1 - i);

            QLGate cp(EBasicOperation::EBO_CP, PI / (1U << (j - i)));
            ret.AppendGate(cp, cpq);
        }
    }

    for (BYTE i = 0; i < qubitCount / 2; ++i)
    {
        TArray<BYTE> swq;
        swq.AddItem(i);
        swq.AddItem(qubitCount - 1 - i);

        QLGate swap = CreateSwapGate();
        ret.AppendGate(swap, swq);
    }
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================