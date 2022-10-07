//=============================================================================
// FILENAME : QLSimulatorMeasure.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [05/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

void QLSimulatorMeasure::Simulate(QLSimulatorParameters * params, QLSimulatorOutput* output) const
{
    //appGeneral(_T("in QLSimulatorMeasure::Simulate\n"));

    QLSimulatorParametersMeasure* param = dynamic_cast<QLSimulatorParametersMeasure*>(params);
    TArray<BYTE> qubits;
    for (BYTE byQ = 0; byQ < param->m_byQubitCount; ++byQ)
    {
        qubits.AddItem(byQ);
    }
    TArray<SBasicOperation> ops = param->m_MasterGate.GetOperation(qubits);
    SIZE_T opssize = ops.Num();

    appGeneral(_T("%d gates to apply!!\n"), opssize);

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(param->m_byQubitCount, evn);
    

    //This is a lazy slow implement, I need to use cuda to improve it
    LONGLONG veclen = 1LL << param->m_byQubitCount;
    QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex)* veclen));
    if (NULL == res)
    {
        appCrucial("buffer not created!");
        return;
    }

    //appGeneral(_T("building zero start\n"));

    if (param->m_lstStart.Num() < static_cast<INT>(veclen))
    {
        param->BuildZeroStart(param->m_byQubitCount);
    }

    syncQuESTEnv(evn);
    for (LONGLONG line2 = 0; line2 < veclen; ++line2)
    {
        vec.stateVec.real[line2] = param->m_lstStart[static_cast<INT>(line2)].x;
        vec.stateVec.imag[line2] = param->m_lstStart[static_cast<INT>(line2)].y;
    }
    copyStateToGPU(vec);
    syncQuESTEnv(evn);

    for (INT i = 0; i < opssize; ++i)
    {
        QLGate::PerformBasicOperation(vec, ops[i]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    for (LONGLONG line2 = 0; line2 < veclen; ++line2)
    {
        res[line2].x = static_cast<Real>(vec.stateVec.real[line2]);
        res[line2].y = static_cast<Real>(vec.stateVec.imag[line2]);
    }

    QLMatrix resmtr(static_cast<UINT>(veclen), 1, res);

    QLSimulatorOutputMeasure* outputMatrix = dynamic_cast<QLSimulatorOutputMeasure*>(output);
    if (NULL != outputMatrix)
    {
        outputMatrix->m_OutputMatrix = resmtr;
    }

    //================ start measure ===============
    TArray<UINT> lstCount;
    for (UINT i = 0; i < (1U << param->m_lstMeasureBits.Num()); ++i)
    {
        lstCount.AddItem(0);
    }

    const QLComplex* hostbuffer = resmtr.HostBuffer();
    for (UINT i = 0; i < param->m_iRepeat; ++i)
    {
        if (0 != i)
        {
            for (LONGLONG line2 = 0; line2 < veclen; ++line2)
            {
                vec.stateVec.real[line2] = hostbuffer[static_cast<INT>(line2)].x;
                vec.stateVec.imag[line2] = hostbuffer[static_cast<INT>(line2)].y;
            }
            copyStateToGPU(vec);
            syncQuESTEnv(evn);
        }

        UINT measureRes = 0;
        for (INT j = 0; j < param->m_lstMeasureBits.Num(); ++j)
        {
            INT out = measure(vec, param->m_lstMeasureBits[j]);
            if (1 == out)
            {
                measureRes = measureRes | (1U << j);
            }
        }
        lstCount[measureRes] = lstCount[measureRes] + 1;
    }

    if (NULL != outputMatrix)
    {
        outputMatrix->m_lstCounts = lstCount;
    }

    if (param->m_bPrint)
    {
        for (INT i = 0; i < lstCount.Num(); ++i)
        {
            appGeneral(_T("%s(%d): %f\n"), Binary(i, param->m_lstMeasureBits.Num()), i, static_cast<Real>(lstCount[i]) / param->m_iRepeat);
        }
    }

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================