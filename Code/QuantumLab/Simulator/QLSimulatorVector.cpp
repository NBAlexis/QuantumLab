//=============================================================================
// FILENAME : QLSimulatorVector.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

void QLSimulatorVector::Simulate(QLSimulatorParameters * params, QLSimulatorOutput* output) const
{
    QLSimulatorParametersVector* param = dynamic_cast<QLSimulatorParametersVector*>(params);
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

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    QLMatrix resmtr(static_cast<UINT>(veclen), 1, res);
    if (param->m_bPrint)
    {
        resmtr.Print();
    }

    QLSimulatorOutputVector* outputVector= dynamic_cast<QLSimulatorOutputVector*>(output);
    if (NULL != outputVector)
    {
        outputVector->m_OutputMatrix = resmtr;
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================