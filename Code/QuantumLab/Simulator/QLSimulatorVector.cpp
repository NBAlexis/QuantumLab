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
    //QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex)* veclen));
    //if (NULL == res)
    //{
    //    appCrucial("buffer not created!");
    //    return;
    //}

    syncQuESTEnv(evn);
    if (param->m_bHasInitial)
    {
        memcpy(vec.stateVec.real, param->m_pRealWavefunction, sizeof(Real) * veclen);
        memcpy(vec.stateVec.imag, param->m_pImageWavefunction, sizeof(Real) * veclen);
    }
    else
    {
        param->BuildZeroStart(param->m_byQubitCount, vec.stateVec.real, vec.stateVec.imag);
    }
    copyStateToGPU(vec);

    for (SIZE_T i = 0; i < opssize; ++i)
    {
        Real fProba = QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        if (NULL != output)
        {
            output->m_fProbability *= fProba;
        }

        if (NULL != param->m_pCallBack)
        {
            (*param->m_pCallBack)(static_cast<UINT>(i), fProba, ops[static_cast<INT>(i)]);
        }
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    QLSimulatorOutputVector* outputVector = dynamic_cast<QLSimulatorOutputVector*>(output);
    QLComplex* res = NULL;
    if (NULL != outputVector && !outputVector->m_bOutputToBuffer)
    {
        res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * veclen));
        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            res[line2].x = static_cast<Real>(vec.stateVec.real[line2]);
            res[line2].y = static_cast<Real>(vec.stateVec.imag[line2]);
        }
    }
    
    if (NULL != outputVector)
    {
        if (outputVector->m_bOutputToBuffer)
        {
            memcpy(outputVector->m_pRealBuffer, vec.stateVec.real, sizeof(Real) * veclen);
            memcpy(outputVector->m_pImageBuffer, vec.stateVec.imag, sizeof(Real) * veclen);
        }
    }

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    if (NULL != outputVector && !outputVector->m_bOutputToBuffer)
    {
        QLMatrix resmtr(static_cast<UINT>(veclen), 1, res);
        if (NULL != outputVector)
        {
            if (!outputVector->m_bOutputToBuffer)
            {
                outputVector->m_OutputMatrix = resmtr;
            }
        }

        if (param->m_bPrint)
        {
            resmtr.Print();
        }
    }
}

QLMatrix QLSimulatorVector::ShowState(const QLGate& gate)
{
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(gate.m_lstQubits.Num(), evn);
    UINT veclen = 1UL << static_cast<UINT>(gate.m_lstQubits.Num());
    memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
    memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
    vec.stateVec.real[0] = F(1.0);
    copyStateToGPU(vec);
    TArray<SBasicOperation> ops = gate.GetOperation(gate.m_lstQubits);
    SIZE_T opssize = ops.Num();
    for (SIZE_T i = 0; i < opssize; ++i)
    {
        QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    QLMatrix state = StateToMatrix(vec);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
    return state;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================