//=============================================================================
// FILENAME : QLSimulatorMatrix.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [10/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

void QLSimulatorMatrix::Simulate(QLSimulatorParameters * params, QLSimulatorOutput* output) const
{
    const QLSimulatorParametersMatrix* param = dynamic_cast<const QLSimulatorParametersMatrix*>(params);
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
    QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex)* veclen * veclen));
    if (NULL == res)
    {
        appCrucial("buffer not created!");
        return;
    }

    for (LONGLONG line = 0; line < veclen; ++line)
    {
        syncQuESTEnv(evn);
        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            if (line2 == line)
            {
                vec.stateVec.real[line2] = 1.0;
                vec.stateVec.imag[line2] = 0.0;
            }
            else
            {
                vec.stateVec.real[line2] = 0.0;
                vec.stateVec.imag[line2] = 0.0;
            }
        }
        copyStateToGPU(vec);

        for (SIZE_T i = 0; i < opssize; ++i)
        {
            const Real fProba = QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);

            if (NULL != param->m_pCallBack)
            {
                (*param->m_pCallBack)(static_cast<UINT>(i), fProba, ops[static_cast<INT>(i)]);
            }
        }
        syncQuESTEnv(evn);
        copyStateFromGPU(vec);

        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            res[line * veclen + line2].x = static_cast<Real>(vec.stateVec.real[line2]);
            res[line * veclen + line2].y = static_cast<Real>(vec.stateVec.imag[line2]);
        }
    }
    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    QLMatrix resmtr(static_cast<UINT>(veclen), static_cast<UINT>(veclen), res);
    if (param->m_bPrint)
    {
        resmtr.Print();
    }

    QLSimulatorOutputMatrix* outputMatrix = dynamic_cast<QLSimulatorOutputMatrix*>(output);
    if (NULL != outputMatrix)
    {
        outputMatrix->m_OutputMatrix = resmtr;
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================