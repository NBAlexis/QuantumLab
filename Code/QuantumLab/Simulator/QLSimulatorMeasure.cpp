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

    if (param->m_bPrint)
    {
        appGeneral(_T("%d gates to apply!!\n"), opssize);
    }

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(param->m_byQubitCount, evn);
    

    //This is a lazy slow implement, I need to use cuda to improve it
    LONGLONG veclen = 1LL << param->m_byQubitCount;

    syncQuESTEnv(evn);
    param->BuildZeroStart(param->m_byQubitCount, vec.stateVec.real, vec.stateVec.imag);
    copyStateToGPU(vec);
    syncQuESTEnv(evn);

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

    QLSimulatorOutputMeasure* outputMatrix = dynamic_cast<QLSimulatorOutputMeasure*>(output);

    //================ start measure ===============
    TArray<UINT> lstCount;
    TArray<UINT> lstHist;

    if (param->m_iMeasureUntil > 0)
    {
        for (UINT i = 0; i < (1U << param->m_lstMeasureBits.Num()); ++i)
        {
            lstCount.AddItem(0);
        }

        UBOOL bContinue = TRUE;
        UINT uiCount = 0;
        while (bContinue)
        {
            if (0 != uiCount)
            {
                copyStateToGPU(vec);
                syncQuESTEnv(evn);
            }

            ++uiCount;

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
            lstHist.AddItem(measureRes);

            if (static_cast<INT>(lstCount[measureRes]) >= param->m_iMeasureUntil)
            {
                if (NULL != outputMatrix)
                {
                    outputMatrix->m_uiMeasureUntilCount = uiCount;
                    outputMatrix->m_uiMeasureUntilRes = measureRes;

                    for (UINT i = 0; i < (1U << param->m_lstMeasureBits.Num()); ++i)
                    {
                        outputMatrix->m_lstMeasureOutcomes[i] = lstCount[i] / static_cast<Real>(uiCount);
                    }
                }

                bContinue = FALSE;
            }
        }
    }
    else
    {
        if (param->m_iRepeat > 0)
        {
            for (UINT i = 0; i < (1U << param->m_lstMeasureBits.Num()); ++i)
            {
                lstCount.AddItem(0);
            }

            for (UINT i = 0; i < param->m_iRepeat; ++i)
            {
                if (0 != i)
                {
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
                lstHist.AddItem(measureRes);
            }

            if (NULL != outputMatrix)
            {
                outputMatrix->m_lstCounts = lstCount;
                outputMatrix->m_lstHist = lstHist;

                for (UINT i = 0; i < (1U << param->m_lstMeasureBits.Num()); ++i)
                {
                    outputMatrix->m_lstMeasureOutcomes[i] = lstCount[i] / static_cast<Real>(param->m_iRepeat);
                }
            }
        }
        else if (NULL != outputMatrix)
        {
            Real* outcome = reinterpret_cast<Real*>(malloc(sizeof(Real) * (1U << param->m_lstMeasureBits.Num())));
            TArray<INT> tobemeasured;
            for (INT i = 0; i < param->m_lstMeasureBits.Num(); ++i)
            {
                tobemeasured.AddItem(param->m_lstMeasureBits[i]);
            }
            calcProbOfAllOutcomes(outcome, vec, tobemeasured.GetData(), param->m_lstMeasureBits.Num());
            outputMatrix->m_lstMeasureOutcomes.Append(outcome, (1U << param->m_lstMeasureBits.Num()));
            appSafeFree(outcome);
        }
    }

    if (param->m_bPrint)
    {
        for (INT i = 0; i < lstCount.Num(); ++i)
        {
            appGeneral(_T("%s(%d): %f\n"), Binary(i, param->m_lstMeasureBits.Num()).c_str(), i, static_cast<Real>(lstCount[i]) / param->m_iRepeat);
        }
    }

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================