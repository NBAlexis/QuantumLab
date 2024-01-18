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

void QLSimulatorDensityMatrix::Simulate(QLSimulatorParameters * params, QLSimulatorOutput* out) const
{
    QLSimulatorParametersDensityMatrix* param = dynamic_cast<QLSimulatorParametersDensityMatrix*>(params);
    QLSimulatorOutputDensityMatrix* output = dynamic_cast<QLSimulatorOutputDensityMatrix*>(out);
    TArray<BYTE> qubits;
    for (BYTE byQ = 0; byQ < param->m_byQubitCount; ++byQ)
    {
        qubits.AddItem(byQ);
    }
    TArray<SBasicOperation> ops = param->m_MasterGate.GetOperation(qubits);
    SIZE_T opssize = ops.Num();

    //appGeneral(_T("%d gates to apply!!\n"), opssize);

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createDensityQureg(param->m_byQubitCount, evn);
    Qureg vec2;
    if (param->m_bMeasureFidelity)
    {
        vec2 = createQureg(param->m_byQubitCount, evn);
    }

    syncQuESTEnv(evn);

    //LONGLONG bufferlength = param->BuildZeroStart(param->m_byQubitCount, vec.stateVec.real, vec.stateVec.imag);

    copyStateToGPU(vec);

    UINT uiSingleQubitGate = 0;
    UINT uiTwoQubitGate = 0;

    for (SIZE_T i = 0; i < opssize; ++i)
    {
        SBasicOperation& op = ops[static_cast<INT>(i)];
        if (op.IsSingleQubitGate())
        {
            ++uiSingleQubitGate;
        }
        if (op.IsTwoQubitGate())
        {
            ++uiTwoQubitGate;
            if (!op.IsBasicGate())
            {
                appParanoiac(_T("there is no basic gate: %d\n"), op.m_eOperation);
            }
        }

        Real fProba = QLGate::PerformBasicOperation(vec, op);
        if (NULL != output)
        {
            output->m_fProbability *= fProba;
        }

        if (NULL != param->m_pCallBack)
        {
            (*param->m_pCallBack)(static_cast<UINT>(i), fProba, ops[static_cast<INT>(i)]);
        }

        if (param->m_bMeasureFidelity)
        {
            QLGate::PerformBasicOperation(vec2, op);
        }

        if (!op.IsNoise() && !op.IsMeasure())
        {
            if (param->m_fDampingAfterGate > _QL_FLT_MIN)
            {
                for (INT qubitIdx = 0; qubitIdx < op.m_lstQubits.Num(); ++qubitIdx)
                {
                    mixDamping(vec, op.m_lstQubits[qubitIdx], param->m_fDampingAfterGate);
                }
            }

            if (param->m_fDephaseAfterGate > _QL_FLT_MIN)
            {
                if (param->m_fTwoDephaseAfterGate < _QL_FLT_MIN)
                {
                    for (INT qubitIdx = 0; qubitIdx < op.m_lstQubits.Num(); ++qubitIdx)
                    {
                        mixDephasing(vec, op.m_lstQubits[qubitIdx], param->m_fDephaseAfterGate);
                    }
                }
                else
                {
                    if (1 == op.m_lstQubits.Num())
                    {
                        mixDephasing(vec, op.m_lstQubits[0], param->m_fDephaseAfterGate);
                    }
                }
            }
            if (param->m_fDepolarisingAfterGate > _QL_FLT_MIN)
            {
                if (param->m_fTwoDepolarisingAfterGate < _QL_FLT_MIN)
                {
                    for (INT qubitIdx = 0; qubitIdx < op.m_lstQubits.Num(); ++qubitIdx)
                    {
                        mixDepolarising(vec, op.m_lstQubits[qubitIdx], param->m_fDepolarisingAfterGate);
                    }
                }
                else
                {
                    if (1 == op.m_lstQubits.Num())
                    {
                        mixDepolarising(vec, op.m_lstQubits[0], param->m_fDepolarisingAfterGate);
                    }
                }
            }
            if (param->m_fMixPauliAfterGate > _QL_FLT_MIN)
            {
                for (INT qubitIdx = 0; qubitIdx < op.m_lstQubits.Num(); ++qubitIdx)
                {
                    mixPauli(vec, op.m_lstQubits[qubitIdx], param->m_fMixPauliAfterGate, param->m_fMixPauliAfterGate, param->m_fMixPauliAfterGate);
                }
            }
            if (param->m_fTwoDepolarisingAfterGate > _QL_FLT_MIN)
            {
                if (2 == op.m_lstQubits.Num())
                {
                    mixTwoQubitDepolarising(vec, op.m_lstQubits[0], op.m_lstQubits[1], param->m_fTwoDepolarisingAfterGate);
                }
            }
            if (param->m_fTwoDephaseAfterGate > _QL_FLT_MIN)
            {
                if (2 == op.m_lstQubits.Num())
                {
                    mixTwoQubitDephasing(vec, op.m_lstQubits[0], op.m_lstQubits[1], param->m_fTwoDephaseAfterGate);
                }
            }
        }
    }
    syncQuESTEnv(evn);

    appGeneral(_T("qubits: %d, operator number: %d, single qubit gate: %d, two qubit gate: %d\n"), param->m_byQubitCount, opssize, uiSingleQubitGate, uiTwoQubitGate);

    if (param->m_iMeasureTimes > 0)
    {
        copyStateFromGPU(vec);
        syncQuESTEnv(evn);

        //Real* resreal = reinterpret_cast<Real*>(malloc(sizeof(Real) * bufferlength));
        //Real* resimag = reinterpret_cast<Real*>(malloc(sizeof(Real) * bufferlength));
        //if (NULL == resreal || NULL == resimag)
        //{
        //    appCrucial("buffer not created!");
        //    return;
        //}

        //memcpy(resreal, vec.stateVec.real, sizeof(Real) * bufferlength);
        //memcpy(resimag, vec.stateVec.imag, sizeof(Real) * bufferlength);

        TArray<UINT> lstCount;
        for (UINT i = 0; i < (1U << param->m_lstMeasureQubits.Num()); ++i)
        {
            lstCount.AddItem(0);
        }

        for (INT i = 0; i < param->m_iMeasureTimes; ++i)
        {
            if (0 != i)
            {
                //memcpy(vec.stateVec.real, resreal, sizeof(Real)* bufferlength);
                //memcpy(vec.stateVec.imag, resimag, sizeof(Real)* bufferlength);
                copyStateToGPU(vec);
                syncQuESTEnv(evn);
            }

            if (param->m_fDampingBeforeMeasure > _QL_FLT_MIN)
            {
                for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
                {
                    mixDamping(vec, qubit, param->m_fDampingBeforeMeasure);
                }
            }
            if (param->m_fDepolarisingBeforeMeasure > _QL_FLT_MIN)
            {
                for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
                {
                    mixDepolarising(vec, qubit, param->m_fDepolarisingBeforeMeasure);
                }
            }
            if (param->m_fDephaseBeforeMeasure > _QL_FLT_MIN)
            {
                for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
                {
                    mixDephasing(vec, qubit, param->m_fDephaseBeforeMeasure);
                }
            }
            if (param->m_fMixPauliBeforeMeasure > _QL_FLT_MIN)
            {
                for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
                {
                    mixPauli(vec, qubit, param->m_fMixPauliBeforeMeasure, param->m_fMixPauliBeforeMeasure, param->m_fMixPauliBeforeMeasure);
                }
            }

            UINT measureRes = 0;
            for (INT j = 0; j < param->m_lstMeasureQubits.Num(); ++j)
            {
                INT outDigit = measure(vec, param->m_lstMeasureQubits[j]);
                if (1 == outDigit)
                {
                    measureRes = measureRes | (1U << j);
                }
            }
            lstCount[measureRes] = lstCount[measureRes] + 1;
        }

        MeasurePurity(param, output, &vec);
        MeasureFidelity(param, output, &vec, &vec2);

        for (UINT i = 0; i < (1U << param->m_lstMeasureQubits.Num()); ++i)
        {
            if (NULL != output)
            {
                output->m_lstMeasureOutcomes.AddItem(lstCount[i] / static_cast<Real>(param->m_iMeasureTimes));
            }
            if (param->m_bPrint)
            {
                appGeneral(_T("%s(%d): %f\n"), 
                    Binary(i, param->m_lstMeasureQubits.Num()).c_str(), 
                    i, 
                    lstCount[i] / static_cast<Real>(param->m_iMeasureTimes));
            }
        }

        //appSafeFree(resreal);
        //appSafeFree(resimag);
    }
    else
    {
        if (param->m_fDampingBeforeMeasure > _QL_FLT_MIN)
        {
            for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
            {
                mixDamping(vec, qubit, param->m_fDampingBeforeMeasure);
            }
        }
        if (param->m_fDepolarisingBeforeMeasure > _QL_FLT_MIN)
        {
            for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
            {
                mixDepolarising(vec, qubit, param->m_fDepolarisingBeforeMeasure);
            }
        }
        if (param->m_fDephaseBeforeMeasure > _QL_FLT_MIN)
        {
            for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
            {
                mixDephasing(vec, qubit, param->m_fDephaseBeforeMeasure);
            }
        }
        if (param->m_fMixPauliBeforeMeasure > _QL_FLT_MIN)
        {
            for (BYTE qubit = 0; qubit < param->m_byQubitCount; ++qubit)
            {
                mixPauli(vec, qubit, param->m_fMixPauliBeforeMeasure, param->m_fMixPauliBeforeMeasure, param->m_fMixPauliBeforeMeasure);
            }
        }

        MeasurePurity(param, output, &vec);
        MeasureFidelity(param, output, &vec, &vec2);

        TArray<INT> tomeasure;
        for (INT i = 0; i < param->m_lstMeasureQubits.Num(); ++i)
        {
            tomeasure.AddItem(static_cast<BYTE>(param->m_lstMeasureQubits[i]));
        }

        Real* res = (Real*)malloc(sizeof(Real) * (1ULL << static_cast<UINT>(param->m_lstMeasureQubits.Num())));
        calcProbOfAllOutcomes(res, vec, tomeasure.GetData(), tomeasure.Num());
        if (NULL != output)
        {
            output->m_lstMeasureOutcomes.Append(res, 1U << param->m_lstMeasureQubits.Num());
        }

        if (param->m_bPrint)
        {
            for (UINT i = 0; i < (1U << param->m_lstMeasureQubits.Num()); ++i)
            {
                appGeneral(_T("%s(%d): %f\n"),
                    Binary(i, param->m_lstMeasureQubits.Num()).c_str(),
                    i,
                    res[i]);
            }
        }

        appSafeFree(res);
    }

    destroyQureg(vec, evn);
    if (param->m_bMeasureFidelity)
    {
        destroyQureg(vec2, evn);
    }
    destroyQuESTEnv(evn);
}

void QLSimulatorDensityMatrix::MeasurePurity(const QLSimulatorParametersDensityMatrix* param, QLSimulatorOutputDensityMatrix* output, const Qureg* vec)
{
    if (param->m_bMeasurePurity)
    {
        Real fPurity = calcPurity(*vec);
        if (NULL != output)
        {
            output->m_lstPurity.AddItem(fPurity);
        }
        if (param->m_bPrint)
        {
            appGeneral(_T("Purity: %f\n"), fPurity);
        }
    }
}

void QLSimulatorDensityMatrix::MeasureFidelity(const QLSimulatorParametersDensityMatrix* param, QLSimulatorOutputDensityMatrix* output, const Qureg* v1, const Qureg* v2)
{
    if (param->m_bMeasureFidelity)
    {
        Real fPurity = calcFidelity(*v1, *v2);
        if (NULL != output)
        {
            output->m_lstFidelity.AddItem(fPurity);
        }
        if (param->m_bPrint)
        {
            appGeneral(_T("Fidelity: %f\n"), fPurity);
        }
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================