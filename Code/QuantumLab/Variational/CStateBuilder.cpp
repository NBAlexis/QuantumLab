//=============================================================================
// FILENAME : CStateBuilder.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

TArray<Real> CStateBuilder::Fit(COptimizer* optimizer, UINT uiMaxStep)
{
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(m_GateToBuild.m_lstQubits.Num(), evn);
    UINT veclen = 1UL << static_cast<UINT>(m_GateToBuild.m_lstQubits.Num());
    memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
    memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
    vec.stateVec.real[0] = F(1.0);
    copyStateToGPU(vec);
    TArray<SBasicOperation> ops = m_GateToBuild.GetOperation(m_GateToBuild.m_lstQubits);
    SIZE_T opssize = ops.Num();
    for (SIZE_T i = 0; i < opssize; ++i)
    {
        QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    m_StateToFit = StateToMatrix(vec);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    optimizer->SetLossFunc(this);
    return optimizer->Optimize(m_fGoal, uiMaxStep);
}

Real CStateBuilder::LossFunction(const QLGate& ansatzGate)
{
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(ansatzGate.m_lstQubits.Num(), evn);
    UINT veclen = 1UL << static_cast<UINT>(m_GateToBuild.m_lstQubits.Num());
    memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
    memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
    vec.stateVec.real[0] = F(1.0);
    copyStateToGPU(vec);
    TArray<SBasicOperation> ops = ansatzGate.GetOperation(ansatzGate.m_lstQubits);
    SIZE_T opssize = ops.Num();
    for (SIZE_T i = 0; i < opssize; ++i)
    {
        QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    QLMatrix ansatzState = StateToMatrix(vec);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    return F(1.0) - _cuCabsf(ansatzState.VectorDot(m_StateToFit));
}

void QLAPI FitAE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, 
    UINT uiLevel, Real fLearnRate, Real fGoal, UINT uiMaxStep, UBOOL bOnlyReal)
{
    UINT w, h;
    TArray<QLComplex> points = ReadCSVA(sPointFile, w, h);
    UINT wpower = MostSignificantPowerTwo(w);
    UINT hpower = MostSignificantPowerTwo(h);

    QLGate aegate = AmplitudeEncodeVectors(points.GetData(), hpower, wpower, FALSE);
    BYTE byQubit = static_cast<BYTE>(aegate.m_lstQubits.Num());

    ESingleLayer eSingle = ESingleLayer::RYRZ;
    ELinkLayer eLinker = ELinkLayer::CZ;
    if (bOnlyReal)
    {
        eSingle = ESingleLayer::RY;
    }

    CTwoLocal ansatz(byQubit, uiLevel, eSingle, eLinker, ELinkStyle::SCA);
    CAdam optimizer(&ansatz, NULL, fLearnRate);
    CStateBuilder builder(aegate, fGoal);
    TArray<Real> history = builder.Fit(&optimizer, uiMaxStep);
    ansatz.SaveParameters(sAnsatzFile);
    SaveCSVAR(history.GetData(), 1, history.Num(), sHistoryFile);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================