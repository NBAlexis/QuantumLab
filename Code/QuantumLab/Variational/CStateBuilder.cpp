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

TArray<Real> CStateBuilder::Fit(COptimizer* optimizer, const CCString& sAnsatzFile, UINT uiMaxStep)
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

    appParanoiac(_T("state was built\n"));
    optimizer->SetLossFunc(this);
    return optimizer->Optimize(m_fGoal, uiMaxStep, sAnsatzFile);
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

    return F(1.0) - (m_bMaeOrMse ? _cuCabsf(ansatzState.VectorDot(m_StateToFit)) : __cuCabsSqf(ansatzState.VectorDot(m_StateToFit)));
}

void QLAPI FitAE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, UINT uiLevel, 
    ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, EAnsatzInitial eAnsatzInitial,
    Real fLearnRate, Real fGoal, UINT uiMaxStep, UBOOL bUseAbsorlute)
{
    UINT w, h;
    TArray<QLComplex> points = ReadCSVA(sPointFile, w, h);
    UINT wpower = MostSignificantPowerTwo(w);
    UINT hpower = MostSignificantPowerTwo(h);

    QLGate aegate = AmplitudeEncodeVectors(points.GetData(), hpower, wpower, FALSE);
    BYTE byQubit = static_cast<BYTE>(aegate.m_lstQubits.Num());

    CTwoLocal ansatz(byQubit, uiLevel, eAnsatzSingleLayer, eAnsatzLayer, eAnsatzStyle, eAnsatzInitial);
    CAdam optimizer(&ansatz, NULL, fLearnRate);
    CStateBuilder builder(aegate, fGoal, bUseAbsorlute);
    TArray<Real> history = builder.Fit(&optimizer, sAnsatzFile, uiMaxStep);
    SaveCSVAR(history.GetData(), 1, history.Num(), sHistoryFile);
}

void QLAPI FitSE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, BYTE byEncodeQubits, UINT uiLevel, 
    ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, EAnsatzInitial eAnsatzInitial,
    ELinkStyle eSimpleEncodeStyle, ELinkLayer eSimpleencodeLayer, Real fLearnRate, Real fGoal, UINT uiMaxStep, UBOOL bUseAbsorlute)
{
    UINT w, h;
    TArray<QLComplex> points = ReadCSVA(sPointFile, w, h);
    //UINT wpower = MostSignificantPowerTwo(w);
    UINT hpower = MostSignificantPowerTwo(h);

    QLGate segate = SimpleEncodeVectorsWithLinkStype(points.GetData(), eSimpleEncodeStyle, eSimpleencodeLayer, static_cast<BYTE>(hpower), byEncodeQubits, w);
    BYTE byQubit = static_cast<BYTE>(segate.m_lstQubits.Num());

    CTwoLocal ansatz(byQubit, uiLevel, eAnsatzSingleLayer, eAnsatzLayer, eAnsatzStyle, eAnsatzInitial);
    CAdam optimizer(&ansatz, NULL, fLearnRate);
    CStateBuilder builder(segate, fGoal, bUseAbsorlute);
    TArray<Real> history = builder.Fit(&optimizer, sAnsatzFile, uiMaxStep);
    SaveCSVAR(history.GetData(), 1, history.Num(), sHistoryFile);
}

void QLAPI FitAE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, 
    ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, EAnsatzInitial eAnsatzInitial,
    Real fLearnRate, Real fGoal, UINT uiMaxStep, UINT uiMaxLayer, UINT uiAdaptiveWait, Real fAdaptiveEps, UBOOL bUseAbsorlute)
{
    UINT w, h;
    TArray<QLComplex> points = ReadCSVA(sPointFile, w, h);
    UINT wpower = MostSignificantPowerTwo(w);
    UINT hpower = MostSignificantPowerTwo(h);

    QLGate aegate = AmplitudeEncodeVectors(points.GetData(), hpower, wpower, FALSE);
    BYTE byQubit = static_cast<BYTE>(aegate.m_lstQubits.Num());

    CTwoLocalAdaptive ansatz(byQubit, eAnsatzSingleLayer, eAnsatzLayer, eAnsatzStyle, eAnsatzInitial);
    ansatz.SetMaxLayer(uiMaxLayer);
    CAdam optimizer(&ansatz, NULL, fLearnRate);
    optimizer.SetAdapetiveParameter(uiAdaptiveWait, fAdaptiveEps);
    CStateBuilder builder(aegate, fGoal, bUseAbsorlute);
    TArray<Real> history = builder.Fit(&optimizer, sAnsatzFile, uiMaxStep);
    SaveCSVAR(history.GetData(), 1, history.Num(), sHistoryFile);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================