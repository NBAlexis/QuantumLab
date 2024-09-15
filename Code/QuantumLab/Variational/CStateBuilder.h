//=============================================================================
// FILENAME : CStateBuilder.h
// 
// DESCRIPTION:
// 
// Train an ansatz to mimic a gate
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#ifndef _CSTATEBUILDER_H_
#define _CSTATEBUILDER_H_

__BEGIN_NAMESPACE

class QLAPI CStateBuilder : public CLossFunction
{
public:

    CStateBuilder(const QLGate& gateToBuild, Real fGoal = F(0.000001), UBOOL bAbsorlute = TRUE) : m_GateToBuild(gateToBuild), m_fGoal(fGoal), m_bMaeOrMse(bAbsorlute)
    {
    }

    TArray<Real> Fit(COptimizer* optimizer, const CCString& sAnsatzFile, UINT uiMaxStep = 10000);

protected:

    Real LossFunction(const QLGate& ansatzGate);

    QLGate m_GateToBuild;
    QLMatrix m_StateToFit;
    Real m_fGoal;
    UBOOL m_bMaeOrMse;
};

extern void QLAPI FitAE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, UINT uiLevel, ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, EAnsatzInitial eAnsatzInitial, Real fLearnRate, Real fGoal, UINT uiMaxStep, UBOOL bUseAbsorlute);

extern void QLAPI FitSE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, BYTE byEncodeQubits, UINT uiLevel, ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, EAnsatzInitial eAnsatzInitial, ELinkStyle eSimpleEncodeStyle, ELinkLayer eSimpleencodeLayer, Real fLearnRate, Real fGoal, UINT uiMaxStep, UBOOL bUseAbsorlute);

extern void QLAPI FitAE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, EAnsatzInitial eAnsatzInitial, Real fLearnRate, Real fGoal, UINT uiMaxStep, UINT uiMaxLayer, UINT adaptiveWait, Real fadaptiveEps, UBOOL bUseAbsorlute);

__END_NAMESPACE


#endif //#ifndef _CSTATEBUILDER_H_

//=============================================================================
// END OF FILE
//=============================================================================