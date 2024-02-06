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

    CStateBuilder(const QLGate& gateToBuild, Real fGoal = F(0.000001)) : m_GateToBuild(gateToBuild), m_fGoal(fGoal)
    {
    }

    TArray<Real> Fit(COptimizer* optimizer, UINT uiMaxStep = 10000);

protected:

    Real LossFunction(const QLGate& ansatzGate);

    QLGate m_GateToBuild;
    QLMatrix m_StateToFit;
    Real m_fGoal;
};

extern void QLAPI FitAE(const CCString& sPointFile, const CCString& sAnsatzFile, const CCString& sHistoryFile, 
    UINT uiLevel, Real fLearnRate, Real fGoal, UINT uiMaxStep, UBOOL bOnlyReal);

__END_NAMESPACE


#endif //#ifndef _CSTATEBUILDER_H_

//=============================================================================
// END OF FILE
//=============================================================================