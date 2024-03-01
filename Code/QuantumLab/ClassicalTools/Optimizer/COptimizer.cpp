//=============================================================================
// FILENAME : COptimizer.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

Real COptimizer::LossFunction() const
{
    QLGate ansatz = m_pAnsatzToOptimize->BuildStateWithParam();
    return m_pLoss->LossFunction(ansatz);
}

Real COptimizer::GradientOnParameter(UINT uiParam) const
{
    QLGate ansatz1 = m_pAnsatzToOptimize->BuildStateOffset(uiParam, m_fGradientOffset);
    Real yp = m_pLoss->LossFunction(ansatz1);
    QLGate ansatz2 = m_pAnsatzToOptimize->BuildStateOffset(uiParam, -m_fGradientOffset);
    Real ym = m_pLoss->LossFunction(ansatz2);
    return (yp - ym) * F(0.5) / m_fGradientOffset;
}

Real COptimizer::GradientOnParameterWithKnown(UINT uiParam, Real fLossNow) const
{
    QLGate ansatz1 = m_pAnsatzToOptimize->BuildStateOffset(uiParam, m_fGradientOffset);
    Real yp = m_pLoss->LossFunction(ansatz1);
    return (yp - fLossNow) / m_fGradientOffset;
}

TArray<Real> COptimizer::GetGradients() const
{
    TArray<Real> ret;
    for (UINT p = 0; p < m_pAnsatzToOptimize->ParameterCount(); ++p)
    {
        ret.AddItem(GradientOnParameter(p));
    }
    return ret;
}

TArray<Real> COptimizer::GetGradientsWithKnown(Real fLossNow) const
{
    TArray<Real> ret;
    for (UINT p = 0; p < m_pAnsatzToOptimize->ParameterCount(); ++p)
    {
        ret.AddItem(GradientOnParameterWithKnown(p, fLossNow));
    }
    return ret;
}

TArray<Real> COptimizer::Optimize(Real fGoal, UINT uiMaxStep, const CCString& smallestFileName)
{
    TArray<Real> history;
    Real fLoss = Start();
    history.AddItem(fLoss);
    UINT step = 0;
    Real fMinLoss = fLoss;
    //TArray<Real> minLossParam;
    while (fLoss > fGoal && step < uiMaxStep)
    {
        ++step;
        appGeneral(_T("optimizer iteration: step = %d, loss = %f\n"), step, fLoss);
        fLoss = Iteration();
        history.AddItem(fLoss);
        if (fLoss < fMinLoss)
        {
            fMinLoss = fLoss;
            //minLossParam = m_pAnsatzToOptimize->GetParameters();
            if (!smallestFileName.IsEmpty())
            {
                m_pAnsatzToOptimize->SaveParameters(smallestFileName);
            }
        }
        CheckAdaptive(fLoss);
    }

    if (fLoss > fGoal)
    {
        appGeneral(_T("optimizer failed with max step reached, loss = %f minloss = %f\n"), fLoss, fMinLoss);
    }
    else
    {
        appGeneral(_T("optimizer finished with loss = %f\n"), fLoss);
    }

    return history;
}

void COptimizer::CheckAdaptive(Real fLoss)
{
    if (m_fSmallestLoss - fLoss < m_fAdaptiveEps)
    {
        ++m_uiAdaptiveIteration;
        if (m_uiAdaptiveIteration >= m_uiAdaptiveIterationMax)
        {
            m_uiAdaptiveIteration = 0;
            m_pAnsatzToOptimize->IncreaseAdaptive();
            Start();
        }
    }
    else
    {
        m_uiAdaptiveIteration = 0;
        m_fSmallestLoss = fLoss;
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================