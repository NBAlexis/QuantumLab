//=============================================================================
// FILENAME : CAdam.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [03/02/2024 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

CAdam::CAdam(class CAnsatz* pAnsatz, CLossFunction* pLoss, Real fLearnRate,
    Real fEps, Real fBeta1, Real fBeta2, Real fGradientOffset)
    : COptimizer(pAnsatz, pLoss, fGradientOffset)
    , m_fEta(fLearnRate)
    , m_fLastLoss(F(0.0))
    , m_fEps(fEps)
    , m_fBeta1(fBeta1)
    , m_fBeta2(fBeta2)
    , m_uiItr(0)
{

}

Real CAdam::Start()
{
    m_fLastLoss = LossFunction();
    m_m.RemoveAll();
    m_v.RemoveAll();

    for (UINT i = 0; i < m_pAnsatzToOptimize->ParameterCount(); ++i)
    {
        m_m.AddItem(F(0.0));
        m_v.AddItem(F(0.0));
    }
    m_uiItr = 0;

    return m_fLastLoss;
}

/**
* t = 0
* m(t) = 0
* v(t) = 0
* while not converge:
*   t = t + 1
*   g = gradient
*   m(t) = beta1 m(t-1) + (1 - beta1) g
*   v(t) = beta2 v(t-1) + (1 - beta1) g * g
*   mhat = m(t) / (1 - beta1^t)
*   vhat = v(t) / (1 - beta2^t)
*   x = x - learnrate * mhat / (sqrt(vhat) + eps)
*   
*/
Real CAdam::Iteration()
{
    ++m_uiItr;
    TArray<Real> gradients = GetGradientsWithKnown(m_fLastLoss);
    TArray<Real> mhat;
    TArray<Real> vhat;
    for (UINT p = 0; p < m_pAnsatzToOptimize->ParameterCount(); ++p)
    {
        m_m[p] = m_fBeta1 * m_m[p] + (1 - m_fBeta1) * gradients[p];
        m_v[p] = m_fBeta2 * m_v[p] + (1 - m_fBeta2) * gradients[p] * gradients[p];
        mhat.AddItem(static_cast<Real>(m_m[p] / (1.0 - pow(m_fBeta1, m_uiItr))));
        vhat.AddItem(static_cast<Real>(m_v[p] / (1.0 - pow(m_fBeta2, m_uiItr))));
    }

    for (UINT p = 0; p < m_pAnsatzToOptimize->ParameterCount(); ++p)
    {
        m_pAnsatzToOptimize->SetParameterWithOffset(p, -m_fEta * mhat[p] / (sqrt(vhat[p]) + m_fEps));
    }
    m_fLastLoss = LossFunction();
    return m_fLastLoss;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================