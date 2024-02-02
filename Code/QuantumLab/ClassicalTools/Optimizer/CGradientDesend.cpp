//=============================================================================
// FILENAME : CGradientDesend.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

CGradientDesend::CGradientDesend(CAnsatz* pAnsatz, CLossFunction* pLoss, Real fLearnRate, Real fGradientOffset)
    : COptimizer(pAnsatz, pLoss, fGradientOffset)
    , m_fEta(fLearnRate)
    , m_fLastLoss(F(0.0))
{

}

Real CGradientDesend::Start()
{
    m_fLastLoss = LossFunction();
    return m_fLastLoss;
}

Real CGradientDesend::Iteration()
{
    TArray<Real> gradients = GetGradientsWithKnown(m_fLastLoss);

    for (UINT p = 0; p < m_pAnsatzToOptimize->ParameterCount(); ++p)
    {
        m_pAnsatzToOptimize->SetParameterWithOffset(p, -gradients[p] * m_fEta);
    }
    m_fLastLoss = LossFunction();
    return m_fLastLoss;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================