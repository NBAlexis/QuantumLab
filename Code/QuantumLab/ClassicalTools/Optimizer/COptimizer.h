//=============================================================================
// FILENAME : COptimizer.h
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#ifndef _COPTIMIZER_H_
#define _COPTIMIZER_H_

__BEGIN_NAMESPACE

/**
* To make it stochastic, change the loss function, when there are many learning samples
* When there are many samples, randomly pick some of them to calculate average gradient is stochastic
*/
class QLAPI CLossFunction
{
public:
    virtual ~CLossFunction() {}
    virtual Real LossFunction(const QLGate& gateOfAnsatz) = 0;
};

class QLAPI COptimizer
{
public:

    COptimizer(class CAnsatz* pAnsatz, CLossFunction* pLoss, Real fGradientOffset)
        : m_pAnsatzToOptimize(pAnsatz)
        , m_pLoss(pLoss)
        , m_fGradientOffset(fGradientOffset) 
    {

    }

    virtual Real Start() = 0;
    virtual Real Iteration() = 0;

    TArray<Real> Optimize(Real fGoal, UINT uiMaxStep = 10000);

    void SetLossFunc(CLossFunction* pLoss)
    {
        m_pLoss = pLoss;
    }

protected:

    Real LossFunction() const;
    Real GradientOnParameter(UINT uiParam) const;
    Real GradientOnParameterWithKnown(UINT uiParam, Real fLossNow) const;
    TArray<Real> GetGradients() const;
    TArray<Real> GetGradientsWithKnown(Real fLossNow) const;

    class CAnsatz* m_pAnsatzToOptimize;
    CLossFunction* m_pLoss;
    Real m_fGradientOffset;
};

__END_NAMESPACE


#endif //#ifndef _COPTIMIZER_H_

//=============================================================================
// END OF FILE
//=============================================================================