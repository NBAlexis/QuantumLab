//=============================================================================
// FILENAME : CAdam.h
// 
// DESCRIPTION:
// 
// 1412.6980
//
// REVISION: [dd/mm/yy]
//  [03/02/2024 nbale]
//=============================================================================

#ifndef _CADAM_H_
#define _CADAM_H_

__BEGIN_NAMESPACE

class QLAPI CAdam : public COptimizer
{
public:

    CAdam(class CAnsatz* pAnsatz, CLossFunction* pLoss, Real fLearnRate = F(0.001), Real fEps = F(1.0e-8),
        Real fBeta1 = F(0.9), Real fBeta2 = F(0.999), Real fGradientOffset = F(0.0001));

    Real Start() override;
    Real Iteration() override;

protected:

    Real m_fEta;
    Real m_fLastLoss;
    Real m_fEps;
    Real m_fBeta1;
    Real m_fBeta2;
    TArray<Real> m_m;
    TArray<Real> m_v;
    UINT m_uiItr;
};

__END_NAMESPACE


#endif //#ifndef _CADAM_H_

//=============================================================================
// END OF FILE
//=============================================================================