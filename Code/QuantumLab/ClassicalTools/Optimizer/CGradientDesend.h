//=============================================================================
// FILENAME : CGradientDesend.h
// 
// DESCRIPTION:
// 
// 
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#ifndef _CGRADIENTDESEND_H_
#define _CGRADIENTDESEND_H_

__BEGIN_NAMESPACE

class QLAPI CGradientDesend : public COptimizer
{
public:

    CGradientDesend(class CAnsatz* pAnsatz, CLossFunction* pLoss, Real fLearnRate, Real fGradientOffset = F(0.0001));

    Real Start() override;
    Real Iteration() override;

protected:

    Real m_fEta;
    Real m_fLastLoss;
};

__END_NAMESPACE


#endif //#ifndef _CGRADIENTDESEND_H_

//=============================================================================
// END OF FILE
//=============================================================================