//=============================================================================
// FILENAME : CHamitonianList.h
// 
// DESCRIPTION:
// This is the sum of hamitonian terms
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CHAMITONIANLIST_H_
#define _CHAMITONIANLIST_H_

__BEGIN_NAMESPACE

enum class ETimeDecomposeType
{
    ETDT_Trotter,
    ETDT_ABCBA,
};

class QLAPI CHamitonianList
{
public:

    CHamitonianList(ETimeDecomposeType eDecompose = ETimeDecomposeType::ETDT_ABCBA) : m_pLattice(NULL), m_eDecompose(eDecompose) {}

    virtual ~CHamitonianList();

    QLGate BuildSimulationCircuit(Real fTime, UINT uiTrotterStep) const;
    Real Measure(const Real* hostWaveFunctionReal, const Real* hostWaveFunctionImagin, INT iRepeat) const;

    void SetDecomposeType(ETimeDecomposeType eDecompose)
    {
        m_eDecompose = eDecompose;
    }

protected:

    TArray<CHamitonianTerm*> m_lstTerms;
    CLattice* m_pLattice;
    ETimeDecomposeType m_eDecompose;

    QLGate BuildSimulationCircuitTrotter(Real fTime, UINT uiTrotterStep) const;
    virtual QLGate BuildSimulationCircuitABCBA(Real fTime, UINT uiTrotterStep) const;
};

/**
* sz(n) sz(n+1) + lambda sx(n)
*/
class QLAPI CHamitonianTIM1D : public CHamitonianList
{
public:
    CHamitonianTIM1D(UINT uiSiteCount, Real fLambda);
    void SetLambda(Real fLambda);
};

/**
* a [psibar (i gamma1 partial1 + m) psi 
* + m psibar psi
* + g0 psibar gamma0 psi
* + g1 psibar gamma1 psi
* + f1sq (psibar psi)^2]
* 
* Note that
* g5 is not implemented
* g5.g0 = -g1
* g5.g1 = -g0
* 
* g5sq = f1sq = -g1sq
* g0sq = f1sq + g0
*/
struct QLAPI SJordanWeigner1DTerms
{
    Real m_fLatticeSpacing;
    Real m_fMass;
    Real m_fG0;
    Real m_fG1;
    Real m_f1Sq;

    UBOOL m_bHasMass;
    UBOOL m_bHasG0;
    UBOOL m_bHasG1;
    UBOOL m_bHasF1sq;
};

class QLAPI CHamitonianFermion1D : public CHamitonianList
{
public:

    CHamitonianFermion1D(UINT uiSiteCount, const SJordanWeigner1DTerms& param);

    SJordanWeigner1DTerms m_sParam;
    Real m_fLatticeSpacing;

protected:

    QLGate BuildSimulationCircuitABCBA(Real fTime, UINT uiTrotterStep) const override;

    void QLGateA(QLGate& gate, Real fTime, UINT uiControllerCount) const;
    void QLGateB(QLGate& gate, Real fTime, UINT uiControllerCount) const;
    void QLGateC(QLGate& gate, Real fTime, UINT uiControllerCount) const;

};


__END_NAMESPACE


#endif //#ifndef _CHAMITONIANLIST_H_

//=============================================================================
// END OF FILE
//=============================================================================