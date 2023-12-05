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
* psibar (i gamma1 partial1 + m) psi - g (psibar psi)^2 - mu psibar gamma0 psi
*/
class QLAPI CHamitonianNJL1D : public CHamitonianList
{
public:

    CHamitonianNJL1D(UINT uiSiteCount, Real fMass, Real fMu, Real fg, Real fLatticeSpacing = F(1.0));
    void SetLatticeSpacing(Real fa);
    void SetMass(Real fm);
    void SetChemicalPotential(Real fMu);
    void SetContactCoupling(Real fg);

protected:

    QLGate BuildSimulationCircuitABCBA(Real fTime, UINT uiTrotterStep) const override;

private:

    Real m_fG;
};


__END_NAMESPACE


#endif //#ifndef _CHAMITONIANLIST_H_

//=============================================================================
// END OF FILE
//=============================================================================