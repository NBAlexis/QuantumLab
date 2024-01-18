//=============================================================================
// FILENAME : CHamitonianStaggeredJordanWigner1D.h
// 
// DESCRIPTION:
// This is for simulation of staggered fermions in 1D
// The fermions are decomposed as Jordan Wigner strings
// 
// Note: only support even site numbers because of staggered fermion
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CHAMITONIANSTAGGEREDJORDANWIGNER1D_H_
#define _CHAMITONIANSTAGGEREDJORDANWIGNER1D_H_

__BEGIN_NAMESPACE

/**
* \bar{psi} i \gamma _1 \partial _1 \psi
* 
* coefficient = 1/a
*/
class QLAPI CHamitonianStaggeredJordanWigner1DKinetic : public CHamitonianTerm
{
public:

    CHamitonianStaggeredJordanWigner1DKinetic(Real fCoeff, EPauliType eType)
        : CHamitonianTerm(fCoeff)
        , m_ePauliType(eType)

    {

    }

    TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const override;

private:

    //Fully term is two uncommutable terms, one with sigmaX and the other with sigmaY
    EPauliType m_ePauliType;
};

/**
* \bar{psi} \psi
* 
* coefficient = m
*/
class QLAPI CHamitonianStaggeredJordanWigner1DPsibarPsi : public CHamitonianTerm
{
public:

    CHamitonianStaggeredJordanWigner1DPsibarPsi(Real fCoeff) : CHamitonianTerm(fCoeff)
    {

    }

    TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const override;
};

/**
* \bar{psi} \gamma _0 \psi
* 
* coefficient = -mu
*/
class QLAPI CHamitonianStaggeredJordanWigner1DPsibarG0Psi : public CHamitonianTerm
{
public:

    CHamitonianStaggeredJordanWigner1DPsibarG0Psi(Real fCoeff) : CHamitonianTerm(fCoeff)
    {

    }

    TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const override;
};

/**
* (\bar{psi} \psi)^2
* 
* coefficient = g/a
*/
class QLAPI CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure : public CHamitonianTerm
{
public:

    CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure(Real fCoeff) : CHamitonianTerm(fCoeff)
    {

    }

    TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const override;
};


__END_NAMESPACE


#endif //#ifndef _CHAMITONIANSTAGGEREDJORDANWIGNER1D_H_

//=============================================================================
// END OF FILE
//=============================================================================