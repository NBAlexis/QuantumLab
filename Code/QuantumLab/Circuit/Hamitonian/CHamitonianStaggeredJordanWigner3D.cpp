//=============================================================================
// FILENAME : CHamitonianStaggeredJordanWigner1D.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

/**
* (1/4a) n=0,...,n-2 [sx(n)sx(n+1)+sy(n)sy(n+1)]
* + (-1)^{N/2} /4a [sx(N-1)sx(0)+sy(N-1)sy(0)] prod _{n=0,...,n-2} sz(n)
*/
TArray<PauliProduct> CHamitonianStaggeredJordanWigner1DKinetic::GetAllTerms(const CLattice* pLattice) const
{
    TArray<CLatticeSiteData> pairs = pLattice->GetNeighbourPairs();
    UINT uiControllerCount = pLattice->GetControllerCount();

    TArray<PauliProduct> ret;

    UBOOL bReverse = (uiControllerCount >> 1) & 1;

    Real fCoeff = m_fCoefficient * F(0.25);

    for (INT i = 0; i < pairs.Num(); ++i)
    {
        TArray<BYTE> pauliProduct;
        for (UINT uiQubit = 0; uiQubit < uiControllerCount; ++uiQubit)
        {
            if (pairs[i].m_bCrossBoundary)
            {
                if (static_cast<UINT>(pairs[i].m_lstSites[0]) == uiQubit)
                {
                    pauliProduct.AddItem(static_cast<BYTE>(m_ePauliType));
                }
                else if (static_cast<UINT>(pairs[i].m_lstSites[1]) == uiQubit)
                {
                    pauliProduct.AddItem(static_cast<BYTE>(m_ePauliType));
                }
                else
                {
                    pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_Z));
                }
            }
            else
            {
                if (static_cast<UINT>(pairs[i].m_lstSites[0]) == uiQubit)
                {
                    pauliProduct.AddItem(static_cast<BYTE>(m_ePauliType));
                }
                else if (static_cast<UINT>(pairs[i].m_lstSites[1]) == uiQubit)
                {
                    pauliProduct.AddItem(static_cast<BYTE>(m_ePauliType));
                }
                else
                {
                    pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_I));
                }
            }
        }

        if (pairs[i].m_bCrossBoundary && bReverse)
        {
            ret.AddItem(PauliProduct(pauliProduct, -fCoeff));
        }
        else
        {
            ret.AddItem(PauliProduct(pauliProduct, fCoeff));
        }
    }

    return ret;
}

/**
* (-1)^n sigma_z / 2
*/
TArray<PauliProduct> CHamitonianStaggeredJordanWigner1DPsibarPsi::GetAllTerms(const CLattice* pLattice) const
{
    TArray<CLatticeSiteData> sites = pLattice->GetSites();
    UINT uiControllerCount = pLattice->GetControllerCount();

    Real fCoeff = m_fCoefficient * F(0.5);

    TArray<PauliProduct> ret;

    for (INT i = 0; i < sites.Num(); ++i)
    {
        TArray<BYTE> pauliProduct;
        for (UINT uiQubit = 0; uiQubit < uiControllerCount; ++uiQubit)
        {
            if (static_cast<UINT>(sites[i].m_lstSites[0]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_Z));
            }
            else
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_I));
            }
        }

        ret.AddItem(PauliProduct(pauliProduct, (sites[i].m_byStaggeredEta & 1) ? -fCoeff : fCoeff));
    }

    return ret;
}

/**
* sigma_z / 2
*/
TArray<PauliProduct> CHamitonianStaggeredJordanWigner1DPsibarG0Psi::GetAllTerms(const CLattice* pLattice) const
{
    TArray<CLatticeSiteData> sites = pLattice->GetSites();
    UINT uiControllerCount = pLattice->GetControllerCount();

    Real fCoeff = m_fCoefficient * F(0.5);

    TArray<PauliProduct> ret;

    for (INT i = 0; i < sites.Num(); ++i)
    {
        TArray<BYTE> pauliProduct;
        for (UINT uiQubit = 0; uiQubit < uiControllerCount; ++uiQubit)
        {
            if (static_cast<UINT>(sites[i].m_lstSites[0]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_Z));
            }
            else
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_I));
            }
        }

        ret.AddItem(PauliProduct(pauliProduct, fCoeff));
    }

    return ret;
}

/**
*   -g/4a[\sum _{n=0,N-2}(1+sz(n))(1+sz(n+1))  + (1+sz(N-1))(1+sz(0)) - 2\sum _{n=0,N-1}(1+sz(n))]
* = -g/4a[\sum _{n=0,N-2}(1+sz(n)+sz(n+1)+sz(n)sz(n+1))  + (1+sz(N-1)+sz(0)+sz(0)sz(N-1)) - 2\sum _{n=0,N-1}(1+sz(n))]
* = -g/4a[\sum _{n=0,N-2}(1+sz(n) + 1+sz(n+1) +sz(n)sz(n+1) - 1)  + (1+sz(N-1) +1+sz(0) +sz(0)sz(N-1) - 1) - 2\sum _{n=0,N-1}(1+sz(n))]
* drop constants
* = -g/4a[\sum _{n=0,N-2}(1+sz(n) + 1+sz(n+1) +sz(n)sz(n+1))  + (1+sz(N-1) +1+sz(0) +sz(0)sz(N-1)) - 2\sum _{n=0,N-1}(1+sz(n))]
* = -g/4a[\sum _{n=0,N-2}sz(n)sz(n+1)  + sz(0)sz(N-1)]
*/
TArray<PauliProduct> CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure::GetAllTerms(const CLattice* pLattice) const
{
    Real fCoeff = m_fCoefficient * F(-0.25);

    TArray<CLatticeSiteData> pairs = pLattice->GetNeighbourPairs();
    UINT uiControllerCount = pLattice->GetControllerCount();

    TArray<PauliProduct> ret;

    for (INT i = 0; i < pairs.Num(); ++i)
    {
        TArray<BYTE> pauliProduct;
        for (UINT uiQubit = 0; uiQubit < uiControllerCount; ++uiQubit)
        {
            if (static_cast<UINT>(pairs[i].m_lstSites[0]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_Z));
            }
            else if (static_cast<UINT>(pairs[i].m_lstSites[1]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_Z));
            }
            else
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_I));
            }
        }

        ret.AddItem(PauliProduct(pauliProduct, fCoeff));
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================