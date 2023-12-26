//=============================================================================
// FILENAME : QLSimulator.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

TArray<PauliProduct> CHamitonianPauli::GetAllTerms(const CLattice* pLattice) const
{
    TArray<PauliProduct> ret;

    TArray<CLatticeSiteData> sites = pLattice->GetSites();
    UINT uiControllerCount = pLattice->GetControllerCount();

    for (INT i = 0; i < sites.Num(); ++i)
    {
        TArray<BYTE> pauliProduct;
        for (UINT uiQubit = 0; uiQubit < uiControllerCount; ++uiQubit)
        {
            if (static_cast<UINT>(sites[i].m_lstSites[0]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(m_eType));
            }
            else
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_I));
            }
        }

        ret.AddItem(PauliProduct(pauliProduct, m_fCoefficient));
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================