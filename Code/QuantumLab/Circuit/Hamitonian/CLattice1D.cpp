//=============================================================================
// FILENAME : CLattice1D.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [05/12/2023 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE


CLattice1D::CLattice1D(UINT numberOfSites, UBOOL bPeriodic)
    : m_uiNumberOfSites(numberOfSites)
    , m_bPeriodic(bPeriodic)
{

}

TArray<CLatticeSiteData> CLattice1D::GetSites() const
{
    TArray<CLatticeSiteData> ret;
    for (UINT i = 0; i < m_uiNumberOfSites; ++i)
    {
        CLatticeSiteData toAdd;
        toAdd.m_lstSites.AddItem(static_cast<BYTE>(i));
        toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
        toAdd.m_byStaggeredEta = static_cast<BYTE>(i & 1);
        ret.AddItem(toAdd);
    }
    return ret;
}

TArray<CLatticeSiteData> CLattice1D::GetNeighbourPairs() const
{
    TArray<CLatticeSiteData> ret;
    for (UINT i = 0; i < m_uiNumberOfSites; ++i)
    {
        CLatticeSiteData toAdd;
        if (m_uiNumberOfSites - 1 == i)
        {
            if (m_bPeriodic)
            {
                toAdd.m_lstSites.AddItem(static_cast<BYTE>(i));
                toAdd.m_lstSites.AddItem(static_cast<BYTE>(0));
                toAdd.m_bCrossBoundary = TRUE;
            }
        }
        else
        {
            toAdd.m_lstSites.AddItem(static_cast<BYTE>(i));
            toAdd.m_lstSites.AddItem(static_cast<BYTE>(i + 1));
        }
        toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
        toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
        ret.AddItem(toAdd);
    }
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================