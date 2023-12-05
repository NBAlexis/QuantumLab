//=============================================================================
// FILENAME : CLattice1D.h
// 
// DESCRIPTION:
// This build the data for lattice
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CLATTICE1D_H_
#define _CLATTICE1D_H_

__BEGIN_NAMESPACE

class QLAPI CLattice1D : public CLattice
{
public:

    CLattice1D(UINT numberOfSites, UBOOL bPeriodic = TRUE);

    TArray<CLatticeSiteData> GetSites() const override;
    TArray<CLatticeSiteData> GetNeighbourPairs() const override;


    TArray<CLatticeSiteData> GetPlaquttes() const override
    {
        appCrucial(_T("1D Lattice does NOT support plaqutte!\n"));
        TArray<CLatticeSiteData> ret;
        return ret;
    }

    void Draw() const
    {
        appWarning(_T("1D Lattice is naive to draw.\n"));
    }

    UINT GetControllerCount() const 
    {
        return m_uiNumberOfSites;
    }

    UINT m_uiNumberOfSites;
    UBOOL m_bPeriodic;

};

__END_NAMESPACE


#endif //#ifndef _CLATTICE1D_H_

//=============================================================================
// END OF FILE
//=============================================================================