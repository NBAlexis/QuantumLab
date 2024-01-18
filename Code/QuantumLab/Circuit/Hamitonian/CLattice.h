//=============================================================================
// FILENAME : CLattice.h
// 
// DESCRIPTION:
// This is the class to map qubits into a lattice
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CLATTICE_H_
#define _CLATTICE_H_

__BEGIN_NAMESPACE

class QLAPI CLatticeSiteData
{
public:

    CLatticeSiteData() 
        : m_bCrossBoundary(FALSE) 
        , m_byStaggeredEta(0) 
    {

    }

    UBOOL operator==(const CLatticeSiteData& other) const
    {
        return m_lstSites == other.m_lstSites;
    }

    //The sites
    TArray<BYTE> m_lstSites;

    //When site is on links, whether it is inverse, note: FALSE = inverse
    TArray<UBOOL> m_lstDirectionOfLinks;

    //When site data is a neighbour pair, whether it cross boundary
    UBOOL m_bCrossBoundary;

    //If (eta & 1) = 1, it need to flip the coefficient for staggered fermion
    BYTE m_byStaggeredEta;
};

class QLAPI CLattice
{
public:

    virtual ~CLattice() {}
    
    virtual TArray<CLatticeSiteData> GetSites() const = 0;
    virtual TArray<CLatticeSiteData> GetNeighbourPairs() const = 0;
    virtual TArray<CLatticeSiteData> GetPlaquttes() const = 0;
    virtual TArray<CLatticeSiteData> GetCrosses() const = 0;
    virtual UINT GetControllerCount() const = 0;

    virtual void Draw() const = 0;

};


__END_NAMESPACE


#endif //#ifndef _CLATTICE_H_

//=============================================================================
// END OF FILE
//=============================================================================