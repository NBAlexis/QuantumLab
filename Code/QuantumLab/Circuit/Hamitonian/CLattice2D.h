//=============================================================================
// FILENAME : CLattice2D.h
// 
// DESCRIPTION:
// This build the data for lattice
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CLATTICE2D_H_
#define _CLATTICE2D_H_

__BEGIN_NAMESPACE

enum class EPeriodicType : INT
{
    EPT_Open,
    EPT_Torus,
    EPT_Klein,
    EPT_ProjectivePlane,
};

class QLAPI CLattice2D : public CLattice
{
public:

    CLattice2D(UINT uiX, UINT uiY, UBOOL bLinkOrSite = FALSE, EPeriodicType eBoundary = EPeriodicType::EPT_Open);

    TArray<CLatticeSiteData> GetSites() const override;
    TArray<CLatticeSiteData> GetNeighbourPairs() const override;


    TArray<CLatticeSiteData> GetPlaquttes() const override;
    TArray<CLatticeSiteData> GetCrosses() const override;

    void Draw() const override;

    UINT GetControllerCount() const override
    {
        return m_uiX * m_uiY;
    }

    UINT m_uiX;
    UINT m_uiY;
    UBOOL m_bLinkOrSite;
    EPeriodicType m_eBoundary;

protected:

    UBOOL GetCrossBoundaryX(UINT x, UINT y, UINT& site) const;
    UBOOL GetCrossBoundaryY(UINT x, UINT y, UINT& site) const;

    /**
    * dir = 0, or 1 means x or y
    */
    UINT GetLinkIndex(UINT x, UINT y, UINT dir) const
    {
        return (x * m_uiX + y) * 2 + dir;
    }

    /**
    * dir = -1, or -2 means minus x, or minus y
    * dir = 1 or 2 means plus x, or plus y
    * 
    * return value:
    * 1, -1: has this link, if is -1, it is reversed.
    * 0: not has this link
    */
    INT MapOutSideLink(INT x, INT y, INT dir, INT iReverse, UINT& linkIndex) const;

    static void PrintNumber(BYTE* grid, INT grididx, BYTE left, BYTE right, const CCString& sNumber);

    void DrawLinkXWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const;
    void DrawLinkYWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const;
    void DrawLinkInverseXWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const;
    void DrawLinkInverseYWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const;
    void DrawLinks(const TArray<CLatticeSiteData>& all) const;

    void DrawSitePrint(BYTE* grid, const CLatticeSiteData& link, UBOOL bXLink) const;
    void DrawSites(const TArray<CLatticeSiteData>& all) const;

    static void PrintGrids(BYTE* grid, UINT uiSize);


};


__END_NAMESPACE


#endif //#ifndef _CLATTICE2D_H_

//=============================================================================
// END OF FILE
//=============================================================================