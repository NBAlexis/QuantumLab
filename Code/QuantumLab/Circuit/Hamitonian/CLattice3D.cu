//=============================================================================
// FILENAME : CLattice2D.cu
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

CLattice2D::CLattice2D(UINT uiX, UINT uiY, UBOOL bLinkOrSite, EPeriodicType eBoundary)
    : CLattice()
    , m_uiX(uiX)
    , m_uiY(uiY)
    , m_bLinkOrSite(bLinkOrSite)
    , m_eBoundary(eBoundary)
{

}

TArray<CLatticeSiteData> CLattice2D::GetSites() const
{
    TArray<CLatticeSiteData> ret;
    for (UINT x = 0; x < m_uiX; ++x)
    {
        for (UINT y = 0; y < m_uiY; ++y)
        {
            CLatticeSiteData toAdd;
            toAdd.m_lstSites.AddItem(static_cast<BYTE>(x * m_uiY + y));
            toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
            toAdd.m_byStaggeredEta = static_cast<BYTE>((x + y) & 1);
            ret.AddItem(toAdd);
        }
    }
    return ret;
}

/**
* torus: 
* 
* (Lx-1, y) - (0,y)
* (x, Ly-1) - (x,0)
* 
* klein:
* 
* (Lx-1, y) - (0,Ly-y-1)
* (x, Ly-1) - (x,0)
* 
* projective plane:
* 
* (Lx-1, y) - (0,Ly-y-1)
* (x, Ly-1) - (Lx-x-1,0)
* 
*/
UBOOL CLattice2D::GetCrossBoundaryX(UINT x, UINT y, UINT& site) const
{
    if (x < m_uiX)
    {
        site = x * m_uiX + y;
        return FALSE;
    }

    switch (m_eBoundary)
    {
    case EPeriodicType::EPT_Open:
    case EPeriodicType::EPT_Torus:
        x = 0;
        site = x * m_uiX + y;
        return TRUE;
    case EPeriodicType::EPT_Klein:
    case EPeriodicType::EPT_ProjectivePlane:
        x = 0;
        y = m_uiY - y - 1;
        site = x * m_uiX + y;
        return TRUE;
    default:
        appWarning(_T("Not implemented!"));
        site = x * m_uiX + y;
        return FALSE;
    }
}

UBOOL CLattice2D::GetCrossBoundaryY(UINT x, UINT y, UINT& site) const
{
    if (y < m_uiY)
    {
        site = x * m_uiX + y;
        return FALSE;
    }

    switch (m_eBoundary)
    {
    case EPeriodicType::EPT_Open:
    case EPeriodicType::EPT_Torus:
    case EPeriodicType::EPT_Klein:
        y = 0;
        site = x * m_uiX + y;
        return TRUE;
    case EPeriodicType::EPT_ProjectivePlane:
        x = m_uiX - x - 1;
        y = 0;
        site = x * m_uiX + y;
        return TRUE;
    default:
        appWarning(_T("Not implemented!"));
        site = x * m_uiX + y;
        return FALSE;
    }
}

TArray<CLatticeSiteData> CLattice2D::GetNeighbourPairs() const
{
    TArray<CLatticeSiteData> ret;
    for (UINT x = 0; x < m_uiX; ++x)
    {
        for (UINT y = 0; y < m_uiY; ++y)
        {
            //Add X Link
            if (! (EPeriodicType::EPT_Open == m_eBoundary && x == m_uiX - 1))
            {
                UINT site1 = x * m_uiX + y;
                UINT site2 = 0;
                CLatticeSiteData toAdd;
                toAdd.m_bCrossBoundary = GetCrossBoundaryX(x + 1, y, site2);
                
                toAdd.m_lstSites.AddItem(static_cast<BYTE>(site1));
                toAdd.m_lstSites.AddItem(static_cast<BYTE>(site2));
                toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
                toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
                ret.AddItem(toAdd);
            }

            //Add Y Link
            if (!(EPeriodicType::EPT_Open == m_eBoundary && y == m_uiY - 1))
            {
                UINT site1 = x * m_uiX + y;
                UINT site2 = 0;
                CLatticeSiteData toAdd;
                toAdd.m_bCrossBoundary = GetCrossBoundaryY(x, y + 1, site2);

                toAdd.m_lstSites.AddItem(static_cast<BYTE>(site1));
                toAdd.m_lstSites.AddItem(static_cast<BYTE>(site2));
                toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
                toAdd.m_lstDirectionOfLinks.AddItem(TRUE);
                ret.AddItem(toAdd);
            }
        }
    }
    return ret;
}

INT CLattice2D::MapOutSideLink(INT x, INT y, INT dir, INT iReverse, UINT& linkIndex) const
{
    if (dir < 0)
    {
        if (-1 == dir)
        {
            x = x - 1;
            dir = -dir;
            iReverse = iReverse * -1;
        }
        if (-2 == dir)
        {
            y = y - 1;
            dir = -dir;
            iReverse = iReverse * -1;
        }
    }
    dir = dir - 1;

    if (x < 0 || x >= m_uiX)
    {
        switch (m_eBoundary)
        {
        case EPeriodicType::EPT_Open:
        case EPeriodicType::EPT_Torus:
            {
                x = (x < 0) ? (x + m_uiX) : (x - m_uiX);
            }
            break;
        case EPeriodicType::EPT_Klein:
        case EPeriodicType::EPT_ProjectivePlane:
            {
                x = (x < 0) ? (x + m_uiX) : (x - m_uiX);
                y = m_uiY - 1 - y;
                if (1 == dir)
                {
                    //y = y - 1;
                    iReverse = iReverse * -1;
                }
            }
            break;
        default:
            appWarning(_T("Not implemented!"));
            break;
        }
    }

    if (y < 0 || y >= m_uiY)
    {
        switch (m_eBoundary)
        {
        case EPeriodicType::EPT_Open:
        case EPeriodicType::EPT_Torus:
        case EPeriodicType::EPT_Klein:
            {
                y = (y < 0) ? (y + m_uiY) : (y - m_uiY);
            }
            break;
        case EPeriodicType::EPT_ProjectivePlane:
            {
                y = (y < 0) ? (y + m_uiY) : (y - m_uiY);
                x = m_uiX - 1 - x;
                if (0 == dir)
                {
                    //x = x - 1;
                    iReverse = iReverse * -1;
                }
            }
            break;
        default:
            appWarning(_T("Not implemented!"));
            break;
        }
    }

    linkIndex = GetLinkIndex(static_cast<UINT>(x), static_cast<UINT>(y), static_cast<UINT>(dir));

    if (0 == x && 1 == dir && EPeriodicType::EPT_Open == m_eBoundary)
    {
        return 0;
    }
    if (0 == y && 0 == dir && EPeriodicType::EPT_Open == m_eBoundary)
    {
        return 0;
    }
    return iReverse;
}

/**
* --->
* |  |
* <--|
*/
TArray<CLatticeSiteData> CLattice2D::GetPlaquttes() const
{
    TArray<CLatticeSiteData> ret;
    for (INT x = 0; x < static_cast<INT>(m_uiX); ++x)
    {
        for (INT y = 0; y < static_cast<INT>(m_uiY); ++y)
        {
            CLatticeSiteData onePlaq;
            UINT site = 0;

            INT ireverse = MapOutSideLink(x, y, 1, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ireverse = MapOutSideLink(x + 1, y, 2, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ireverse = MapOutSideLink(x + 1, y + 1, -1, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ireverse = MapOutSideLink(x, y + 1, -2, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ret.AddItem(onePlaq);
        }
    }
    return ret;
}

TArray<CLatticeSiteData> CLattice2D::GetCrosses() const
{
    TArray<CLatticeSiteData> ret;
    for (INT x = 0; x < static_cast<INT>(m_uiX); ++x)
    {
        for (INT y = 0; y < static_cast<INT>(m_uiY); ++y)
        {
            CLatticeSiteData onePlaq;
            UINT site = 0;

            INT ireverse = MapOutSideLink(x, y, 1, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ireverse = MapOutSideLink(x, y, 2, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ireverse = MapOutSideLink(x, y, -1, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ireverse = MapOutSideLink(x, y, -2, 1, site);
            onePlaq.m_lstSites.AddItem(site);
            onePlaq.m_lstDirectionOfLinks.AddItem(ireverse);

            ret.AddItem(onePlaq);
        }
    }
    return ret;
}

void CLattice2D::Draw() const
{
    if (m_bLinkOrSite)
    {
        DrawLinks(GetPlaquttes());
    }
    else
    {
        DrawSites(GetNeighbourPairs());
    }
}

void CLattice2D::PrintNumber(BYTE* grid, INT grididx, BYTE left, BYTE right, const CCString& sNumber)
{
    assert(sNumber.GetLength() > 0 && sNumber.GetLength() <= 3);

    if (3 == sNumber.GetLength())
    {
        grid[grididx] = sNumber.GetAt(0);
        grid[grididx + 1] = sNumber.GetAt(1);
        grid[grididx + 2] = sNumber.GetAt(2);
    }
    else if (2 == sNumber.GetLength())
    {
        grid[grididx] = left;
        grid[grididx + 1] = sNumber.GetAt(0);
        grid[grididx + 2] = sNumber.GetAt(1);
    }
    else
    {
        grid[grididx] = left;
        grid[grididx + 1] = sNumber.GetAt(0);
        grid[grididx + 2] = right;
    }
}

/**
* 321--123-->
*/
void CLattice2D::DrawLinkXWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const
{
    UINT siteIdx = idx >> 1;
    CCString sSiteIdx;
    sSiteIdx.Format(_T("%d"), siteIdx);
    CCString sLinkIdx;
    sLinkIdx.Format(_T("%d"), idx);

    BYTE byCode = (1 == iDir ? 252 : 246);

    UINT uiSiteNumber = y * 8 * (12 * m_uiX + 1) + x * 12;
    PrintNumber(grid, uiSiteNumber, 255, byCode, sSiteIdx);
    PrintNumber(grid, uiSiteNumber + 5, byCode, byCode, sLinkIdx);

    grid[uiSiteNumber + 3] = byCode;
    grid[uiSiteNumber + 4] = byCode;
    grid[uiSiteNumber + 8] = byCode;
    grid[uiSiteNumber + 9] = byCode;
    grid[uiSiteNumber + 10] = 251;

}

/**
*
* 321--123-->
*  ^       321
*  |        |
*  |       123
* 123       |
*  |        |
* 321       v
*  <--123--321
*
*/
void CLattice2D::DrawLinkYWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const
{
    UINT siteIdx = idx >> 1;
    CCString sSiteIdx;
    sSiteIdx.Format(_T("%d"), siteIdx);
    CCString sLinkIdx;
    sLinkIdx.Format(_T("%d"), idx);

    UINT uiSiteNumber = (y * 8 + 1) * (12 * m_uiX + 1) + x * 12 + 9;
    PrintNumber(grid, uiSiteNumber, 255, 255, sSiteIdx);
    uiSiteNumber = (y * 8 + 3) * (12 * m_uiX + 1) + x * 12 + 9;
    PrintNumber(grid, uiSiteNumber, 255, 255, sLinkIdx);

    grid[(y * 8 + 2) * (12 * m_uiX + 1) + x * 12 + 10] = (1 == iDir ? 253 : 247);
    grid[(y * 8 + 4) * (12 * m_uiX + 1) + x * 12 + 10] = (1 == iDir ? 253 : 247);
    grid[(y * 8 + 5) * (12 * m_uiX + 1) + x * 12 + 10] = (1 == iDir ? 253 : 247);
    grid[(y * 8 + 6) * (12 * m_uiX + 1) + x * 12 + 10] = 250;
}

/**
*
*  <--123--321
*
*/
void CLattice2D::DrawLinkInverseXWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const
{
    UINT siteIdx = idx >> 1;
    CCString sSiteIdx;
    sSiteIdx.Format(_T("%d"), siteIdx);
    CCString sLinkIdx;
    sLinkIdx.Format(_T("%d"), idx);

    BYTE byCode = (1 == iDir ? 252 : 246);

    UINT uiSiteNumber = (y * 8 + 7) * (12 * m_uiX + 1) + x * 12 + 9;
    PrintNumber(grid, uiSiteNumber, byCode, 255, sSiteIdx);
    PrintNumber(grid, uiSiteNumber - 5, byCode, byCode, sLinkIdx);

    grid[uiSiteNumber - 1] = byCode;
    grid[uiSiteNumber - 2] = byCode;
    grid[uiSiteNumber - 6] = byCode;
    grid[uiSiteNumber - 7] = byCode;
    grid[uiSiteNumber - 8] = 249;
}

/**
*
* 321--123-->
*  ^       321
*  |        |
*  |       123
* 123       |
*  |        |
* 321       v
*  <--123--321
*
*/
void CLattice2D::DrawLinkInverseYWithIndex(BYTE* grid, INT x, INT y, INT iDir, UINT idx) const
{
    UINT siteIdx = idx >> 1;
    CCString sSiteIdx;
    sSiteIdx.Format(_T("%d"), siteIdx);
    CCString sLinkIdx;
    sLinkIdx.Format(_T("%d"), idx);

    UINT uiSiteNumber = (y * 8 + 6) * (12 * m_uiX + 1) + x * 12;
    PrintNumber(grid, uiSiteNumber, 255, 255, sSiteIdx);
    uiSiteNumber = (y * 8 + 4) * (12 * m_uiX + 1) + x * 12;
    PrintNumber(grid, uiSiteNumber, 255, 255, sLinkIdx);

    grid[(y * 8 + 5) * (12 * m_uiX + 1) + x * 12 + 1] = (1 == iDir ? 253 : 247);
    grid[(y * 8 + 3) * (12 * m_uiX + 1) + x * 12 + 1] = (1 == iDir ? 253 : 247);
    grid[(y * 8 + 2) * (12 * m_uiX + 1) + x * 12 + 1] = (1 == iDir ? 253 : 247);
    grid[(y * 8 + 1) * (12 * m_uiX + 1) + x * 12 + 1] = 248;
}

/**
* 
* 321--123-->
*  ^       321
*  |        |
*  |       123
* 123       |
*  |        |
* 321       v
*  <--123--321
* 
*/
void CLattice2D::DrawLinks(const TArray<CLatticeSiteData>& all) const
{
    BYTE* grids = (BYTE*)malloc((8 * m_uiY) * (12 * m_uiX + 1) * sizeof(BYTE));
    for (UINT y = 0; y < 8 * m_uiY; ++y)
    {
        for (UINT x = 0; x <= 12 * m_uiX; ++x)
        {
            grids[x + y * (12 * m_uiX + 1)] = 255;
            if (x == 12 * m_uiX)
            {
                grids[x + y * (12 * m_uiX + 1)] = 254;
            }
        }
    }

    for (INT i = 0; i < all.Num(); ++i)
    {
        INT px = i / m_uiY;
        INT py = i % m_uiY;
        if (0 != all[i].m_lstDirectionOfLinks[0])
        {
            DrawLinkXWithIndex(grids, px, py, all[i].m_lstDirectionOfLinks[0], all[i].m_lstSites[0]);
        }

        if (0 != all[i].m_lstDirectionOfLinks[1])
        {
            DrawLinkYWithIndex(grids, px, py, all[i].m_lstDirectionOfLinks[1], all[i].m_lstSites[1]);
        }

        if (0 != all[i].m_lstDirectionOfLinks[2])
        {
            DrawLinkInverseXWithIndex(grids, px, py, all[i].m_lstDirectionOfLinks[2], all[i].m_lstSites[2]);
        }
        
        if (0 != all[i].m_lstDirectionOfLinks[3])
        {
            DrawLinkInverseYWithIndex(grids, px, py, all[i].m_lstDirectionOfLinks[3], all[i].m_lstSites[3]);
        }
    }

    PrintGrids(grids, (8 * m_uiY) * (12 * m_uiX + 1));
    appSafeFree(grids);
}

void CLattice2D::DrawSites(const TArray<CLatticeSiteData>& all) const
{
    BYTE* grids = (BYTE*)malloc((3 * m_uiY + 1) * (5 * m_uiX + 4) * sizeof(BYTE));
    for (UINT y = 0; y <= 3 * m_uiY; ++y)
    {
        for (UINT x = 0; x <= 5 * m_uiX + 3; ++x)
        {
            grids[x + y * (5 * m_uiX + 4)] = 255;
            if (x == 5 * m_uiX + 3)
            {
                grids[x + y * (5 * m_uiX + 4)] = 254;
            }
        }
    }

    UBOOL bXLink = TRUE;
    for (INT i = 0; i < all.Num(); ++i)
    {
        DrawSitePrint(grids, all[i], bXLink);
        bXLink = !bXLink;
    }

    PrintGrids(grids, (3 * m_uiY + 1) * (5 * m_uiX + 4));
    appSafeFree(grids);
}

void CLattice2D::DrawSitePrint(BYTE* grid, const CLatticeSiteData& link, UBOOL bXLink) const
{
    UINT uiX1 = link.m_lstSites[0] / m_uiY;
    UINT uiY1 = link.m_lstSites[0] % m_uiY;

    UINT uiX2 = link.m_lstSites[1] / m_uiY;
    UINT uiY2 = link.m_lstSites[1] % m_uiY;
    CCString sSiteNumber1 = _T("");
    sSiteNumber1.Format(_T("%d"), link.m_lstSites[0]);
    CCString sSiteNumber2 = _T("");
    sSiteNumber2.Format(_T("%d"), link.m_lstSites[1]);
    assert(sSiteNumber1.GetLength() >= 0 && sSiteNumber1.GetLength() <= 3);
    assert(sSiteNumber2.GetLength() >= 0 && sSiteNumber2.GetLength() <= 3);

    UBOOL bYLink = uiX1 == uiX2;
    if (EPeriodicType::EPT_ProjectivePlane == m_eBoundary)
    {
        bYLink = !bXLink;
    }

    //appGeneral(_T("x1:%d, y1:%d to x2:%d, y2:%d\n"), uiX1, uiY1, uiX2, uiY2);

    if (bYLink)
    {
        UINT y1 = 3 * uiY1;
        UINT y2 = 3 * uiY1 + 1;
        UINT y3 = 3 * uiY1 + 2;
        UINT y4 = 3 * uiY1 + 3;

        PrintNumber(grid, y1 * (5 * m_uiX + 4) + 5 * uiX1, 255, 255, sSiteNumber1);
        PrintNumber(grid, y4 * (5 * m_uiX + 4) + 5 * uiX1, 255, 255, sSiteNumber2);

        grid[y2 * (5 * m_uiX + 4) + 5 * uiX1 + 1] = 253;
        grid[y3 * (5 * m_uiX + 4) + 5 * uiX1 + 1] = 253;
    }
    else
    {
        UINT x1 = 5 * uiX1;

        UINT x4 = 5 * uiX1 + 3;
        UINT x5 = 5 * uiX1 + 4;

        UINT x6 = 5 * uiX1 + 5;

        PrintNumber(grid, 3 * uiY1 * (5 * m_uiX + 4) + x1, 255, 255, sSiteNumber1);
        PrintNumber(grid, 3 * uiY1 * (5 * m_uiX + 4) + x6, 255, 255, sSiteNumber2);

        grid[3 * uiY1 * (5 * m_uiX + 4) + x4] = 252;
        grid[3 * uiY1 * (5 * m_uiX + 4) + x5] = 252;
    }
}

void CLattice2D::PrintGrids(BYTE* grid, UINT uiSize)
{
    appPushLogDate(FALSE);
    appGeneral(_T("\n"));
    for (UINT i = 0; i < uiSize; ++i)
    {
        switch (grid[i])
        {
        case 255:
            appGeneral(_T(" "));
            break;
        case 254:
            appGeneral(_T("\n"));
            break;
        case 253:
            appGeneral(_T("│"));
            break;
        case 252:
            appGeneral(_T("─"));
            break;
        case 251:
            appGeneral(_T(">"));
            break;
        case 250:
            appGeneral(_T("▼"));
            break;
        case 249:
            appGeneral(_T("<"));
            break;
        case 248:
            appGeneral(_T("▲"));
            break;
        case 247:
            appGeneral(_T("┆"));
            break;
        case 246:
            appGeneral(_T("┄"));
            break;
        default:
            appGeneral(_T("%d"), grid[i] - '0');
            break;
        }
    }
    appGeneral(_T("\n"));
    appPopLogDate();
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================