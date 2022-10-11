//=============================================================================
// FILENAME : CSV.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [08/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLMatrix QLAPI ReadCSV(const CCString& fileName)
{
    std::ifstream file(fileName);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), fileName.c_str());
        return QLMatrix();
    }

    const SIZE_T buf_size = 1U << 20;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    INT leny = 0;
    INT lenx = 0;
    TArray<QLComplex> all;
    while (file.getline(buf, buf_size))
    {
        CCString sLine(buf);
        TArray<CCString> sep = appGetStringList(sLine, _T(','), EGSLF_IgnorTabSpace | EGSLF_IgnorTabSpaceInSide);

        if (0 != (sep.Num() & 1))
        {
            appCrucial(_T("CSV file format not good!\n"));
            file.close();
            return QLMatrix();
        }
        if (0 != lenx)
        {
            if (leny != sep.Num())
            {
                appCrucial(_T("CSV file format not good!\n"));
                file.close();
                return QLMatrix();
            }
        }

        if (0 == lenx)
        {
            leny = sep.Num() / 2;
        }

        for (INT i = 0; i < leny; ++i)
        {
#if _QL_DOUBLEFLOAT
            Real real = appStoD(sep[2 * i]);
            Real img = appStoD(sep[2 * i + 1]);
#else
            Real real = appStoF(sep[2 * i]);
            Real img = appStoF(sep[2 * i + 1]);
#endif
            all.AddItem(_make_cuComplex(real, img));
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }
    
    file.close();

    QLMatrix ret = QLMatrix::CopyCreate(static_cast<UINT>(lenx), static_cast<UINT>(leny), all.GetData());
    ret.Transpose();

    return ret;
}

void QLAPI SaveCSV(const QLMatrix& m, const CCString& fileName)
{
    std::ofstream file(fileName);
    CCString sOut;

    for (UINT i = 0; i < m.Y(); ++i)
    {
        for (UINT j = 0; j < m.X(); ++j)
        {
            QLComplex c = m.Get(j, i);
            sOut.Format(_T("%f, %f"), c.x, c.y);
            file << sOut.c_str();
            if (j != m.X() - 1)
            {
                file << _T(",");
            }
            else
            {
                file << std::endl;
            }
        }
    }
}

QLMatrix QLAPI ReadCSVR(const CCString& fileName)
{
    std::ifstream file(fileName);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), fileName.c_str());
        return QLMatrix();
    }

    const SIZE_T buf_size = 1U << 20;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    INT leny = 0;
    INT lenx = 0;
    TArray<QLComplex> all;
    while (file.getline(buf, buf_size))
    {
        CCString sLine(buf);
        TArray<CCString> sep = appGetStringList(sLine, _T(','), EGSLF_IgnorTabSpace | EGSLF_IgnorTabSpaceInSide);

        if (0 != lenx)
        {
            if (leny != sep.Num())
            {
                appCrucial(_T("CSV file format not good!\n"));
                file.close();
                return QLMatrix();
            }
        }

        if (0 == lenx)
        {
            leny = sep.Num();
        }

        for (INT i = 0; i < leny; ++i)
        {
#if _QL_DOUBLEFLOAT
            Real real = appStoD(sep[i]);
#else
            Real real = appStoF(sep[i]);
#endif
            all.AddItem(_make_cuComplex(real, F(0.0)));
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }

    file.close();

    QLMatrix ret = QLMatrix::CopyCreate(static_cast<UINT>(lenx), static_cast<UINT>(leny), all.GetData());
    ret.Transpose();

    return ret;
}

void QLAPI SaveCSVR(const QLMatrix& m, const CCString& fileName)
{
    std::ofstream file(fileName);
    CCString sOut;

    for (UINT i = 0; i < m.Y(); ++i)
    {
        for (UINT j = 0; j < m.X(); ++j)
        {
            QLComplex c = m.Get(j, i);
            sOut.Format(_T("%f"), c.x);
            file << sOut.c_str();
            if (j != m.X() - 1)
            {
                file << _T(",");
            }
            else
            {
                file << std::endl;
            }
        }
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================