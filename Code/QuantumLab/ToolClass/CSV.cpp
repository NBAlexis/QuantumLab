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

        //if (0 != (sep.Num() & 1))
        //{
        //    appCrucial(_T("CSV file format not good!\n"));
        //    file.close();
        //    return QLMatrix();
        //}

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

        INT halfleny = (leny + 1) / 2;
        for (INT i = 0; i < halfleny; ++i)
        {
#if _QL_DOUBLEFLOAT
            Real real = appStoD(sep[2 * i]);
            Real img = ((2 * i + 1) < leny) ? appStoD(sep[2 * i + 1]) : 0.0;
#else
            Real real = appStoF(sep[2 * i]);
            Real img = ((2 * i + 1) < leny) ? appStoF(sep[2 * i + 1]) : 0.0f;
#endif
            all.AddItem(_make_cuComplex(real, img));
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }
    
    file.close();

    QLMatrix ret = QLMatrix::CopyCreate(static_cast<UINT>(lenx), static_cast<UINT>((leny + 1) / 2), all.GetData());
    ret.Transpose();

    return ret;
}

TArray<QLComplex> QLAPI ReadCSVA(const CCString& fileName, UINT &w, UINT &h)
{
    TArray<QLComplex> all;
    std::ifstream file(fileName);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), fileName.c_str());
        w = 0;
        h = 0;
        return all;
    }

    const SIZE_T buf_size = 1U << 20;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    INT leny = 0;
    INT lenx = 0;
    
    while (file.getline(buf, buf_size))
    {
        CCString sLine(buf);
        TArray<CCString> sep = appGetStringList(sLine, _T(','), EGSLF_IgnorTabSpace | EGSLF_IgnorTabSpaceInSide);

        //if (0 != (sep.Num() & 1))
        //{
        //    appCrucial(_T("CSV file format not good!\n"));
        //    file.close();
        //    w = leny;
        //    h = lenx;
        //    return all;
        //}
        if (0 != lenx)
        {
            if (leny != sep.Num())
            {
                appCrucial(_T("CSV file format not good!\n"));
                file.close();
                w = (leny + 1) / 2;
                h = lenx;
                return all;
            }
        }

        if (0 == lenx)
        {
            leny = sep.Num();
        }
        INT halfleny = (leny + 1) / 2;
        for (INT i = 0; i < halfleny; ++i)
        {
#if _QL_DOUBLEFLOAT
            Real real = appStoD(sep[2 * i]);
            Real img = ((2 * i + 1) < leny) ? appStoD(sep[2 * i + 1]) : 0.0;
#else
            Real real = appStoF(sep[2 * i]);
            Real img = ((2 * i + 1) < leny) ? appStoF(sep[2 * i + 1]) : 0.0f;
#endif
            all.AddItem(_make_cuComplex(real, img));
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }

    file.close();

    w = (leny + 1) / 2;
    h = lenx;
    return all;
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

    file.flush();
    file.close();
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

    file.flush();
    file.close();
}

TArray<Real> QLAPI ReadCSVAR(const CCString& fileName, UINT& w, UINT& h)
{
    std::ifstream file(fileName);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), fileName.c_str());
        w = 0;
        h = 0;
        return TArray<Real>();
    }

    const SIZE_T buf_size = 1U << 20;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    INT leny = 0;
    INT lenx = 0;
    TArray<Real> all;
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
                w = 0;
                h = 0;
                return TArray<Real>();
            }
        }

        if (0 == lenx)
        {
            leny = sep.Num();
        }

        TArray<Real> oneLine;
        for (INT i = 0; i < leny; ++i)
        {
#if _QL_DOUBLEFLOAT
            Real real = appStoD(sep[i]);
#else
            Real real = appStoF(sep[i]);
#endif
            all.AddItem(real);
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }

    file.close();

    w = static_cast<UINT>(leny);
    h = static_cast<UINT>(lenx);

    return all;
}

TArray<FLOAT> QLAPI ReadCSVAF(const CCString& fileName, UINT& w, UINT& h)
{
    std::ifstream file(fileName);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), fileName.c_str());
        w = 0;
        h = 0;
        return TArray<FLOAT>();
    }

    const SIZE_T buf_size = 1U << 20;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    INT leny = 0;
    INT lenx = 0;
    TArray<FLOAT> all;
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
                w = 0;
                h = 0;
                return TArray<FLOAT>();
            }
        }

        if (0 == lenx)
        {
            leny = sep.Num();
        }

        TArray<FLOAT> oneLine;
        for (INT i = 0; i < leny; ++i)
        {
            FLOAT real = appStoF(sep[i]);
            all.AddItem(real);
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }

    file.close();

    w = static_cast<UINT>(leny);
    h = static_cast<UINT>(lenx);

    return all;
}

TArray<INT> QLAPI ReadCSVAI(const CCString& fileName, UINT& w, UINT& h)
{
    std::ifstream file(fileName);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), fileName.c_str());
        w = 0;
        h = 0;
        return TArray<INT>();
    }

    const SIZE_T buf_size = 1U << 20;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    INT leny = 0;
    INT lenx = 0;
    TArray<INT> all;
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
                w = 0;
                h = 0;
                return TArray<INT>();
            }
        }

        if (0 == lenx)
        {
            leny = sep.Num();
        }

        TArray<INT> oneLine;
        for (INT i = 0; i < leny; ++i)
        {
            INT real = appStoI(sep[i]);
            all.AddItem(real);
        }

        ++lenx;
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }

    file.close();

    w = static_cast<UINT>(leny);
    h = static_cast<UINT>(lenx);

    return all;
}

void QLAPI SaveCSVAR(const Real* m, UINT w, UINT h, const CCString& fileName)
{
    std::ofstream file(fileName);
    CCString sOut;

    for (UINT i = 0; i < h; ++i)
    {
        for (UINT j = 0; j < w; ++j)
        {
            sOut.Format(_T("%f"), m[i * w + j]);
            file << sOut.c_str();
            if (j != w - 1)
            {
                file << _T(",");
            }
            else
            {
                file << std::endl;
            }
        }
    }

    file.flush();
    file.close();
}

template<class T>
void _SaveCSVAI(const T* m, UINT w, UINT h, const CCString& fileName)
{
    std::ofstream file(fileName);
    CCString sOut;

    for (UINT i = 0; i < h; ++i)
    {
        for (UINT j = 0; j < w; ++j)
        {
            sOut.Format(_T("%d"), m[i * w + j]);
            file << sOut.c_str();
            if (j != w - 1)
            {
                file << _T(",");
            }
            else
            {
                file << std::endl;
            }
        }
    }

    file.flush();
    file.close();
}

void QLAPI SaveCSVAI(const INT* m, UINT w, UINT h, const CCString& fileName)
{
    _SaveCSVAI(m, w, h, fileName);
}

void QLAPI SaveCSVAUI(const UINT* m, UINT w, UINT h, const CCString& fileName)
{
    _SaveCSVAI(m, w, h, fileName);
}

void QLAPI SaveCSVAB(const BYTE* m, UINT w, UINT h, const CCString& fileName)
{
    _SaveCSVAI(m, w, h, fileName);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================