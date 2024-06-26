//=============================================================================
// FILENAME : CCString.h
// 
// DESCRIPTION:
//  Use this string class instead of STL, to get rid of the warnings with dll-export
//
// REVISION: [dd/mm/yy]
//  [13/09/2022 nbale]
//=============================================================================
#ifndef _CCSTRING_H_
#define _CCSTRING_H_

#define defineToStr(type) \
template<> \
inline CCString appToString(const type& content) \
{ \
    return CCString(std::to_string(content).c_str()); \
}

__BEGIN_NAMESPACE

inline void appStrTrimLeft(TCHAR*& InStr)
{
    while (*InStr && (' ' == (*InStr) || '\t' == (*InStr)))
        ++InStr;
}

inline INT appStrToINT(const TCHAR* s)
{
    INT base;
    TCHAR* p = const_cast<TCHAR*>(s);
    appStrTrimLeft(p);
    if ('0' == p[0] && ('x' == p[1] || 'X' == p[1]))
        base = 16;
    else
        base = 10;
    return appStoI(p, base);
}

inline UINT appStrToUINT(const TCHAR* s)
{
    INT base;
    TCHAR* p = const_cast<TCHAR*>(s);
    appStrTrimLeft(p);
    if ('0' == p[0] && ('x' == p[1] || 'X' == p[1]))
        base = 16;
    else
        base = 10;
    return appStoUI(p, base);
}

inline FLOAT appStrToFLOAT(const TCHAR* s)
{
    return appStoF(s);
}

inline DOUBLE appStrToDOUBLE(const TCHAR* s)
{
    return appStoD(s);
}

inline Real appStrToReal(const TCHAR* s)
{
    return static_cast<Real>(appStoD(s));
}

//=====================================================
//    CCStringData
//=====================================================
#define __SafeStrlen(s)        ((INT)((s)?appStrlen(s):0))

class CCStringData
{
private:
    std::atomic<INT> m_nRefs;             // reference count

public:

    friend class CCString;

    INT m_nDataLength;        // length of data (including terminator)
    INT m_nAllocLength;       // length of allocation

    TCHAR* Data()           // TCHAR* to managed data
    {
        //skip myself
        return (TCHAR*)(this + 1);
    }
};

//=====================================================
//    CCString
//=====================================================
class QLAPI CCString
{
public:
    enum { MAX_PATH = 8192, };

    CCString();
    CCString(const CCString& stringSrc);
    CCString(TCHAR ch, INT nRepeat = 1);
    CCString(const TCHAR* lpsz);
    CCString(const TCHAR* lpch, INT nLength);
    ~CCString();

    CCStringData* GetData() const
    {
        assert(m_pchData != NULL);
        //skip myself
        return (CCStringData*)(m_pchData)-1;
    }

protected:

    // implementation helpers
    void Init();
    void Release();
    static void Release(CCStringData* pData);
    static void FreeData(CCStringData* pData)
    {
        delete[](BYTE*)(pData);
    }
    /**
    * always allocate one extra character for '\0' termination
    * assumes [optimistically] that data length will equal allocation length
    * pData = [...CCStringData...][...RealString...][\0]
    */
    void AllocBuffer(INT nLen)
    {
        assert(nLen >= 0);
        assert(nLen <= INT_MAX - 1);    // max size (enough room for 1 extra)

        if (nLen == 0)
            Init();
        else
        {
            CCStringData* pData;
            pData = (CCStringData*)(new BYTE[sizeof(CCStringData) + (nLen + 1) * sizeof(TCHAR)]);
            pData->m_nAllocLength = nLen;
            pData->m_nRefs = 1;
            pData->Data()[nLen] = _T('\0');
            pData->m_nDataLength = nLen;
            m_pchData = pData->Data();
        }
    }
    void AllocBeforeWrite(INT nLen)
    {
        if (GetData()->m_nRefs > 1 || nLen > GetData()->m_nAllocLength)
        {
            Release();
            AllocBuffer(nLen);
        }
        assert(GetData()->m_nRefs <= 1);
    }
    /**
    * will clone the data attached to this string
    * allocating 'nExtraLen' characters
    * Places results in uninitialized string 'dest'
    * Will copy the part or all of original data to start of new string
    */
    void AllocCopy(CCString& dest, INT nCopyLen, INT nCopyIndex, INT nExtraLen) const
    {
        const INT nNewLen = nCopyLen + nExtraLen;
        if (nNewLen == 0)
        {
            dest.Init();
        }
        else
        {
            dest.AllocBuffer(nNewLen);
            memcpy(dest.m_pchData, m_pchData + nCopyIndex, nCopyLen * sizeof(TCHAR));
        }
    }
    void AssignCopy(INT nSrcLen, const TCHAR* lpszSrcData)
    {
        AllocBeforeWrite(nSrcLen);
        memcpy(m_pchData, lpszSrcData, nSrcLen * sizeof(TCHAR));
        GetData()->m_nDataLength = nSrcLen;
        m_pchData[nSrcLen] = _T('\0');
    }

    /**
    * The buffer is going to change.
    * There might be multiple CCString instance pointing to a same string buffer
    * Therefor, copy it, so it will not affect other strings.
    * For example, a="X", b=a;
    * when doing a.MakeLower(),
    * if not copy the buffer, a and b will become "x" at the same time
    */
    void CopyBeforeWrite()
    {
        if (GetData()->m_nRefs > 1)
        {
            CCStringData* pData = GetData();
            Release();
            AllocBuffer(pData->m_nDataLength);
            memcpy(m_pchData, pData->Data(), (pData->m_nDataLength + 1) * sizeof(TCHAR));
        }
        assert(GetData()->m_nRefs <= 1);
    }
    /**
    * -- master concatenation routine
    * Concatenate two sources
    * -- assume that 'this' is a new CCStringData object
    */
    void ConcatCopy(INT nSrc1Len, const TCHAR* lpszSrc1Data, INT nSrc2Len, const TCHAR* lpszSrc2Data)
    {
        const INT nNewLen = nSrc1Len + nSrc2Len;
        if (nNewLen != 0)
        {
            AllocBuffer(nNewLen);
            memcpy(m_pchData, lpszSrc1Data, nSrc1Len * sizeof(TCHAR));
            memcpy(m_pchData + nSrc1Len, lpszSrc2Data, nSrc2Len * sizeof(TCHAR));
        }
    }
    /**
    *  -- the main routine for += operators
    */
    void ConcatInPlace(INT nSrcLen, const TCHAR* lpszSrcData)
    {
        // concatenating an empty string is a no-op!
        if (nSrcLen == 0)
            return;

        // if the buffer is too small, or we have a width mis-match, just
        //   allocate a new buffer (slow but sure)
        if (GetData()->m_nRefs > 1 || GetData()->m_nDataLength + nSrcLen > GetData()->m_nAllocLength)
        {
            // we have to grow the buffer, use the ConcatCopy routine
            CCStringData* pOldData = GetData();
            ConcatCopy(GetData()->m_nDataLength, m_pchData, nSrcLen, lpszSrcData);
            assert(pOldData != NULL);
            CCString::Release(pOldData);
        }
        else
        {
            // fast concatenation when buffer big enough
            memcpy(m_pchData + GetData()->m_nDataLength, lpszSrcData, nSrcLen * sizeof(TCHAR));
            GetData()->m_nDataLength += nSrcLen;
            assert(GetData()->m_nDataLength <= GetData()->m_nAllocLength);
            m_pchData[GetData()->m_nDataLength] = _T('\0');
        }
    }

public:

    // Attributes & Operations
    INT GetLength() const { return GetData()->m_nDataLength; }
    INT GetAllocLength() const { return GetData()->m_nAllocLength; }
    UBOOL IsEmpty() const { return 0 == GetData()->m_nDataLength; }
    void Empty(); //out side

    TCHAR GetAt(INT nIndex) const
    {
        assert(nIndex >= 0);
        assert(nIndex < GetData()->m_nDataLength);
        return m_pchData[nIndex];
    }

    TCHAR operator[](INT nIndex) const { return GetAt(nIndex); }
    operator const TCHAR* () const { return m_pchData; }

    const TCHAR* c_str() const { return m_pchData; }

#if _CLG_UNICODE
    //attention here!
    //need to free it yourself!
    operator const ANSICHAR* () const
    {
        INT buf_len = GetLength() * 2 + 1;
        ANSICHAR* new_data = (ANSICHAR*)(appMalloc(buf_len));
        INT ret_len = appUnicodeToAnsi(new_data, (UNICHAR*)m_pchData/*with terminator*/, buf_len/*out buffer bytes*/);
        assert(ret_len <= buf_len);
        return (const ANSICHAR*)new_data;
    }
#endif

    void SetAt(INT nIndex, TCHAR ch)
    {
        assert(nIndex >= 0);
        assert(nIndex < GetData()->m_nDataLength);
        CopyBeforeWrite();
        m_pchData[nIndex] = ch;
    }

    // overloaded assignment
    const CCString& operator=(const CCString& stringSrc); //outside
    const CCString& operator=(TCHAR ch)
    {
        AssignCopy(1, &ch);
        return *this;
    }
    const CCString& operator=(const TCHAR* lpsz)
    {
        AssignCopy(__SafeStrlen(lpsz), lpsz);
        return *this;
    }

    // concatenation
    const CCString& operator+=(const CCString& string)
    {
        ConcatInPlace(string.GetData()->m_nDataLength, string.m_pchData);
        return *this;
    }
    const CCString& operator+=(TCHAR ch)
    {
        ConcatInPlace(1, &ch);
        return *this;
    }
    const CCString& operator+=(const TCHAR* lpsz)
    {
        ConcatInPlace(__SafeStrlen(lpsz), lpsz);
        return *this;
    }

    QLAPI friend CCString operator+(const CCString& string1, const CCString& string2);
    QLAPI friend CCString operator+(const CCString& string, TCHAR ch);
    QLAPI friend CCString operator+(TCHAR ch, const CCString& string);
    QLAPI friend CCString operator+(const CCString& string, const TCHAR* lpsz);
    QLAPI friend CCString operator+(const TCHAR* lpsz, const CCString& string);

    // string comparison
    INT Compare(const TCHAR* lpsz) const { return appStrcmp(m_pchData, lpsz); }
    INT CompareNoCase(const TCHAR* lpsz) const { return appStricmp(m_pchData, lpsz); }

    // simple sub-string extraction
    CCString Mid(INT nFirst, INT nCount) const
    {
        // out-of-bounds requests return sensible things
        if (nFirst < 0)
            nFirst = 0;
        if (nCount < 0)
            nCount = 0;

        if (nFirst + nCount > GetData()->m_nDataLength)
            nCount = GetData()->m_nDataLength - nFirst;
        if (nFirst > GetData()->m_nDataLength)
            nCount = 0;

        assert(nFirst >= 0);
        assert(nFirst + nCount <= GetData()->m_nDataLength);

        // optimize case of returning entire string
        if (nFirst == 0 && nFirst + nCount == GetData()->m_nDataLength)
            return *this;

        CCString dest;
        AllocCopy(dest, nCount, nFirst, 0);
        return dest;
    }
    CCString Mid(INT nFirst) const { return Mid(nFirst, GetData()->m_nDataLength - nFirst); }
    CCString Left(INT nCount) const
    {
        if (nCount < 0)
            nCount = 0;
        if (nCount >= GetData()->m_nDataLength)
            return *this;

        CCString dest;
        AllocCopy(dest, nCount, 0, 0);
        return dest;
    }
    CCString Right(INT nCount) const
    {
        if (nCount < 0)
            nCount = 0;
        if (nCount >= GetData()->m_nDataLength)
            return *this;

        CCString dest;
        AllocCopy(dest, nCount, GetData()->m_nDataLength - nCount, 0);
        return dest;
    }

    // upper/lower/reverse conversion
    void MakeUpper()
    {
        CopyBeforeWrite();
        appStrupr(m_pchData, GetLength() + 1, m_pchData);
    }
    void MakeLower()
    {
        CopyBeforeWrite();
        appStrlwr(m_pchData, GetLength() + 1, m_pchData);
    }
    void MakeReverse()
    {
        CopyBeforeWrite();
        appStrRev(m_pchData);
    }

    // trimming whitespace (either side)
    /**
    * find beginning of trailing matches
    * by starting at beginning (DBCS aware)
    */
    void TrimRight(const TCHAR* lpszTargetList)
    {
        CopyBeforeWrite();
        TCHAR* lpsz = m_pchData;
        TCHAR* lpszLast = NULL;

        while (*lpsz != _T('\0'))
        {
            if (appStrchr(lpszTargetList, *lpsz) != NULL)
            {
                if (lpszLast == NULL)
                    lpszLast = lpsz;
            }
            else
                lpszLast = NULL;
            lpsz = appStrInc(lpsz);
        }

        if (lpszLast != NULL)
        {
            // truncate at left-most matching character
            *lpszLast = _T('\0');
            GetData()->m_nDataLength = (INT)(lpszLast - m_pchData);
        }
    }
    /**
    * find beginning of trailing matches
    * by starting at beginning (DBCS aware)
    */
    void TrimRight(TCHAR chTarget)
    {
        CopyBeforeWrite();
        TCHAR* lpsz = m_pchData;
        TCHAR* lpszLast = NULL;

        while (*lpsz != _T('\0'))
        {
            if (*lpsz == chTarget)
            {
                if (lpszLast == NULL)
                    lpszLast = lpsz;
            }
            else
                lpszLast = NULL;
            lpsz = appStrInc(lpsz);
        }

        if (lpszLast != NULL)
        {
            // truncate at left-most matching character
            *lpszLast = _T('\0');
            GetData()->m_nDataLength = (INT)(lpszLast - m_pchData);
        }
    }
    /**
    * find beginning of trailing spaces by starting at beginning (DBCS aware)
    */
    void TrimRight()
    {
        CopyBeforeWrite();
        TCHAR* lpsz = m_pchData;
        TCHAR* lpszLast = NULL;

        while (*lpsz != _T('\0'))
        {
            if (appIsSpace(*lpsz))
            {
                if (lpszLast == NULL)
                    lpszLast = lpsz;
            }
            else
                lpszLast = NULL;
            lpsz = appStrInc(lpsz);
        }

        if (lpszLast != NULL)
        {
            // truncate at trailing space start
            *lpszLast = _T('\0');
            GetData()->m_nDataLength = (INT)(lpszLast - m_pchData);
        }
    }
    void TrimLeft(const TCHAR* lpszTargets)
    {
        // if we're not trimming anything, we're not doing any work
        if (__SafeStrlen(lpszTargets) == 0)
            return;

        CopyBeforeWrite();
        const TCHAR* lpsz = m_pchData;

        while (*lpsz != _T('\0'))
        {
            if (appStrchr(lpszTargets, *lpsz) == NULL)
                break;
            lpsz = appStrInc(lpsz);
        }

        if (lpsz != m_pchData)
        {
            // fix up data and length
            const INT nDataLength = GetData()->m_nDataLength - (INT)(lpsz - m_pchData);
            memmove(m_pchData, lpsz, (nDataLength + 1) * sizeof(TCHAR));
            GetData()->m_nDataLength = nDataLength;
        }
    }
    /**
    * find first non-matching character
    */
    void TrimLeft(TCHAR chTarget)
    {
        CopyBeforeWrite();
        const TCHAR* lpsz = m_pchData;

        while (chTarget == *lpsz)
            lpsz = appStrInc(lpsz);

        if (lpsz != m_pchData)
        {
            // fix up data and length
            const INT nDataLength = GetData()->m_nDataLength - (INT)(lpsz - m_pchData);
            memmove(m_pchData, lpsz, (nDataLength + 1) * sizeof(TCHAR));
            GetData()->m_nDataLength = nDataLength;
        }
    }
    void TrimLeft(); //out side

    // advanced manipulation
    //broken
    //INT Replace(TCHAR chOld, TCHAR chNew);
    //INT Replace(const TCHAR* lpszOld, const TCHAR* lpszNew);
    CCString Replace(TCHAR chOld, TCHAR chNew) const
    {
        TCHAR o[1];
        TCHAR n[1];
        o[0] = chOld;
        n[0] = chNew;
        return Replace(o, n);
    }
    CCString Replace(const TCHAR* lpszOld, const TCHAR* lpszNew) const;

    INT Remove(TCHAR chRemove);
    INT Insert(INT nIndex, TCHAR ch);
    INT Insert(INT nIndex, const TCHAR* pstr);
    INT Delete(INT nIndex, INT nCount = 1);

    // searching
    INT Find(const TCHAR* lpszSub, INT nStart) const
    {
        const INT nLength = GetData()->m_nDataLength;
        if (nStart > nLength)
            return -1;

        // find first matching substring
        TCHAR* lpsz = appStrstr(m_pchData + nStart, lpszSub);

        // return -1 for not found, distance from beginning otherwise
        return (NULL == lpsz) ? -1 : (INT)(lpsz - m_pchData);
    }
    INT Find(TCHAR ch, INT nStart) const
    {
        const INT nLength = GetData()->m_nDataLength;
        if (nStart >= nLength)
            return -1;

        // find first single character
        TCHAR* lpsz = appStrchr(m_pchData + nStart, ch);

        // return -1 if not found and index otherwise
        return (NULL == lpsz) ? -1 : (INT)(lpsz - m_pchData);
    }
    INT Find(TCHAR ch) const { return Find(ch, 0); }
    INT Find(const TCHAR* lpszSub) const { return Find(lpszSub, 0); }
    INT ReverseFind(TCHAR ch) const
    {
        TCHAR* lpsz = appStrrchr(m_pchData, ch);
        return (lpsz == NULL) ? -1 : (INT)(lpsz - m_pchData);
    }
    INT FindOneOf(const TCHAR* lpszCharSet) const
    {
        const TCHAR* lpsz = appStrpbrk(m_pchData, lpszCharSet);
        return (NULL == lpsz) ? -1 : (INT)(lpsz - m_pchData);
    }

    // simple formatting
    void __cdecl Format(const TCHAR* lpszFormat, ...);
    void FormatV(const TCHAR* lpszFormat, va_list argList);

    static CCString FormatS(const TCHAR* lpszFormat, ...);
    static CCString FormatVS(const TCHAR* lpszFormat, va_list argList);

    // Access to string implementation buffer as "C" character array
    TCHAR* GetBuffer(INT nMinBufLength);
    void ReleaseBuffer(INT nNewLength = -1);
    TCHAR* GetBufferSetLength(INT nNewLength);
    void FreeExtra();
    // Use LockBuffer/UnlockBuffer to turn ref counting off
    TCHAR* LockBuffer()
    {
        TCHAR* lpsz = GetBuffer(0);
        GetData()->m_nRefs = -1;
        return lpsz;
    }
    void UnlockBuffer() const; //out side

protected:

    TCHAR* m_pchData;   // pointer to ref counted string data

};

// Compare helpers
inline UBOOL operator==(const CCString& s1, const CCString& s2)
{
    return 0 == s1.Compare(s2);
}
inline UBOOL operator==(const CCString& s1, const TCHAR* s2)
{
    return 0 == s1.Compare(s2);
}
inline UBOOL operator==(const TCHAR* s1, const CCString& s2)
{
    return 0 == s2.Compare(s1);
}
inline UBOOL operator!=(const CCString& s1, const CCString& s2)
{
    return s1.Compare(s2) != 0;
}
inline UBOOL operator!=(const CCString& s1, const TCHAR* s2)
{
    return s1.Compare(s2) != 0;
}
inline UBOOL operator!=(const TCHAR* s1, const CCString& s2)
{
    return s2.Compare(s1) != 0;
}
inline UBOOL operator<(const CCString& s1, const CCString& s2)
{
    return s1.Compare(s2) < 0;
}
inline UBOOL operator<(const CCString& s1, const TCHAR* s2)
{
    return s1.Compare(s2) < 0;
}
inline UBOOL operator<(const TCHAR* s1, const CCString& s2)
{
    return s2.Compare(s1) > 0;
}
inline UBOOL operator>(const CCString& s1, const CCString& s2)
{
    return s1.Compare(s2) > 0;
}
inline UBOOL operator>(const CCString& s1, const TCHAR* s2)
{
    return s1.Compare(s2) > 0;
}
inline UBOOL operator>(const TCHAR* s1, const CCString& s2)
{
    return s2.Compare(s1) < 0;
}
inline UBOOL operator<=(const CCString& s1, const CCString& s2)
{
    return s1.Compare(s2) <= 0;
}
inline UBOOL operator<=(const CCString& s1, const TCHAR* s2)
{
    return s1.Compare(s2) <= 0;
}
inline UBOOL operator<=(const TCHAR* s1, const CCString& s2)
{
    return s2.Compare(s1) >= 0;
}
inline UBOOL operator>=(const CCString& s1, const CCString& s2)
{
    return s1.Compare(s2) >= 0;
}
inline UBOOL operator>=(const CCString& s1, const TCHAR* s2)
{
    return s1.Compare(s2) >= 0;
}
inline UBOOL operator>=(const TCHAR* s1, const CCString& s2)
{
    return s2.Compare(s1) <= 0;
}

// Globals
extern QLAPI const TCHAR* __GEmptyString;
#define __EmptyString ( (__NAMESPACE::CCString&) *(__NAMESPACE::CCString*)&(__NAMESPACE::__GEmptyString) ) 
inline const CCString& appGetEmptyString() { return __EmptyString; }

inline CCString::CCString() { Init(); }

inline void CCString::Init()
{
    m_pchData = __EmptyString.m_pchData;
}

template<class T>
inline CCString appToString(const T& content)
{
    CCString res;
    res.Format(_T("%d"), content);
    return res;
}

defineToStr(INT)
defineToStr(UINT)
defineToStr(WORD)
defineToStr(SWORD)
defineToStr(BYTE)
defineToStr(SBYTE)
defineToStr(LONGLONG)
defineToStr(ULONGLONG)

defineToStr(Real)
#if _CLG_DOUBLEFLOAT
defineToStr(FLOAT)
#else
defineToStr(DOUBLE)
#endif

template<>
inline CCString appToString(const CCString& content)
{
    return content;
}

template<>
inline CCString appToString(const QLComplex& content)
{
    CCString sret;
    sret.Format(_T("%2.12f %s %2.12f I"),
        content.x,
        content.y < F(0.0) ? _T("-") : _T("+"),
        abs(content.y));
    return sret;
}

#if _CLG_DOUBLEFLOAT
template<>
inline CCString appToString(const cuComplex& content)
{
    CCString sret;
    sret.Format(_T("%2.12f %s %2.12f I"),
        content.x,
        content.y < F(0.0) ? _T("-") : _T("+"),
        abs(content.y));
    return sret;
}
#else
template<>
inline CCString appToString(const cuDoubleComplex& content)
{
    CCString sret;
    sret.Format(_T("%2.12f %s %2.12f I"),
        content.x,
        content.y < F(0.0) ? _T("-") : _T("+"),
        abs(content.y));
    return sret;
}
#endif

template <class T>
inline CCString appToString(const TArray<T>& arr)
{
    CCString sret = _T("{");
    for (INT i = 0; i < arr.Num(); ++i)
    {
        if (i == arr.Num() - 1)
        {
            sret = sret + appToString(arr[i]);
        }
        else
        {
            sret = sret + appToString(arr[i]) + _T(", ");
        }
    }
    sret = sret + _T("}");
    return sret;
}

/**
* The thing is, since we don't have custum types (like vector, rotation or so) to sprint
* Just use the standard one
*/
inline CCString appStringFormatV(const TCHAR* lpszFormat, va_list argList)
{
    CCString sRet;
    sRet.FormatV(lpszFormat, argList);
    return sRet;
}

/**
*
*/
inline CCString appStringFormat(const TCHAR* lpszFormat, ...)
{
    CCString outString;
    va_list argList;
    va_start(argList, lpszFormat);
    outString = appStringFormatV(lpszFormat, argList);
    va_end(argList);
    return outString;
}

/**
*
*/
enum EGetStringListFlag
{
    EGSLF_CutHead = 0x00000001,
    EGSLF_CutTail = 0x00000002,
    EGSLF_IgnorEmety = 0x00000004,
    EGSLF_IgnorTabSpace = 0x0000008,
    EGSLF_IgnorTabSpaceInSide = 0x00000010,

    EGSLF_ForceDWORD = 0x7fffffff,
};

/**
* for example, when use appGetStringList( _T("asdad.asdasda..gagaga", '.', EGSLF_IgnorEmety);
* there is "asdad"  "asdasda"  "gagaga" if not EGSLF_IgnorEmety, it's  "asdad"  "asdasda" ""  "gagaga"
* appGetStringList( _T(".asdad.asdasda.gagaga", '.', EGSLF_CutHead);
* there is(or use EGSLF_IgnorEmety) "asdad"  "asdasda"  "gagaga", if not , it's "" "asdad"  "asdasda"  "gagaga"
* appGetStringList( _T(".asdad.asdasda.   .gagaga", '.', EGSLF_IgnorTabSpace);
* there is "asdad"  "asdasda"  "gagaga" if not EGSLF_IgnorTabSpace, it's  "asdad"  "asdasda" "  "  "gagaga"
* appGetStringList( _T(".asdad.asdasda.   .    gag   aga ", '.', EGSLF_IgnorTabSpaceInSide);
* there is "asdad"  "asdasda" ""  "gagaga" if not EGSLF_IgnorTabSpaceInSide, it's  "asdad"  "asdasda" "  "  "    gaga  ga"
*/
inline TArray<CCString> appGetStringList(const CCString& orignString, TArray<INT> seperate, UINT dwFlag = 0)
{
    TArray<CCString> outList;
    const CCBitFlag flag(dwFlag);
    const TCHAR* lpsz = (const TCHAR*)(orignString);
    TCHAR lp[CCString::MAX_PATH];
    memset(lp, 0, CCString::MAX_PATH * sizeof(TCHAR));
    INT currentIndex = 0;

    while (1)
    {
        UBOOL bIsSeperate = FALSE;
        for (INT i = 0; i < seperate.Num(); ++i)
            bIsSeperate = bIsSeperate || (*lpsz == seperate[i]);

        bIsSeperate = bIsSeperate || (*lpsz == 0) || (*lpsz == '\n') || (*lpsz == '\r');

        const UBOOL bTail = (*lpsz == 0);

        if (bIsSeperate)
        {
            UBOOL bRecord = TRUE;
            if (0 == currentIndex)
            {
                if (flag.HasFlag(EGSLF_IgnorEmety) || flag.HasFlag(EGSLF_IgnorTabSpace))
                    bRecord = FALSE;

                if (flag.HasFlag(EGSLF_CutHead) && 0 == outList.Num())
                    bRecord = FALSE;

                if (flag.HasFlag(EGSLF_CutTail) && bTail)
                    bRecord = FALSE;
            }
            else
            {
                if (flag.HasFlag(EGSLF_IgnorTabSpace))
                {
                    UBOOL bOnlyTabSpace = TRUE;
                    for (INT i = 0; i < currentIndex; ++i)
                        bOnlyTabSpace = bOnlyTabSpace && (lp[i] == ' ' || lp[i] == '\t' || lp[i] == 0);

                    bRecord = !bOnlyTabSpace;
                }
            }

            if (bRecord)
            {
                CCString sNewWord(lp);
                outList.AddItem(sNewWord);
            }

            memset(lp, 0, CCString::MAX_PATH * sizeof(TCHAR));
            currentIndex = 0;
            ++lpsz;
        }
        else
        {
            if (!flag.HasFlag(EGSLF_IgnorTabSpaceInSide) || (!(' ' == (*lpsz) || '\t' == (*lpsz))))
            {
                lp[currentIndex] = *lpsz;
                ++currentIndex;
            }
            ++lpsz;
        }

        if (bTail)
            break;
    }

    return outList;
}

inline TArray<CCString> appGetStringList(const CCString& orignString, INT seperate, UINT dwFlag = 0)
{
    TArray <INT> inSep;
    inSep.AddItem(seperate);
    return appGetStringList(orignString, inSep, dwFlag);
}

/**
* for example, when use appGetCutTail( _T("asdasda.bmp", '.', TRUE);
* there is "asdasda"
* when use appGetCutTail( _T("data\asdasda.bmp", '\\');
* there is "data\"
* when cut head, it's "\asdasda.bmp"
*/
inline CCString appGetCutTail(const CCString& orignString, TArray<INT> seperate, UBOOL bCutSep = FALSE, UBOOL bCutHead = FALSE)
{
    const TCHAR* lpsz = (const TCHAR*)(orignString);
    TCHAR lp[CCString::MAX_PATH];
    memset(lp, 0, CCString::MAX_PATH * sizeof(TCHAR));
    appStrcpy(lp, CCString::MAX_PATH, lpsz);
    INT location = -1;
    for (INT i = 0; i < seperate.Num(); ++i)
    {
        const INT location1 = bCutHead ? orignString.Find(static_cast<TCHAR>(seperate[i])) : orignString.ReverseFind(static_cast<TCHAR>(seperate[i]));

        //no need to do this
//         if (-1 == location1)
//             continue;

        if (-1 == location)
            location = location1;
        else
            location = bCutHead ? std::min(location, location1) : std::max(location, location1);
    }

    if (location < 0)
        return orignString;

    if (bCutHead)
    {
        if (bCutSep)
            ++location;
        memset(lp, 0, CCString::MAX_PATH * sizeof(TCHAR));
        appStrcpy(lp, CCString::MAX_PATH, (lpsz + location));
    }
    else
    {
        if (bCutSep)
            *(lp + location) = 0;
        else
            *(lp + location + 1) = 0;
    }

    CCString outString(lp);
    return outString;
}

inline CCString appGetCutTail(const CCString& orignString, INT seperate, UBOOL bCutSep = FALSE, UBOOL bCutHead = FALSE)
{
    TArray <INT> inSep;
    inSep.AddItem(seperate);
    return appGetCutTail(orignString, inSep, bCutSep, bCutHead);
}

inline CCString appGetHex(FLOAT f)
{
    UINT* pV = (UINT*) &f;
    CCString sRet;
    for (UINT i = 0; i < 8; ++i)
    {
        const BYTE byCurrent = static_cast<BYTE>( ((*pV) & (0xf << (i * 4))) >>  (i * 4) );
        const TCHAR sThisChar = byCurrent > 9 ? byCurrent - 10 + _T('A') : byCurrent + _T('0');
        sRet += (TCHAR)sThisChar;
    }
    printf("%f %s\n", f, sRet.c_str());

    return sRet;
}

__END_NAMESPACE

#endif //#ifndef _CCSTRING_H_

//=============================================================================
// END OF FILE
//=============================================================================

