//=============================================================================
// FILENAME : Tracer.h
// 
// DESCRIPTION:
// This is class for messages
//
// REVISION: [dd/mm/yy]
//  [13/09/2022 nbale]
//=============================================================================

#ifndef _TRACER_H_
#define _TRACER_H_

__BEGIN_NAMESPACE

enum class EVerboseLevel : UINT
{
    CRUCIAL,
    WARNING,
    GENERAL,
    DETAILED,
    PARANOIAC,
};

enum 
{
    _kTraceBuffSize = 4096,
};

class QLAPI CTracer
{
public:
    CTracer(void)
        : m_eLevel(EVerboseLevel::PARANOIAC)
        , m_pStream(NULL)
        , m_pStdStream(NULL)
        , m_bLogDate(TRUE)
    {
        Initial(EVerboseLevel::PARANOIAC);
    }

    ~CTracer(void)
    {
        if (NULL != m_pStream)
        {
            m_pStream->flush();
        }
        appSafeDelete(m_pStream);
    }

    inline void SetVerboseLevel(EVerboseLevel eLevel) { m_eLevel = eLevel; }

    inline void SetOutStream(const CCString& filename = _T("datetime"))
    {
        appSafeDelete(m_pStdStream);
        if (NULL != m_pStream)
        {
            m_pStream->flush();
            appSafeDelete(m_pStream);
        }

        m_pStdStream = new std::ostream(std::cout.rdbuf());
        UBOOL bShowHasFile = FALSE;
        if (filename == _T("stdout"))
        {
            m_pStream = NULL;
        }
        else if (filename == _T("timestamp"))
        {
            const CCString sRealFile = appStringFormat(_T("%d.log"), appGetTimeStamp());
            m_pStream = new std::ofstream(sRealFile);
            bShowHasFile = TRUE;
        }
        else if (filename == _T("datetime"))
        {
            static TCHAR datetime[256];
            appGetTimeNow(datetime, 256);
            const CCString sRealFile = appStringFormat(_T("%s.log"), datetime);
            m_pStream = new std::ofstream(sRealFile);
            bShowHasFile = TRUE;
        }
        else
        {
            m_pStream = new std::ofstream(filename);
            bShowHasFile = TRUE;
        }

        if (NULL == m_pStdStream || (bShowHasFile && NULL == m_pStream))
        {
            printf(_T("ERROR: CTracer: no output stream."));
            if (NULL != m_pStream)
            {
                m_pStream->flush();
            }
            exit(EXIT_FAILURE);
        }
    }

    inline void Initial(EVerboseLevel eLevel = EVerboseLevel::PARANOIAC, const CCString& filename = _T("datetime"), const CCString& sFloatFormat = _T("%.4f"))
    {
        m_eLevel = eLevel;
        m_sFloatFormat = sFloatFormat;
        m_pStdStream = new std::ostream(std::cout.rdbuf());
        UBOOL bShowHasFile = FALSE;
        if (filename == _T("stdout"))
        {
            m_pStream = NULL;
        }
        else if (filename == _T("timestamp"))
        {
            const CCString sRealFile = appStringFormat(_T("%d.log"), appGetTimeStamp());
            m_pStream = new std::ofstream(sRealFile);
            bShowHasFile = TRUE;
        }
        else if (filename == _T("datetime"))
        {
            static TCHAR datetime[256];
            appGetTimeNow(datetime, 256);
            const CCString sRealFile = appStringFormat(_T("%s.log"), datetime);
            m_pStream = new std::ofstream(sRealFile);
            bShowHasFile = TRUE;
            printf("initialled file: %s \n", sRealFile.c_str());
        }
        else
        {
            m_pStream = new std::ofstream(filename);
            bShowHasFile = TRUE;
        }

        if (NULL == m_pStdStream || (bShowHasFile && NULL == m_pStream))
        {
            printf(_T("ERROR: CTracer: no output stream."));
            if (NULL != m_pStream)
            {
                m_pStream->flush();
            }
            exit(EXIT_FAILURE);
        }
    }

    inline void Print(EVerboseLevel level, const TCHAR *format, va_list& arg)
    {
        if ((level <= m_eLevel))
        {
            //assert(NULL != m_pStdStream);
            if (NULL == m_pStdStream)
            {
                //Maybe the first initial is not entered?
            }

            if (EVerboseLevel::CRUCIAL == level && NULL != m_pStdStream)
            {
                //red bold
                *m_pStdStream << _T("\033[31;1m");
            }
            else if (EVerboseLevel::WARNING == level && NULL != m_pStdStream)
            {
                *m_pStdStream << _T("\033[35m");
            }
            else if (EVerboseLevel::DETAILED == level && NULL != m_pStdStream)
            {
                //green
                *m_pStdStream << _T("\033[32m");
            }
            else if (EVerboseLevel::PARANOIAC == level && NULL != m_pStdStream)
            {
                //dark yellow
                *m_pStdStream << _T("\033[33m");
            }

            if (m_bLogDate)
            {
                static TCHAR timeBuffer[256];
                if (level <= EVerboseLevel::GENERAL)
                {
                    appGetTimeNow(timeBuffer, 256);
                    if (NULL != m_pStdStream)
                    {
                        *m_pStdStream << _T("[") << timeBuffer << "|" << m_sTraceHeader.c_str() << _T("]");
                    }
                    if (NULL != m_pStream)
                    {
                        *m_pStream << _T("[") << timeBuffer << "|" << m_sTraceHeader.c_str() << _T("]");
                    }
                }
            }

            appVsnprintf(m_cBuff, _kTraceBuffSize - 1, format, arg);
            if (NULL != m_pStdStream)
            {
                *m_pStdStream << m_cBuff;
            }
            
            if (NULL != m_pStream)
            {
                *m_pStream << m_cBuff;
#ifdef _QL_DEBUG
                *m_pStream << std::flush;
#endif
            }

            if ((EVerboseLevel::CRUCIAL == level || EVerboseLevel::WARNING == level || EVerboseLevel::PARANOIAC == level || EVerboseLevel::DETAILED == level) && NULL != m_pStdStream)
            {
                *m_pStdStream << _T("\033[0m");
            }
        }
    }

    inline void Flush() const
    {
        if (NULL != m_pStream)
        {
            m_pStream->flush();
        }
    }

    inline void SetLogDate(UBOOL bLog)
    {
        m_bLogDate = bLog;
        m_bLogDateHist.RemoveAll();
    }
    inline void PushLogDate(UBOOL bLog)
    {
        m_bLogDateHist.AddItem(m_bLogDate);
        m_bLogDate = bLog;
    }
    inline void PopLogDate()
    {
        if (m_bLogDateHist.Num() > 0)
        {
            m_bLogDate = m_bLogDateHist.Pop();
        }
    }
    inline void SetLogLevel(EVerboseLevel eLevel)
    {
        m_eLevel = eLevel;
        m_eLogLevelHist.RemoveAll();
    }
    inline void PushLogLevel(EVerboseLevel eLevel)
    {
        m_eLogLevelHist.AddItem(m_eLevel);
        m_eLevel = eLevel;
    }
    inline void PopLogLevel()
    {
        if (m_eLogLevelHist.Num() > 0)
        {
            m_eLevel = m_eLogLevelHist.Pop();
        }
    }

    inline void SetLogHeader(const CCString& sHeader) { m_sTraceHeader = sHeader; }

    inline void SetFloatFormat(const CCString& sFormat) { m_sFloatFormat = sFormat; }

    static CCString PrintComplex(Real fReal, Real fImg, const CCString& sFloatFormat, UINT length = 0);

    CCString PrintComplex(Real fReal, Real fImg, UINT length = 0) const
    {
        return PrintComplex(fReal, fImg, m_sFloatFormat, length);
    }

    CCString GetFloatFormat() const { return m_sFloatFormat; }

private:

    EVerboseLevel m_eLevel;
    std::ostream* m_pStream;
    std::ostream* m_pStdStream;
    TCHAR m_cBuff[_kTraceBuffSize];
    UBOOL m_bLogDate;
    TArray<UBOOL> m_bLogDateHist;
    TArray<EVerboseLevel> m_eLogLevelHist;
    CCString m_sTraceHeader;
    CCString m_sFloatFormat;
};

extern QLAPI void appInitialTracer(EVerboseLevel eLevel = EVerboseLevel::PARANOIAC, const CCString& filename = _T("datetime"), const CCString& sFloatFormat = _T("%.6f"));
extern QLAPI void appVOut(EVerboseLevel eLevel, const TCHAR *format, ...);
extern QLAPI void _appCrucial(const TCHAR *format, ...);
extern QLAPI void _appWarning(const TCHAR* format, ...);
extern QLAPI void appGeneral(const TCHAR *format, ...);
extern QLAPI void appDetailed(const TCHAR *format, ...);
extern QLAPI void appParanoiac(const TCHAR *format, ...);

#define appAssert(exp) { if (!(exp)) { appCrucial(_T("assert failed %s\n"), _T(#exp)); } }
#define appCrucial(...) {TCHAR ___msg[1024];appSprintf(___msg, 1024, __VA_ARGS__);_appCrucial(_T("%s(%d): Error: %s\n"), _T(__FILE__), __LINE__, ___msg);}
#define appWarning(...) {TCHAR ___msg[1024];appSprintf(___msg, 1024, __VA_ARGS__);_appWarning(_T("%s(%d): Warning: %s\n"), _T(__FILE__), __LINE__, ___msg);}

extern QLAPI CTracer GTracer;

inline void appSetTracer(EVerboseLevel eLevel, const CCString& filename, const CCString& sFloatFormat = _T("%.6f"))
{
    GTracer.SetVerboseLevel(eLevel);
    GTracer.SetOutStream(filename);
    GTracer.SetFloatFormat(sFloatFormat);
}

inline void appSetFloatFormat(const CCString& sFloatFormat)
{
    GTracer.SetFloatFormat(sFloatFormat);
}

inline void appFlushLog()
{
    GTracer.Flush();
}

inline void appSetLogDate(UBOOL bLog)
{
    GTracer.SetLogDate(bLog);
}

inline void appPushLogDate(UBOOL bLog)
{
    GTracer.PushLogDate(bLog);
}

inline void appPopLogDate()
{
    GTracer.PopLogDate();
}

inline void appSetLogLevel(EVerboseLevel eLevel)
{
    GTracer.SetLogLevel(eLevel);
}

inline void appPushLogLevel(EVerboseLevel eLevel)
{
    GTracer.PushLogLevel(eLevel);
}

inline void appPopLogLevel()
{
    GTracer.PopLogLevel();
}

inline void appSetLogHeader(const CCString& sHeader)
{
    GTracer.SetLogHeader(sHeader);
}

inline CCString appPrintComplex(Real fReal, Real fImg, UINT length = 0)
{
    return GTracer.PrintComplex(fReal, fImg, length);
}

__END_NAMESPACE

#endif //_TRACER_H_

//=============================================================================
// END OF FILE
//=============================================================================
