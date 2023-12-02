//=============================================================================
// FILENAME : Tracer.cpp
// 
// DESCRIPTION:
// This is class for messages
//
// REVISION: [dd/mm/yy]
//  [13/09/2022 nbale]
//=============================================================================
#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLAPI CTracer GTracer;

/**
*
*
*/
QLAPI void appInitialTracer(EVerboseLevel eLevel, const CCString& filename, const CCString& sFloatFormat)
{
    GTracer.Initial(eLevel, filename, sFloatFormat);
}

/**
*
*
*/
QLAPI void appVOut(EVerboseLevel level, const TCHAR *format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(level, format, arg);
        va_end(arg);
    }
}

/**
*
*
*/
QLAPI void _appCrucial(const TCHAR *format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(EVerboseLevel::CRUCIAL, format, arg);
        GTracer.Flush();
        va_end(arg);
    }

#if _QL_DEBUG
#if _QL_WIN
    __debugbreak();
#elif defined(SIGTRAP)
    raise(SIGTRAP);
#else
    __builtin_trap();
#endif
#endif

}

QLAPI void _appWarning(const TCHAR* format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(EVerboseLevel::WARNING, format, arg);
        GTracer.Flush();
        va_end(arg);
    }
}

/**
*
*
*/
QLAPI void appGeneral(const TCHAR *format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(EVerboseLevel::GENERAL, format, arg);
        va_end(arg);
    }
}

/**
*
*
*/
QLAPI void appDetailed(const TCHAR *format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(EVerboseLevel::DETAILED, format, arg);
        va_end(arg);
    }
}

/**
*
*
*/
QLAPI void appParanoiac(const TCHAR *format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(EVerboseLevel::PARANOIAC, format, arg);
        va_end(arg);
    }
}

CCString CTracer::PrintComplex(Real fReal, Real fImg, const CCString& sFloatFormat, UINT length)
{
    static TCHAR bufferf[256];
    static TCHAR buffer[256];
    if (absReal(fImg) > REAL_EPS)
    {
        if (fImg > 0)
        {
            sprintf_(bufferf, 256, "%s + %s I", sFloatFormat.c_str(), sFloatFormat.c_str());
        }
        else
        {
            sprintf_(bufferf, 256, "%s - %s I", sFloatFormat.c_str(), sFloatFormat.c_str());
        }
        sprintf_(buffer, 256, bufferf, fReal, abs(fImg));
    }
    else
    {
        sprintf_(bufferf, 256, "%s", sFloatFormat.c_str());
        sprintf_(buffer, 256, bufferf, fReal);
    }
    if (length > 0)
    {
        INT bufferL = static_cast<INT>(appStrlen(buffer));
        CCString sSpaces;
        for (INT i = 0; i < static_cast<INT>(length) - bufferL; ++i)
        {
            sSpaces = sSpaces + _T(" ");
        }
        return sSpaces + CCString(buffer);
    }

    return CCString(buffer);
}

__END_NAMESPACE

//====================================================================
//====================================================================
