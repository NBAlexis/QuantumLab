//=============================================================================
// FILENAME : QLSimulator.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

std::string QLSimulator::PrintComplex(Real fReal, Real fImg) const
{
    static TCHAR bufferf[256];
    static TCHAR buffer[256];
    if (absReal(fImg) > REAL_EPS)
    {
        if (fImg > 0)
        {
            sprintf_(bufferf, 256, "%s + %s I", m_sFloatFormat.c_str(), m_sFloatFormat.c_str());
        }
        else
        {
            sprintf_(bufferf, 256, "%s - %s I", m_sFloatFormat.c_str(), m_sFloatFormat.c_str());
        }
        sprintf_(buffer, 256, bufferf, fReal, abs(fImg));
    }
    else
    {
        sprintf_(bufferf, 256, "%s", m_sFloatFormat.c_str());
        sprintf_(buffer, 256, bufferf, fReal);
    }
    return std::string(buffer);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================