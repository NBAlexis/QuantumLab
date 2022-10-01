//=============================================================================
// FILENAME : AmplitudeEncode.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#ifndef _AMPLITUDEENCODE_H_
#define _AMPLITUDEENCODE_H_

__BEGIN_NAMESPACE

static inline UINT MostSignificantPowerTwo(UINT n)
{
    UINT k = 0;
    while ((1U << k) < n)
    {
        ++k;
    }
    return k;
}

extern TArray<QLComplex> QLAPI NormalizeV(const TArray<QLComplex>& v, UINT& lenPower);

extern QLGate QLAPI AmplitudeEncode(const TArray<QLComplex>& v);

__END_NAMESPACE


#endif //#ifndef _AMPLITUDEENCODE_H_

//=============================================================================
// END OF FILE
//=============================================================================