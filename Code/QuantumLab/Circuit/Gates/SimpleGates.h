//=============================================================================
// FILENAME : SimpleGate.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [12/09/2022 nbale]
//=============================================================================

#ifndef _SIMPLEGATE_H_
#define _SIMPLEGATE_H_

__BEGIN_NAMESPACE

static inline CCString Binary(UINT num, INT max)
{
    CCString sret;
    for (INT i = max - 1; i >= 0; --i)
    {
        if (num & (1U << i))
        {
            sret = sret + _T("1");
        }
        else
        {
            sret = sret + _T("0");
        }
    }
    return sret;
}

extern QLGate QLAPI CreateZYZGate(const QLMatrix& u, UBOOL bNormalize = TRUE);

extern QLGate QLAPI CreateControlledZYZGate(const QLMatrix& u, UBOOL bNormalize = TRUE);

extern QLGate QLAPI CreateSwapGate();

extern QLGate QLAPI CreateControlledHadamardGate();

__END_NAMESPACE


#endif //#ifndef _SIMPLEGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================