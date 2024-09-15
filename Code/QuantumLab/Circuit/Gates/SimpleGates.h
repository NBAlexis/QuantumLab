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

/**
* U = [u11 u12]
*     [u21 u22]
*
* U = exp(i delta) [ exp(i (alpha+beta)/2) cos(theta/2)  exp(i (alpha-beta)/2) sin(theta/2)] = exp(i delta) Rz(-alpha) Ry(-theta) Rz(-beta)
*                  [-exp(-i(alpha-beta)/2) sin(theta/2)  exp(-i(alpha+beta)/2) cos(theta/2)]
*
* for OpenQASM
* 
* 
* 
* U(a,b,c)=Rz(b)Ry(a)Rz(c)|psi>, where a=-theta, b=-alpha, c=-beta
* a,b,c = degree(0), degree(2), degree(1)
*/
extern TArray<Real> QLAPI GetZYZDecompose(const QLMatrix& u);

extern QLGate QLAPI CreateZYZGate(const QLMatrix& u);

extern QLGate QLAPI CreateControlledZYZGate(const QLMatrix& u);

extern QLGate QLAPI CreateSwapGate();

extern QLGate QLAPI CreateControlledSwap(BYTE numOfController);

extern QLGate QLAPI CreateControlledHadamardGate();

__END_NAMESPACE


#endif //#ifndef _SIMPLEGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================