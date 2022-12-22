//=============================================================================
// FILENAME : FRGate.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [30/09/2022 nbale]
//=============================================================================

#ifndef _FRGATE_H_
#define _FRGATE_H_

__BEGIN_NAMESPACE

inline UINT GrayCode(UINT n)
{
    return n ^ (n >> 1);
}

/**
* see http://graphics.stanford.edu/~seander/bithacks.html
* n is expected to be a 32-bit unsigned number and is power-2
*/
inline UINT Log2(INT n)
{
    static const int MultiplyDeBruijnBitPosition[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };
    return MultiplyDeBruijnBitPosition[((UINT)((n & -n) * 0x077CB531U)) >> 27];
}

/**
* see http://graphics.stanford.edu/~seander/bithacks.html
*/
inline UINT Counting1Bit(UINT n)
{
    UINT c;
    for (c = 0; n; ++c)
    {
        n &= n - 1;
    }
    return c;
}

inline UINT GrayCodeDifferent(UINT n, UINT maxN)
{
    UINT a = GrayCode(n);
    UINT b = (n == maxN - 1) ? GrayCode(0) : GrayCode(n + 1);
    UINT c = a ^ b;
    return Log2(c);
}

inline UINT BitWiseInnerProduct(UINT a, UINT b)
{
    return Counting1Bit(a & b);
}

extern TArray<Real> QLAPI SpliteAngles(const TArray<Real>& angles, UINT length);
extern TArray<Real> QLAPI SpliteAngles(const Real* angles, UINT length);

/**
* Uniform controlled rotation in 10.1103/PhysRevLett.93.130502
*/
extern QLGate QLAPI FRy(const TArray<Real>& angles, UINT numberOfQubits);

extern QLGate QLAPI FRz(const TArray<Real>& angles, UINT numberOfQubits);

/**
* Used in Cosin-Sine decomposition, when FRz-dagger is appended after FRy, two CNots can be removed
*/
extern QLGate QLAPI FRyz(const TArray<Real>& anglesY, const TArray<Real>& anglesZ, UINT numberOfQubits);

extern QLGate QLAPI FRyz(const Real* anglesY, const Real* anglesZ, UINT numberOfQubits);

extern QLGate QLAPI FRp(const TArray<Real>& angles, UINT numberOfQubits);

extern QLGate QLAPI FRPh(const TArray<Real>& angles, UINT numberOfQubits);

__END_NAMESPACE


#endif //#ifndef _FRGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================