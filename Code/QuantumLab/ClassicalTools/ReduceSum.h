//=============================================================================
// FILENAME : ReduceSum.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/12/2022 nbale]
//=============================================================================

#ifndef _REDUCESUM_H_
#define _REDUCESUM_H_

__BEGIN_NAMESPACE

static inline UINT GetReduceDim(UINT uiLength)
{
    UINT iRet = 0;
    while ((1U << iRet) < uiLength)
    {
        ++iRet;
    }
    return iRet;
}

/**
* reduce sum
* 'value' will be changed, the first element is the result
*/
extern QLAPI Real ReduceSum(Real* value, UINT count);

extern QLAPI INT ReduceSum(INT* value, UINT count);

extern QLAPI UINT ReduceSum(UINT* value, UINT count);

extern QLAPI FLOAT ReduceSum(FLOAT* value, UINT count);

/**
* reduce sum
* 'value' will be changed, the first element is the result
*/
extern QLAPI QLComplex ReduceSum(QLComplex* value, UINT count);

/**
* sum _{if cndition[n] = conditionEqual} value[n * stride + offset]
* work space must has length >= (len(v) + 1) / 2
*/
extern QLAPI FLOAT ConditionalSum(const FLOAT* value, BYTE byStride, BYTE offset, const BYTE* condition, BYTE conditionEqual, UINT count, FLOAT* workSpace);

/**
* sum _{if cndition[n] = conditionEqual} 1
* work space must has length >= (len(condition) + 1) / 2
*/
extern QLAPI UINT ConditionalCount(const INT* condition, INT conditionEqual, UINT count, UINT* workSpace);

extern QLAPI UINT ConditionalCount(const BYTE* condition, BYTE conditionEqual, UINT count, UINT* workSpace);

/**
* find index of the max
*/
extern QLAPI UINT Max(const INT* v, UINT* workspace, UINT count);

__END_NAMESPACE


#endif //#ifndef _REDUCESUM_H_

//=============================================================================
// END OF FILE
//=============================================================================