//=============================================================================
// FILENAME : QuantumDataViewer.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [08/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

#pragma region kernel

template<class T>
__global__ void
_QL_LAUNCH_BOUND
_kernelReduceReal(T* arr, UINT uiJump, UINT uiMax)
{
    //for length 16 array
    //for jump = 1, this is 1->0, 3->2, 5->4, 7->6, 9->10, 11->10, 13->12, 15->14 
    //for jump = 2, this is 2->0, 6->4, 10->8, 14->12 
    //for jump = 4, this is 4->0, 12->8 
    //for jump = 8, this is 8->0, and is finished.

    //id target = idx * (jump << 1)
    //id from = target + jump
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] += arr[uiIdFrom];
    }
}

template<class T>
__global__ void
_QL_LAUNCH_BOUND
_kernelReduceMax(const T * __restrict__ arr, UINT* idxbuffer, UINT uiJump, UINT uiMax)
{
    const UINT uiIdTo = (threadIdx.x + blockIdx.x * blockDim.x)* (uiJump << 1);
    const UINT uiIdFrom = uiIdTo + uiJump;
    if (uiIdFrom < uiMax)
    {
        if (1 == uiJump)
        {
            idxbuffer[uiIdFrom] = uiIdFrom;
            idxbuffer[uiIdTo] = uiIdTo;
        }

        if (arr[idxbuffer[uiIdFrom]] > arr[idxbuffer[uiIdTo]])
        {
            idxbuffer[uiIdTo] = idxbuffer[uiIdFrom];
        }
    }
}

template<class T>
__global__ void
_QL_LAUNCH_BOUND
_kernelReduceMin(const T* __restrict__ arr, UINT* idxbuffer, UINT uiJump, UINT uiMax)
{
    const UINT uiIdTo = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1);
    const UINT uiIdFrom = uiIdTo + uiJump;
    if (uiIdFrom < uiMax)
    {
        if (1 == uiJump)
        {
            idxbuffer[uiIdFrom] = uiIdFrom;
            idxbuffer[uiIdTo] = uiIdTo;
        }

        if (arr[idxbuffer[uiIdFrom]] < arr[idxbuffer[uiIdTo]])
        {
            idxbuffer[uiIdTo] = idxbuffer[uiIdFrom];
        }
    }
}

template<class T>
__global__ void
_QL_LAUNCH_BOUND
_kernelReduceMax(T* arr, UINT uiJump, UINT uiMax)
{
    const UINT uiIdTo = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1);
    const UINT uiIdFrom = uiIdTo + uiJump;
    if (uiIdFrom < uiMax)
    {
        if (arr[uiIdTo] < arr[uiIdFrom])
        {
            arr[uiIdTo] = arr[uiIdFrom];
        }
    }
}

template<class T>
__global__ void
_QL_LAUNCH_BOUND
_kernelReduceMin(T* arr, UINT uiJump, UINT uiMax)
{
    const UINT uiIdTo = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1);
    const UINT uiIdFrom = uiIdTo + uiJump;
    if (uiIdFrom < uiMax)
    {
        if (arr[uiIdTo] > arr[uiIdFrom])
        {
            arr[uiIdTo] = arr[uiIdFrom];
        }
    }
}

__global__ void
_QL_LAUNCH_BOUND
_kernelReduceComp(QLComplex* arr, UINT uiJump, UINT uiMax)
{
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] = _cuCaddf(arr[uiIdFrom - uiJump], arr[uiIdFrom]);
    }
}

__global__ void
_QL_LAUNCH_BOUND
_kernelReduceCompDB(cuDoubleComplex* arr, UINT uiJump, UINT uiMax)
{
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] = cuCadd(arr[uiIdFrom - uiJump], arr[uiIdFrom]);
    }
}

template<class T, class Tc>
__global__ void
_QL_LAUNCH_BOUND
_kernelConditionalReduceReal(const T * __restrict__ arr, T* res, BYTE stride, BYTE offset, UINT uiMax, const Tc * __restrict__ condition, Tc condEqual)
{
    const UINT uiId = (threadIdx.x + blockIdx.x * blockDim.x);
    UINT uiId2 = (uiId << 1);
    res[uiId] = (condition[uiId2] == condEqual) ? arr[uiId2 * stride + offset] : static_cast<T>(0);
    ++uiId2;
    if (uiId2 < uiMax && condition[uiId2] == condEqual)
    {
        res[uiId] += arr[uiId2 * stride + offset];
    }
}

template<class T, class Tc>
__global__ void
_QL_LAUNCH_BOUND
_kernelConditionalCount(T* res, UINT uiMax, const Tc * __restrict__ condition, Tc condEqual)
{
    const UINT uiId = (threadIdx.x + blockIdx.x * blockDim.x);
    UINT uiId2 = (uiId << 1);
    res[uiId] = (condition[uiId2] == condEqual) ? static_cast<T>(1) : static_cast<T>(0);
    ++uiId2;
    if (uiId2 < uiMax && condition[uiId2] == condEqual)
    {
        res[uiId] += static_cast<T>(1);
    }
}

#pragma endregion

template<class T>
T ReduceSumT(T* value, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceReal << <iBlock, iThread >> > (value, iJump, count);
    }
    T result[1];
    cudaMemcpy(result, value, sizeof(T), cudaMemcpyDeviceToHost);
    return result[0];
}

template<class T>
void ReduceSumT(T* res, T* value, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceReal << <iBlock, iThread >> > (value, iJump, count);
    }
    cudaMemcpy(res, value, sizeof(T), cudaMemcpyDeviceToDevice);
}

template<class T>
UINT ReduceMaxT(const T* value, UINT* idxBuffer, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceMax << <iBlock, iThread >> > (value, idxBuffer, iJump, count);
    }
    UINT result[1];
    cudaMemcpy(result, idxBuffer, sizeof(UINT), cudaMemcpyDeviceToHost);
    return result[0];
}

template<class T>
UINT ReduceMinT(const T* value, UINT* idxBuffer, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceMin << <iBlock, iThread >> > (value, idxBuffer, iJump, count);
    }
    UINT result[1];
    cudaMemcpy(result, idxBuffer, sizeof(UINT), cudaMemcpyDeviceToHost);
    return result[0];
}

template<class T>
T ReduceMaxT(T* value, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceMax << <iBlock, iThread >> > (value, iJump, count);
    }
    T result[1];
    cudaMemcpy(result, value, sizeof(T), cudaMemcpyDeviceToHost);
    return result[0];
}

template<class T>
T ReduceMinT(T* value, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceMin << <iBlock, iThread >> > (value, iJump, count);
    }
    T result[1];
    cudaMemcpy(result, value, sizeof(T), cudaMemcpyDeviceToHost);
    return result[0];
}

Real ReduceSum(Real* value, UINT count)
{
    return ReduceSumT(value, count);
}

extern QLAPI INT ReduceSum(INT* value, UINT count)
{
    return ReduceSumT(value, count);
}

extern QLAPI UINT ReduceSum(UINT* value, UINT count)
{
    return ReduceSumT(value, count);
}

extern QLAPI Real ReduceMin(Real* value, UINT count)
{
    return ReduceMinT(value, count);
}

extern QLAPI Real ReduceMax(Real* value, UINT count)
{
    return ReduceMaxT(value, count);
}

#if _QL_DOUBLEFLOAT
extern QLAPI FLOAT ReduceSum(FLOAT* value, UINT count)
{
    return ReduceSumT(value, count);
}
#else
extern QLAPI DOUBLE ReduceSum(DOUBLE* value, UINT count)
{
    return ReduceSumT(value, count);
}
#endif

extern QLAPI void ReduceSum(Real* res, Real* value, UINT count)
{
    ReduceSumT(res, value, count);
}

/**
* reduce sum
* 'value' will be changed, the first element is the result
*/
QLComplex ReduceSum(QLComplex* value, UINT count)
{
    const UINT iRequiredDim = (count + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        const UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;
        _kernelReduceComp << <iBlock, iThread >> > (value, iJump, count);
    }
    QLComplex result[1];
    cudaMemcpy(result, value, sizeof(QLComplex), cudaMemcpyDeviceToHost);
    return result[0];
}

/**
* sum _{if cndition[n] = conditionEqual} value[n * stride + offset]
* work space must has length >= (len(v) + 1) / 2
*/
Real ConditionalSum(const Real* value, BYTE byStride, BYTE offset, const BYTE* condition, BYTE conditionEqual, UINT count, Real* workSpace)
{
    UINT iThreadNeeded = Ceil(count, 2);
    UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;

    _kernelConditionalReduceReal << <iBlock, iThread >> > (value, workSpace, byStride, offset, count, condition, conditionEqual);

    return ReduceSum(workSpace, iThreadNeeded);
}

void ConditionalSum(Real* pDeviceRes, const Real* value, BYTE byStride, BYTE offset, const BYTE* condition, BYTE conditionEqual, UINT count, Real* workSpace)
{
    UINT iThreadNeeded = Ceil(count, 2);
    UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;

    _kernelConditionalReduceReal << <iBlock, iThread >> > (value, workSpace, byStride, offset, count, condition, conditionEqual);

    ReduceSum(pDeviceRes, workSpace, iThreadNeeded);
}

/**
* sum _{if cndition[n] = conditionEqual} 1
* work space must has length >= (len(condition) + 1) / 2
*/
UINT ConditionalCount(const INT* condition, INT conditionEqual, UINT count, UINT* workSpace)
{
    UINT iThreadNeeded = Ceil(count, 2);
    UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;

    _kernelConditionalCount << <iBlock, iThread >> > (workSpace, count, condition, conditionEqual);

    return ReduceSum(workSpace, iThreadNeeded);
}

UINT ConditionalCount(const BYTE* condition, BYTE conditionEqual, UINT count, UINT* workSpace)
{
    UINT iThreadNeeded = Ceil(count, 2);
    UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock) : iThreadNeeded;

    _kernelConditionalCount << <iBlock, iThread >> > (workSpace, count, condition, conditionEqual);

    return ReduceSum(workSpace, iThreadNeeded);
}

UINT Max(const INT* v, UINT* workspace, UINT count)
{
    return ReduceMaxT(v, workspace, count);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================