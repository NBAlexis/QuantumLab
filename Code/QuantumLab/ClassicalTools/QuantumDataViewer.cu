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

__global__ void _QL_LAUNCH_BOUND
_kernelStateVectorToAmp(QLComplex* buffer, const Real* __restrict__ real, const Real* __restrict__ imag, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (idx < uiMax)
    {
        buffer[idx] = _make_cuComplex(real[idx], imag[idx]);
    }
}


#pragma endregion


/**
* upLimit: the full list of upper
* lowerLimit: the full list of lower
* for example, at some qubit, 
*    if it want index = 1, it should be lower=1, upper=2, 
*    if it want index = 0, it should be lower=0, upper=1, 
*    if it want all, it should be lower=0, upper=2, 
* addedIndex: start position in stateData
* qubitIndex: the qubit this recursive is reading
*/
TArray<QLComplex> RecursivelyGetElement(const QLComplex* stateData, TArray<BYTE> upLimit, TArray<BYTE> lowerLimit, UINT addedIndex, BYTE qubitIndex)
{
    TArray<QLComplex> ret;
    //factor is stride
    UINT factor = 1U << (upLimit.Num() - 1 - qubitIndex);

    if (1 == factor)
    {
        //this is qubit 0
        for (BYTE i = lowerLimit[qubitIndex]; i < upLimit[qubitIndex]; ++i)
        {
            ret.AddItem(stateData[addedIndex + i]);
        }
        return ret;
    }

    for (BYTE i = lowerLimit[qubitIndex]; i < upLimit[qubitIndex]; ++i)
    {
        ret.Append(RecursivelyGetElement(stateData, upLimit, lowerLimit, addedIndex + i * factor, qubitIndex + 1));
    }
    return ret;
}

QLMatrix QLAPI ShowStateVectorDetail(const QLComplex* statedata, TArray<BYTE> index, UBOOL bNormalize)
{
    TArray<BYTE> indexUpperLimit;
    TArray<BYTE> indexLowerLimit;
    INT expectedCount = 1;
    for (INT i = 0; i < index.Num(); ++i)
    {
        if (2 == index[i])
        {
            indexUpperLimit.InsertAt(0, 2);
            indexLowerLimit.InsertAt(0, static_cast<BYTE>(0));
            expectedCount = (expectedCount << 1);
        }
        else
        {
            indexUpperLimit.InsertAt(0, index[i] + 1);
            indexLowerLimit.InsertAt(0, index[i]);
        }
    }

    TArray<QLComplex> res = RecursivelyGetElement(statedata, indexUpperLimit, indexLowerLimit, 0, 0);
    if (expectedCount != res.Num())
    {
        appCrucial(_T("bad thing happens in ShowStateVectorDetail!\n"));
        return QLMatrix();
    }
    if (bNormalize)
    {
        UINT len;
        res = NormalizeV(res, len);
    }
    QLMatrix ret = QLMatrix::CopyCreate(static_cast<UINT>(expectedCount), 1, res.GetData());
    return ret;
}

QLMatrix QLAPI ShowStateVectorDetail(const Qureg& vec, TArray<BYTE> index, UBOOL bNormalize)
{
    UINT uiLen = 1UL << static_cast<UINT>(index.Num());
    QLComplex* stateData = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex)* uiLen));

    for (UINT idx = 0; idx < uiLen; ++idx)
    {
        stateData[idx] = _make_cuComplex(vec.stateVec.real[idx], vec.stateVec.imag[idx]);
    }
    
    QLMatrix ret = ShowStateVectorDetail(stateData, index, bNormalize);
    appSafeFree(stateData);
    return ret;
}

QLMatrix QLAPI ShowStateVectorDetail(const QLComplex* statedata, BYTE count, BYTE idx0, ...)
{
    va_list arg;
    {
        va_start(arg, idx0);
        TArray<BYTE> idx;
        idx.AddItem(idx0);
        for (BYTE i = 1; i < count; ++i)
        {
            idx.AddItem(va_arg(arg, BYTE));
        }

        QLMatrix ret = ShowStateVectorDetail(statedata, idx);
        va_end(arg);

        return ret;
    }
}

QLMatrix QLAPI ShowStateVectorDetail(const Qureg& vec, BYTE count, BYTE idx0, ...)
{
    va_list arg;
    {
        va_start(arg, idx0);
        TArray<BYTE> idx;
        idx.AddItem(idx0);
        for (BYTE i = 1; i < count; ++i)
        {
            idx.AddItem(va_arg(arg, BYTE));
        }

        QLMatrix ret = ShowStateVectorDetail(vec, idx);
        va_end(arg);

        return ret;
    }
}

void QLAPI HostBufferViewer(const Real* buffer, UINT w, UINT h)
{
    appPushLogDate(FALSE);

    for (UINT y = 0; y < h; ++y)
    {
        for (UINT x = 0; x < w; ++x)
        {
            appGeneral(_T("%f\t"), buffer[y * w + x]);
        }
        appGeneral(_T("\n"));
    }

    appPopLogDate();
}

void QLAPI HostBufferViewer(const QLComplex* buffer, UINT w, UINT h)
{
    appPushLogDate(FALSE);

    for (UINT y = 0; y < h; ++y)
    {
        for (UINT x = 0; x < w; ++x)
        {
            appGeneral(_T("%s,\t"), appPrintComplex(buffer[y * w + x].x, buffer[y * w + x].y).c_str());
        }
        appGeneral(_T("\n"));
    }

    appPopLogDate();
}

void QLAPI DeviceBufferViewer(const Real* buffer, UINT w, UINT h)
{
    Real* pHostBuffer = reinterpret_cast<Real*>(malloc(sizeof(Real) * w * h));
    checkCudaErrors(cudaMemcpy(pHostBuffer, buffer, sizeof(Real) * w * h, cudaMemcpyDeviceToHost));

    HostBufferViewer(pHostBuffer, w, h);

    appSafeFree(pHostBuffer);
}

QLMatrix QLAPI StateToMatrix(const struct Qureg& vec)
{
    UINT veclen = 1UL << static_cast<UINT>(vec.numQubitsInStateVec);

    Real* pReal = NULL;
    Real* pImag = NULL;
    QLComplex* pAmp = NULL;
    checkCudaErrors(cudaMalloc((void**)&pReal, sizeof(Real) * veclen));
    checkCudaErrors(cudaMalloc((void**)&pImag, sizeof(Real) * veclen));
    checkCudaErrors(cudaMalloc((void**)&pAmp, sizeof(QLComplex) * veclen));
    checkCudaErrors(cudaMemcpy(pReal, vec.stateVec.real, sizeof(Real) * veclen, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(pImag, vec.stateVec.imag, sizeof(Real) * veclen, cudaMemcpyHostToDevice));
    
    UINT block = 0;
    UINT thread = 0;
    __DECOMPOSE(veclen, block, thread);
    _kernelStateVectorToAmp << <block, thread >> > (pAmp, pReal, pImag, veclen);

    QLComplex* pHostAmp = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * veclen));
    checkCudaErrors(cudaMemcpy(pHostAmp, pAmp, sizeof(QLComplex) * veclen, cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(pReal));
    checkCudaErrors(cudaFree(pImag));
    checkCudaErrors(cudaFree(pAmp));

    return QLMatrix(veclen, 1, pHostAmp);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================