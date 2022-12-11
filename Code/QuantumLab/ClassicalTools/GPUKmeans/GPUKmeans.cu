//=============================================================================
// FILENAME : GPUKmeans.cu
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [09/12/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void
_QL_LAUNCH_BOUND
_kernelInitialK(BYTE* kBuffer, UINT uiMax, BYTE kMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        kBuffer[idx] = static_cast<BYTE>(__r->_deviceRandomI(threadIdx.x, kMax));
    }
}

__global__ void
_QL_LAUNCH_BOUND
_kernelSpliteK(BYTE* kBuffer, UINT uiMax, BYTE kToSplit, BYTE* splitTo, BYTE countOfSplitTo)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax && kBuffer[idx] == kToSplit)
    {
        kBuffer[idx] = splitTo[__r->_deviceRandomI(threadIdx.x, countOfSplitTo)];
    }
}

__global__ void
_QL_LAUNCH_BOUND
_kernelReclassify(
    const FLOAT * __restrict__ pDataPoints,
    const FLOAT * const __restrict__ * __restrict__ pCenters, BYTE* kBuffer, UINT* countBuffer, UINT uiMax, BYTE kMax, BYTE kStride, BYTE kDim)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx >= uiMax)
    {
        return;
    }

    BYTE byMin = 0;
    FLOAT dist = 0.0f;
    for (BYTE kD = 0; kD < kDim; ++kD)
    {
        const UINT pIdx = idx * kStride + kD;
        const FLOAT fDelta = pDataPoints[pIdx] - pCenters[0][kD];
        dist += fDelta * fDelta;
    }

    for (BYTE k = 1; k < kMax; ++k)
    {
        FLOAT dist2 = 0.0f;
        for (BYTE kD = 0; kD < kDim; ++kD)
        {
            const UINT pIdx = idx * kStride + kD;
            const FLOAT fDelta = pDataPoints[pIdx] - pCenters[k][kD];
            dist2 += fDelta * fDelta;
        }

        if (dist2 < dist)
        {
            dist = dist2;
            byMin = k;
        }
    }

    if (byMin != kBuffer[idx])
    {
        kBuffer[idx] = byMin;
        countBuffer[idx] = 1;
    }
    else
    {
        countBuffer[idx] = 0;
    }
}

#pragma endregion

QLGPUKmeans::QLGPUKmeans(const CCString & fileName, BYTE dim, BYTE k, UBOOL bDebugChange)
    : m_bDebugChange(bDebugChange)
    , m_byK(k)
    , m_uiN(0)
    , m_byDim(dim)
    , m_byStride(dim)
    , m_uiBlock(1)
    , m_uiThread(1)
    , m_pDeviceData(NULL)
    , m_pDeviceK(NULL)
    , m_pHostK(NULL)
    , m_pDeviceDeviceCenters(NULL)
    , m_pDeviceFLOATWorkingBuffer(NULL)
    , m_pDeviceUINTWorkingBuffer(NULL)
    , m_pDeviceTempKBuffer(NULL)
{
    UINT w, h;
    TArray<FLOAT> data = ReadCSVAF(fileName, w, h);

    if (m_bDebugChange)
    {
        appGeneral(_T("File loaded.\n"));
    }

    m_byStride = static_cast<BYTE>(w);
    m_uiN = h;

    m_uiBlock = m_uiN > _QL_LAUNCH_MAX_THREAD ? Ceil(m_uiN, _QL_LAUNCH_MAX_THREAD) : 1;
    m_uiThread = m_uiN > _QL_LAUNCH_MAX_THREAD ? Ceil(m_uiN, m_uiBlock) : m_uiN;

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(FLOAT) * w * h));
    checkCudaErrors(cudaMemcpy(m_pDeviceData, data.GetData(), sizeof(FLOAT) * w * h, cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceK, sizeof(BYTE) * h));
    m_pHostK = (BYTE*)malloc(sizeof(BYTE) * h);

    for (BYTE byK = 0; byK < k; ++byK)
    {
        TArray<FLOAT> onecenter;
        for (BYTE byD = 0; byD < dim; ++byD)
        {
            onecenter.AddItem(0.0f);
        }
        m_pHostCenters.AddItem(onecenter);

        FLOAT* onecenterDevice;
        checkCudaErrors(cudaMalloc((void**)&onecenterDevice, sizeof(FLOAT) * dim));
        m_pDeviceCenters.AddItem(onecenterDevice);
    }

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDeviceCenters, sizeof(FLOAT*) * k));
    checkCudaErrors(cudaMemcpy(m_pDeviceDeviceCenters, m_pDeviceCenters.GetData(), sizeof(FLOAT*) * k, cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceFLOATWorkingBuffer, sizeof(FLOAT) * h));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceUINTWorkingBuffer, sizeof(UINT) * h));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceTempKBuffer, sizeof(BYTE) * k));
}

QLGPUKmeans::~QLGPUKmeans()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
    checkCudaErrors(cudaFree(m_pDeviceK));
    checkCudaErrors(cudaFree(m_pDeviceDeviceCenters));
    checkCudaErrors(cudaFree(m_pDeviceFLOATWorkingBuffer));
    checkCudaErrors(cudaFree(m_pDeviceUINTWorkingBuffer));
    checkCudaErrors(cudaFree(m_pDeviceTempKBuffer));

    for (BYTE k = 0; k < m_byK; ++k)
    {
        checkCudaErrors(cudaFree(m_pDeviceCenters[k]));
    }

    appSafeFree(m_pHostK);
}

void QLGPUKmeans::InitialK()
{
    _kernelInitialK << <m_uiBlock, m_uiThread >> > (m_pDeviceK, m_uiN, m_byK);

    UBOOL split = CalculateCountAndSplit();
    while (split)
    {
        split = CalculateCountAndSplit();
    }
}

UBOOL QLGPUKmeans::CalculateCountAndSplit()
{
    UINT iMax = 0;
    BYTE maxK = 0;
    TArray<BYTE> needToSplit;
    m_HostCount.RemoveAll();
    for (BYTE i = 0; i < m_byK; ++i)
    {
        UINT uiCount = ConditionalCount(m_pDeviceK, i, m_uiN, m_pDeviceUINTWorkingBuffer);
        if (0 == uiCount)
        {
            needToSplit.AddItem(i);
        }
        else if (uiCount > iMax)
        {
            maxK = i;
            iMax = uiCount;
        }
        m_HostCount.AddItem(uiCount);
    }

    if (needToSplit.Num() > 0)
    {
        if (m_bDebugChange)
        {
            appGeneral(_T("splited.\n"));
        }
        needToSplit.AddItem(maxK);

        checkCudaErrors(cudaMemcpy(m_pDeviceTempKBuffer, needToSplit.GetData(), sizeof(BYTE) * needToSplit.Num(), cudaMemcpyHostToDevice));

        _kernelSpliteK << <m_uiBlock, m_uiThread >> > (m_pDeviceK, m_uiN, maxK, m_pDeviceTempKBuffer, static_cast<BYTE>(needToSplit.Num()));
        return TRUE;
    }
    return FALSE;
}

void QLGPUKmeans::CalculateCenter()
{
    for (BYTE byK = 0; byK < m_byK; ++byK)
    {
        for (BYTE byDim = 0; byDim < m_byDim; ++byDim)
        {
            m_pHostCenters[byK][byDim] = ConditionalSum(m_pDeviceData, m_byStride, byDim, m_pDeviceK, byK, m_uiN, m_pDeviceFLOATWorkingBuffer) / m_HostCount[byK];
            //appGeneral(_T("%d-%d %f\n"), byK, byDim, m_pHostCenters[byK][byDim]);
        }

        checkCudaErrors(cudaMemcpy(m_pDeviceCenters[byK], m_pHostCenters[byK].GetData(), sizeof(FLOAT) * m_byDim, cudaMemcpyHostToDevice));
    }
}

UINT QLGPUKmeans::ReClassify()
{
    _kernelReclassify << <m_uiBlock, m_uiThread >> > (m_pDeviceData, m_pDeviceDeviceCenters, m_pDeviceK, m_pDeviceUINTWorkingBuffer, m_uiN, m_byK, m_byStride, m_byDim);
    return ReduceSum(m_pDeviceUINTWorkingBuffer, m_uiN);
}

void QLGPUKmeans::Build(UINT minChange)
{
    InitialK();
    CalculateCenter();
    UINT changed = ReClassify();
    while (changed > minChange)
    {
        UBOOL split = CalculateCountAndSplit();
        while (split)
        {
            split = CalculateCountAndSplit();
        }

        CalculateCenter();
        changed = ReClassify();
        if (m_bDebugChange)
        {
            appGeneral(_T("%d Changed.\n"), changed);
        }
    }

    checkCudaErrors(cudaMemcpy(m_pHostK, m_pDeviceK, sizeof(BYTE) * m_uiN, cudaMemcpyDeviceToHost));
    if (m_bDebugChange)
    {
        appGeneral(_T("Finished.\n"));
    }
}

void QLGPUKmeans::Save(const CCString& fileName)
{
    SaveCSVAB(m_pHostK, 1, m_uiN, fileName);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================