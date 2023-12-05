//=============================================================================
// FILENAME : QKMeans.cu
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [20/12/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

#pragma region kernel

__global__ void _QL_LAUNCH_BOUND
_kernelQKMADVectorForm(Real* resAbs, Real* resPhase, const Real* __restrict__ v, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        const QLComplex cmp = _make_cuComplex(v[2 * idx], v[2 * idx + 1]);
        resAbs[idx] = _cuCabsf(cmp);
        resPhase[idx] = __cuCargf(cmp);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelQKMADVectorFormC(Real* resAbs, Real* resPhase, const QLComplex* __restrict__ v, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        resAbs[idx] = _cuCabsf(v[idx]);
        resPhase[idx] = __cuCargf(v[idx]);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelQKMADPickVector(Real* targetAbs, Real* targetPhase, const Real* __restrict__ sourceAbs, const Real* __restrict__ sourcePhase, UINT uiDim, UINT uiMax, UINT uiTotal)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    
    if (idx < uiMax)
    {
        UINT toPick = _deviceRandomUI(threadIdx.x, uiTotal);
        for (UINT i = 0; i < uiDim; ++i)
        {
            targetAbs[uiDim * idx + i] = sourceAbs[uiDim * toPick + i];
            targetPhase[uiDim * idx + i] = sourcePhase[uiDim * toPick + i];
        }
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelQKMADPickVectorAndCheckDot(
    Real* targetAbs, Real* targetPhase, 
    QLComplex* dotres, 
    const Real* __restrict__ sourceAbs, 
    const Real* __restrict__ sourcePhase, 
    const Real* __restrict__ uAbs, 
    const Real* __restrict__ uPhase, 
    UINT uiDim, UINT uiMax, UINT uiTotal)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        UINT toPick = _deviceRandomUI(threadIdx.x, uiTotal);
        QLComplex c = _zeroc;
        for (UINT i = 0; i < uiDim; ++i)
        {
            const Real fAbs = sourceAbs[uiDim * toPick + i];
            const Real fPhase = sourcePhase[uiDim * toPick + i];
            targetAbs[uiDim * idx + i] = fAbs;
            targetPhase[uiDim * idx + i] = fPhase;

            const QLComplex c1 = _make_cuComplex(uAbs[i] * _cos(uPhase[i]), -uAbs[i] * _sin(uPhase[i]));
            const QLComplex c2 = _make_cuComplex(fAbs * _cos(fPhase), fAbs * _sin(fPhase));
            c = _cuCaddf(_cuCmulf(c1, c2), c);
        }
        dotres[idx] = _make_cuComplex(_cuCabsf(c), F(0.0));
    }
}

#if 0
/**
* resv = [v[0], v[1] + v[2] i, v[3] + v[4] i, v[5] + v[6] i]
*/
__global__ void
_QL_LAUNCH_BOUND
_kernelQKMADVectorToNormalizedVector(QLComplex* resv, const Real * __restrict__ v, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        resv[4 * idx] = _make_cuComplex(v[7 * idx], F(0.0));
        resv[4 * idx + 1] = _make_cuComplex(v[7 * idx + 1], v[7 * idx + 4]);
        resv[4 * idx + 2] = _make_cuComplex(v[7 * idx + 2], v[7 * idx + 5]);
        resv[4 * idx + 3] = _make_cuComplex(v[7 * idx + 3], v[7 * idx + 6]);
    }
}

/**
* 
*/
__global__ void
_QL_LAUNCH_BOUND
_kernelQKMADVectorToAngle(Real* y1, Real* y2, Real* z1, Real* z2, QLComplex* v, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        const Real ph1 = __cuCargf(v[4 * idx + 1]);
        const Real ph2 = __cuCargf(v[4 * idx + 2]);
        const Real ph3 = __cuCargf(v[4 * idx + 3]);
        const Real a3 = ph1 * F(0.5);
        const Real a4 = (ph3 - ph2) * F(0.5);
        const Real a6 = (ph2 + ph3 - ph1) * F(0.25);

        z1[idx] = a6 * F(-2.0);
        z2[2 * idx] = a3 * F(-2.0);
        z2[2 * idx + 1] = a4 * F(-2.0);

        const BYTE lengthPower = 2;
        Real YdegreeList[3];
        UINT degreelst = 0;
        #pragma unroll
        for (UINT i = 0; i < lengthPower; ++i)
        {
            const UINT degreeCount = 1U << (lengthPower - i - 1U);
            const UINT skip = 1U << i;
            const UINT stride = 1U << (i + 1U);
            for (UINT j = 0; j < degreeCount; ++j)
            {
                const Real cs = _cuCabsf(v[4 * idx + j * stride]);
                const Real sn = _cuCabsf(v[4 * idx + j * stride + skip]);
                const Real degree = _atan2(sn, cs);
                YdegreeList[degreelst] = degree;
                ++degreelst;
                const Real csd = _cos(degree);
                const Real snd = _sin(degree);
                if (csd > _QL_FLT_MIN_)
                {
                    for (UINT k = 0; k < skip; ++k)
                    {
                        v[4 * idx + j * stride + k] = cuCdivf_cr(v[4 * idx + j * stride + k], csd);
                    }
                }
                if (snd > _QL_FLT_MIN_)
                {
                    for (UINT k = 0; k < skip; ++k)
                    {
                        v[4 * idx + j * stride + skip + k] = cuCdivf_cr(v[4 * idx + j * stride + skip + k], snd);
                    }
                }
            }
        }

        y1[idx] = YdegreeList[2] * F(2.0);
        y2[2 * idx] = YdegreeList[0] * F(2.0);
        y2[2 * idx + 1] = YdegreeList[1] * F(2.0);
    }
}

template<class T, class Tc>
__global__ void
_QL_LAUNCH_BOUND
_kernelCenterListStep1(const T* __restrict__ arr, T* res, UINT* rescount, UINT uiMax, const Tc* __restrict__ condition, Tc condEqual)
{
    const UINT uiId = (threadIdx.x + blockIdx.x * blockDim.x);
    UINT uiId2 = (uiId << 1);
    if (condition[uiId2] == condEqual)
    {
        res[uiId * 7] = arr[uiId2 * 7];
        res[uiId * 7 + 1] = arr[uiId2 * 7 + 1];
        res[uiId * 7 + 2] = arr[uiId2 * 7 + 2];
        res[uiId * 7 + 3] = arr[uiId2 * 7 + 3];
        res[uiId * 7 + 4] = arr[uiId2 * 7 + 4];
        res[uiId * 7 + 5] = arr[uiId2 * 7 + 5];
        res[uiId * 7 + 6] = arr[uiId2 * 7 + 6];

        rescount[uiId] = 1;
    }
    else
    {
        res[uiId * 7] = static_cast<T>(0);
        res[uiId * 7 + 1] = static_cast<T>(0);
        res[uiId * 7 + 2] = static_cast<T>(0);
        res[uiId * 7 + 3] = static_cast<T>(0);
        res[uiId * 7 + 4] = static_cast<T>(0);
        res[uiId * 7 + 5] = static_cast<T>(0);
        res[uiId * 7 + 6] = static_cast<T>(0);

        rescount[uiId] = 0;
    }

    ++uiId2;
    if (uiId2 < uiMax && condition[uiId2] == condEqual)
    {
        res[uiId * 7] += arr[uiId2 * 7];
        res[uiId * 7 + 1] += arr[uiId2 * 7 + 1];
        res[uiId * 7 + 2] += arr[uiId2 * 7 + 2];
        res[uiId * 7 + 3] += arr[uiId2 * 7 + 3];
        res[uiId * 7 + 4] += arr[uiId2 * 7 + 4];
        res[uiId * 7 + 5] += arr[uiId2 * 7 + 5];
        res[uiId * 7 + 6] += arr[uiId2 * 7 + 6];

        ++rescount[uiId];
    }
}

template<class T>
__global__ void
_QL_LAUNCH_BOUND
_kernelCenterListStep2(T* arr, UINT* arrcount, UINT uiJump, UINT uiMax)
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
        #pragma unroll
        for (BYTE i = 0; i < 7; ++i)
        {
            arr[(uiIdFrom - uiJump) * 7 + i] += arr[uiIdFrom * 7 + i];
        }
        arrcount[uiIdFrom - uiJump] += arrcount[uiIdFrom];
    }
}

template<class T>
__global__ void
_QL_LAUNCH_BOUND_SINGLE
_kernelCenterListStep3(T* res, UINT* count, const T* __restrict__ arr, const UINT* __restrict__ arrcount)
{
    if (arrcount[0] > 0)
    {
        #pragma unroll
        for (BYTE i = 0; i < 7; ++i)
        {
            res[i] = arr[i] / arrcount[0];
        }
    }
    count[0] = arrcount[0];
}

//__global__ void
//_QL_LAUNCH_BOUND
//_kernelQKMeansInitialK(BYTE* kBuffer, UINT uiMax, BYTE kMax)
//{
//    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
//
//    if (idx < uiMax)
//    {
//        kBuffer[idx] = static_cast<BYTE>(__r->_deviceRandomI(threadIdx.x, kMax));
//    }
//}

__global__ void
_QL_LAUNCH_BOUND
_kernelQKMeansSpliteK(UINT* kBuffer, UINT uiMax, UINT kToSplit, UINT* splitTo, UINT countOfSplitTo)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax && kBuffer[idx] == kToSplit)
    {
        kBuffer[idx] = splitTo[__r->_deviceRandomI(threadIdx.x, countOfSplitTo)];
    }
}

__global__ void
_QL_LAUNCH_BOUND
_kernelQKMeansCalcDist(
    const UINT * __restrict__ kBuffer, 
    const QLComplex* __restrict__ points, 
    const QLComplex* __restrict__ centers,
    Real* res, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        const UINT pidx = kBuffer[idx];

        QLComplex resc = _cuCmulf(points[4 * idx], centers[4 * pidx]);
        resc = _cuCaddf(resc, _cuCmulf(_cuConjf(points[4 * idx + 1]), centers[4 * pidx + 1]));
        resc = _cuCaddf(resc, _cuCmulf(_cuConjf(points[4 * idx + 2]), centers[4 * pidx + 2]));
        resc = _cuCaddf(resc, _cuCmulf(_cuConjf(points[4 * idx + 3]), centers[4 * pidx + 3]));

        res[idx] = _cuCabsf(resc);
    }
}


#endif

#pragma endregion

QLQuantumKmeans::QLQuantumKmeans()
    : m_uiVectorDim(0)
    , m_uiReferenceVectorCount(0)
    , m_uiTestVectorCount(0)
    , m_pDeviceReferenceVectorsAbs(NULL)
    , m_pDeviceReferenceVectorsPhase(NULL)
    , m_pDeviceTestVectorsAbs(NULL)
    , m_pDeviceTestVectorsPhase(NULL)
{
    
}

void QLQuantumKmeans::LoadFile(const CCString& sReferenceCSV, const CCString& sTestCSV)
{

}

QLQuantumKmeans::~QLQuantumKmeans()
{
    if (NULL != m_pDeviceReferenceVectorsAbs)
    {
        checkCudaErrors(cudaFree(m_pDeviceReferenceVectorsAbs));
    }

    if (NULL != m_pDeviceReferenceVectorsAbs)
    {
        checkCudaErrors(cudaFree(m_pDeviceReferenceVectorsAbs));
    }

    if (NULL != m_pDeviceTestVectorsAbs)
    {
        checkCudaErrors(cudaFree(m_pDeviceTestVectorsAbs));
    }

    if (NULL != m_pDeviceTestVectorsPhase)
    {
        checkCudaErrors(cudaFree(m_pDeviceTestVectorsPhase));
    }
}

void QLQuantumKmeans::LoadFile(const CCString& sReferenceCSV, Real** targetAbs, Real** targetPhase, UINT& uiDim, UINT& uiCount)
{
    UINT w, h;
    TArray<Real> data = ReadCSVAR(sReferenceCSV, w, h);
    appGeneral(_T("CSV %s \n loaded %d x %d.\n"), sReferenceCSV.c_str(), w, h);

    if (0 != (w & 1))
    {
        appGeneral(_T("w must be power of 2\n"));
        return;
    }

    uiDim = w / 2;
    uiCount = h;
    Real* pDeviceVectors = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceVectors, sizeof(Real) * w * h));
    checkCudaErrors(cudaMalloc((void**)targetAbs, sizeof(Real) * uiDim * uiCount));
    checkCudaErrors(cudaMalloc((void**)targetPhase, sizeof(Real) * uiDim * uiCount));
    checkCudaErrors(cudaMemcpy(pDeviceVectors, data.GetData(), sizeof(Real) * w * h, cudaMemcpyHostToDevice));
    UINT elementCount = uiDim * uiCount;
    UINT uiBlock = Ceil(elementCount, _QL_LAUNCH_MAX_THREAD);
    _kernelQKMADVectorForm << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (*targetAbs, *targetPhase, pDeviceVectors, uiDim * uiCount);
    checkCudaErrors(cudaFree(pDeviceVectors));
}

void QLQuantumKmeans::TestCircuitBuildState(const CCString& sReferenceCSV, const CCString& sAmplitudeSave, const CCString& sMeasureRate, UINT vectorCount, UINT testRepeat)
{
    UINT uiVectorDim = 0;
    UINT uiReferenceVectorCount = 0;
    Real* pDeviceAbs = NULL;
    Real* pDevicePhase = NULL;
    LoadFile(sReferenceCSV, &pDeviceAbs, &pDevicePhase, uiVectorDim, uiReferenceVectorCount);
    UINT uiVectorPower = MostSignificantPowerTwo(uiVectorDim);
    UINT uiVectorCountPower = MostSignificantPowerTwo(vectorCount);

    Real* pDeviceWorkSpaceAbsU = NULL;
    Real* pDeviceWorkSpacePhaseU = NULL;
    QLComplex* pDeviceDotRes = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceAbsU, sizeof(Real) * uiVectorDim * testRepeat));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpacePhaseU, sizeof(Real) * uiVectorDim * testRepeat));
    checkCudaErrors(cudaMalloc((void**)&pDeviceDotRes, sizeof(QLComplex) * vectorCount));

    UINT uiBlock = Ceil(testRepeat, _QL_LAUNCH_MAX_THREAD);
    _kernelQKMADPickVector << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (pDeviceWorkSpaceAbsU, pDeviceWorkSpacePhaseU, pDeviceAbs, pDevicePhase, uiVectorDim, testRepeat, uiReferenceVectorCount);

    Real* pDeviceWorkSpaceAbsV = NULL;
    Real* pDeviceWorkSpacePhaseV = NULL;
    Real* pDeviceWorkSpaceY = NULL;
    Real* pDeviceWorkSpaceZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceAbsV, sizeof(Real) * uiVectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpacePhaseV, sizeof(Real) * uiVectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceY, sizeof(Real) * uiVectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceZ, sizeof(Real) * uiVectorDim * vectorCount));
    Real* pHostWorkSpaceY = reinterpret_cast<Real*>(malloc(sizeof(Real) * uiVectorDim * vectorCount));
    Real* pHostWorkSpaceZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * uiVectorDim * vectorCount));

    LONGLONG veclen = 1LL << (uiVectorPower + uiVectorCountPower);
    QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * veclen));

    TArray<BYTE> ampGate;
    TArray<BYTE> ampDaggerGate;
    for (BYTE byToAdd = 0; byToAdd < static_cast<BYTE>(uiVectorCountPower + uiVectorPower); ++byToAdd)
    {
        ampGate.AddItem(byToAdd);
        if (byToAdd < uiVectorPower)
        {
            ampDaggerGate.AddItem(byToAdd);
        }
    }
    TArray<BYTE> idxs;
    for (UINT i = 0; i < uiVectorPower + uiVectorCountPower; ++i)
    {
        if (i < uiVectorPower)
        {
            idxs.AddItem(static_cast<BYTE>(0));
            //idxs.AddItem(static_cast<BYTE>(2));
        }
        else
        {
            idxs.AddItem(static_cast<BYTE>(2));
        }
    }

    TArray<QLComplex> finalAbs;
    TArray<QLComplex> measureRate;

    appPushLogDate(FALSE);
    for (UINT i = 0; i < testRepeat; ++i)
    {
        uiBlock = Ceil(vectorCount, _QL_LAUNCH_MAX_THREAD);
        _kernelQKMADPickVectorAndCheckDot << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (
            pDeviceWorkSpaceAbsV, 
            pDeviceWorkSpacePhaseV,
            pDeviceDotRes,
            pDeviceAbs, 
            pDevicePhase, 
            pDeviceWorkSpaceAbsU + uiVectorDim * i, 
            pDeviceWorkSpacePhaseU + uiVectorDim * i,
            uiVectorDim, 
            vectorCount, 
            uiReferenceVectorCount);

        QLComplex* pHostDotRes = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * vectorCount));
        checkCudaErrors(cudaMemcpy(pHostDotRes, pDeviceDotRes, sizeof(QLComplex) * vectorCount, cudaMemcpyDeviceToHost));
        QLMatrix expectedRes(vectorCount, 1, pHostDotRes);

        //==================================
        // 1 - Build the amplitude encoder for V
        CalculateDegrees(pDeviceWorkSpaceAbsV, pDeviceWorkSpacePhaseV, vectorCount, uiVectorPower, pDeviceWorkSpaceY, pDeviceWorkSpaceZ);
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceY, pDeviceWorkSpaceY, sizeof(Real) * uiVectorDim * vectorCount, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceZ, pDeviceWorkSpaceZ, sizeof(Real) * uiVectorDim * vectorCount, cudaMemcpyDeviceToHost));
        QLGate ampGate1 = ExchangeToYZGate(uiVectorCountPower, uiVectorPower, pHostWorkSpaceY, pHostWorkSpaceZ, FALSE);

        // 2 - Build the amplitude encoder for U
        CalculateDegrees(pDeviceWorkSpaceAbsU + uiVectorDim * i, pDeviceWorkSpacePhaseU + uiVectorDim * i, 1, uiVectorPower, pDeviceWorkSpaceY, pDeviceWorkSpaceZ);
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceY, pDeviceWorkSpaceY, sizeof(Real) * uiVectorDim, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceZ, pDeviceWorkSpaceZ, sizeof(Real) * uiVectorDim, cudaMemcpyDeviceToHost));
        QLGate ampGate2 = ExchangeToYZGate(uiVectorPower, pHostWorkSpaceY, pHostWorkSpaceZ, FALSE);
        ampGate2.Dagger();

        // 3 - controlled collapse
        QLGate totalCircuit;
        totalCircuit.AddQubits(static_cast<BYTE>(uiVectorCountPower + uiVectorPower));

        totalCircuit.AppendGate(ampGate1, ampGate);
        totalCircuit.AppendGate(ampGate2, ampDaggerGate);
        QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
        for (UINT j = 0; j < uiVectorPower; ++j)
        {
            TArray<BYTE> bitsctr;
            bitsctr.AddItem(static_cast<BYTE>(j));
            totalCircuit.AppendGate(ctrCollapse, bitsctr);
        }

        // 4 - Simulate and get the measure rate
        TArray<SBasicOperation> ops = totalCircuit.GetOperation(totalCircuit.m_lstQubits);
        SIZE_T opssize = ops.Num();

        QuESTEnv evn = createQuESTEnv();
        Qureg vec = createQureg(totalCircuit.GetQubitCount(), evn);

        syncQuESTEnv(evn);
        vec.stateVec.real[0] = F(1.0);
        vec.stateVec.imag[0] = F(0.0);
        for (LONGLONG line2 = 1; line2 < veclen; ++line2)
        {
            vec.stateVec.real[line2] = F(0.0);
            vec.stateVec.imag[line2] = F(0.0);
        }
        copyStateToGPU(vec);

        Real probToBuild = F(1.0);
        for (SIZE_T j = 0; j < opssize; ++j)
        {
            probToBuild = probToBuild * QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(j)]);
        }
        measureRate.AddItem(_make_cuComplex(probToBuild, F(0.0)));
        syncQuESTEnv(evn);
        copyStateFromGPU(vec);
        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            res[line2].x = static_cast<Real>(vec.stateVec.real[line2]);
            res[line2].y = static_cast<Real>(vec.stateVec.imag[line2]);
        }
        QLMatrix finalstate = ShowStateVectorDetail(res, idxs, FALSE);
        finalstate.ElementAbs();
        finalstate.ReShape(uiVectorDim, vectorCount);
        finalstate.Print();

        // 5 - Compare the results
        Real norm = expectedRes.Norm2();
        expectedRes = expectedRes / norm;
        finalstate.Sub(expectedRes);
        Real diff = finalstate.Norm2();
        if (diff > F(0.001))
        {
            appGeneral(_T("difference is too large: %f\n"), diff);
        }

        destroyQureg(vec, evn);
        destroyQuESTEnv(evn);

        // 6 - Record the results
        finalAbs.Append(expectedRes.ToVector().GetData(), vectorCount);
        if (49U == (i % 50U))
        {
            appGeneral(_T("="));
            appPopLogDate();
            appGeneral(_T("%d\n"), i);
            appPushLogDate(FALSE);
        }
        else
        {
            appGeneral(_T("="));
        }
    }
    appPopLogDate();

    QLMatrix finalAmplitudes = QLMatrix::CopyCreate(vectorCount, testRepeat, finalAbs.GetData());
    SaveCSVR(finalAmplitudes, sAmplitudeSave);

    QLMatrix finalRates = QLMatrix::CopyCreate(1, testRepeat, measureRate.GetData());
    SaveCSVR(finalRates, sMeasureRate);

    checkCudaErrors(cudaFree(pDeviceWorkSpaceAbsU));
    checkCudaErrors(cudaFree(pDeviceWorkSpacePhaseU));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceAbsV));
    checkCudaErrors(cudaFree(pDeviceWorkSpacePhaseV));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceY));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceZ));

    appSafeFree(pHostWorkSpaceY);
    appSafeFree(pHostWorkSpaceZ);
    appSafeFree(res);
}

void QLQuantumKmeans::TestCircuitBuildStateOnce(const QLMatrix& hostVi, const QLMatrix& hostU, UINT vectorCount, UINT vectorDim)
{
    #pragma region build circuit
    //================= Build Circuit ===================
    UINT uiVectorPower = MostSignificantPowerTwo(vectorDim);
    UINT uiVectorCountPower = MostSignificantPowerTwo(vectorCount);

    QLComplex* deviceWorkSpaceComplex = NULL;
    Real* pDeviceWorkSpaceAbs = NULL;
    Real* pDeviceWorkSpacePhase = NULL;
    Real* pDeviceWorkSpaceY = NULL;
    Real* pDeviceWorkSpaceZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceWorkSpaceComplex, sizeof(QLComplex) * vectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceAbs, sizeof(Real) * vectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpacePhase, sizeof(Real) * vectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceY, sizeof(Real) * vectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceZ, sizeof(Real) * vectorDim * vectorCount));
    Real* pHostWorkSpaceY = reinterpret_cast<Real*>(malloc(sizeof(Real) * vectorDim * vectorCount));
    Real* pHostWorkSpaceZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * vectorDim * vectorCount));

    checkCudaErrors(cudaMemcpy(deviceWorkSpaceComplex, hostVi.HostBuffer(), sizeof(QLComplex) * vectorDim * vectorCount, cudaMemcpyHostToDevice));
    UINT elementCount = vectorDim * vectorCount;
    UINT uiBlock = Ceil(elementCount, _QL_LAUNCH_MAX_THREAD);
    _kernelQKMADVectorFormC << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (pDeviceWorkSpaceAbs, pDeviceWorkSpacePhase, deviceWorkSpaceComplex, vectorDim * vectorCount);
    CalculateDegrees(pDeviceWorkSpaceAbs, pDeviceWorkSpacePhase, vectorCount, uiVectorPower, pDeviceWorkSpaceY, pDeviceWorkSpaceZ);
    checkCudaErrors(cudaMemcpy(pHostWorkSpaceY, pDeviceWorkSpaceY, sizeof(Real) * vectorDim * vectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostWorkSpaceZ, pDeviceWorkSpaceZ, sizeof(Real) * vectorDim * vectorCount, cudaMemcpyDeviceToHost));
    QLGate ampGate1 = ExchangeToYZGate(uiVectorCountPower, uiVectorPower, pHostWorkSpaceY, pHostWorkSpaceZ, FALSE);

    checkCudaErrors(cudaMemcpy(deviceWorkSpaceComplex, hostU.HostBuffer(), sizeof(QLComplex) * vectorDim, cudaMemcpyHostToDevice));
    elementCount = vectorDim;
    uiBlock = Ceil(elementCount, _QL_LAUNCH_MAX_THREAD);
    _kernelQKMADVectorFormC << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (pDeviceWorkSpaceAbs, pDeviceWorkSpacePhase, deviceWorkSpaceComplex, vectorDim);
    CalculateDegrees(pDeviceWorkSpaceAbs, pDeviceWorkSpacePhase, 1, uiVectorPower, pDeviceWorkSpaceY, pDeviceWorkSpaceZ);
    checkCudaErrors(cudaMemcpy(pHostWorkSpaceY, pDeviceWorkSpaceY, sizeof(Real) * vectorDim, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostWorkSpaceZ, pDeviceWorkSpaceZ, sizeof(Real) * vectorDim, cudaMemcpyDeviceToHost));
    QLGate ampGate2 = ExchangeToYZGate(uiVectorPower, pHostWorkSpaceY, pHostWorkSpaceZ, FALSE);
    ampGate2.Dagger();

    QLGate totalCircuit;
    totalCircuit.AddQubits(static_cast<BYTE>(uiVectorCountPower + uiVectorPower));

    TArray<BYTE> ampGateQubits;
    TArray<BYTE> ampDaggerGateQubits;
    for (BYTE byToAdd = 0; byToAdd < static_cast<BYTE>(uiVectorCountPower + uiVectorPower); ++byToAdd)
    {
        ampGateQubits.AddItem(byToAdd);
        if (byToAdd < uiVectorPower)
        {
            ampDaggerGateQubits.AddItem(byToAdd);
        }
    }
    totalCircuit.AppendGate(ampGate1, ampGateQubits);
    totalCircuit.AppendGate(ampGate2, ampDaggerGateQubits);

    #pragma endregion

    //================= Add Amplitude amplifier ===================
    QLGate totalWithAA;
    UINT totalQubits = uiVectorCountPower + uiVectorPower + 1;
    //UINT totalQubits = uiVectorCountPower + uiVectorPower;
    totalWithAA.AddQubits(static_cast<BYTE>(totalQubits));
    TArray<BYTE> addAQubits;
    addAQubits.Append(ByteSequnce + 1, uiVectorCountPower + uiVectorPower);
    TArray<BYTE> sxQubits;
    sxQubits.Append(ByteSequnce + 1, uiVectorPower);
    sxQubits.AddItem(0);
    //sxQubits.AddItem(uiVectorCountPower + uiVectorPower + 1);
    QLGate s0Gate = GroverSXGate(static_cast<BYTE>(uiVectorPower), 0, FALSE);

    QLGate x(EBasicOperation::EBO_X);
    QLGate h(EBasicOperation::EBO_H);
    TArray<BYTE> ancillaQubit;
    ancillaQubit.AddItem(0);
    totalWithAA.AppendGate(totalCircuit, addAQubits);
    totalWithAA.AppendGate(x, ancillaQubit);
    totalWithAA.AppendGate(h, ancillaQubit);

    for (UINT i = 0; i < 1; ++i)
    {
        totalWithAA.AppendGate(s0Gate, sxQubits);
        totalCircuit.Dagger();
        totalWithAA.AppendGate(totalCircuit, addAQubits);
        totalWithAA.AppendGate(s0Gate, sxQubits);
        totalCircuit.Dagger();
        totalWithAA.AppendGate(totalCircuit, addAQubits);
    }

    totalWithAA.AppendGate(h, ancillaQubit);
    totalWithAA.AppendGate(x, ancillaQubit);

    //totalWithAA.DebugPrint(2);


    //================= Calculate Classical result ===================
    //The hostVi is:
    //v1x v2x
    //v1y v2y
    //v1z v2z
    //so we need to dagger it first
    QLMatrix hostViCopy = hostVi;
    hostViCopy.Dagger();
    hostViCopy = hostViCopy * hostU;
    Real norm = hostViCopy.Norm2();
    hostViCopy = hostViCopy / norm;
    hostViCopy.ReShape(vectorCount, 1);
    hostViCopy.ElementAbs();
    hostViCopy.Print("expected");


    //================= Compare Classical result =====================
    //QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
    //for (UINT j = 0; j < uiVectorPower; ++j)
    //{
    //    TArray<BYTE> bitsctr;
    //    bitsctr.AddItem(static_cast<BYTE>(j + 1));
    //    totalWithAA.AppendGate(ctrCollapse, bitsctr);
    //}
    LONGLONG veclen = 1LL << totalQubits;
    QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * veclen));
    TArray<SBasicOperation> ops = totalWithAA.GetOperation(totalWithAA.m_lstQubits);
    SIZE_T opssize = ops.Num();

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(totalWithAA.GetQubitCount(), evn);

    syncQuESTEnv(evn);
    vec.stateVec.real[0] = F(1.0);
    vec.stateVec.imag[0] = F(0.0);
    for (LONGLONG line2 = 1; line2 < veclen; ++line2)
    {
        vec.stateVec.real[line2] = F(0.0);
        vec.stateVec.imag[line2] = F(0.0);
    }
    copyStateToGPU(vec);

    Real probToBuild = F(1.0);
    for (SIZE_T i = 0; i < opssize; ++i)
    {
        probToBuild = probToBuild * QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);
    for (LONGLONG line2 = 0; line2 < veclen; ++line2)
    {
        res[line2].x = static_cast<Real>(vec.stateVec.real[line2]);
        res[line2].y = static_cast<Real>(vec.stateVec.imag[line2]);
    }
    TArray<BYTE> idxs;
    for (UINT i = 0; i < totalQubits; ++i)
    {
        if (i < uiVectorPower)
        {
            //idxs.AddItem(static_cast<BYTE>(0));
            idxs.AddItem(static_cast<BYTE>(2));
        }
        else
        {
            idxs.AddItem(static_cast<BYTE>(2));
        }
    }
    QLMatrix finalstate = ShowStateVectorDetail(res, idxs, FALSE);
    //finalstate.ElementAbs();
    finalstate.ReShape(vectorDim, vectorCount * 2);
    finalstate.Print("res");

    //================= Show Measure probability =====================
    appGeneral(_T("measure rate = %f\n"), probToBuild);

    //================= Free the resources =====================
    checkCudaErrors(cudaFree(deviceWorkSpaceComplex));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceAbs));
    checkCudaErrors(cudaFree(pDeviceWorkSpacePhase));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceY));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceZ));
    appSafeFree(pHostWorkSpaceY);
    appSafeFree(pHostWorkSpaceZ);
    appSafeFree(res);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================