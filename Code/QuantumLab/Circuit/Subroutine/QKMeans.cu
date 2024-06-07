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

/**
* generate random numbers
* uiMax = number of vectors + 1
*/
__global__ void _QL_LAUNCH_BOUND
_kernelQKMADRandomVector_Step1(QLComplex* buffer, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        buffer[idx] = _deviceRandomC(threadIdx.x);
    }
}

/**
* normalize
* 
* uiMax = number of vectors + 1
*/
__global__ void _QL_LAUNCH_BOUND
_kernelQKMADRandomVector_Step2(QLComplex* buffer, UINT uiDim, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        Real fLength = F(0.0);
        for (UINT i = 0; i < uiDim; ++i)
        {
            fLength += __cuCabsSqf(buffer[idx * uiDim + i]);
        }
        fLength = _sqrt(fLength);

        for (UINT i = 0; i < uiDim; ++i)
        {
            buffer[idx * uiDim + i] = cuCdivf_cr(buffer[idx * uiDim + i], fLength);
        }
    }
}

/**
* to abs and phase
* uiMax = number of vectors
*/
__global__ void _QL_LAUNCH_BOUND
_kernelQKMADRandomVector_Step3(
    Real* absBuffer,
    Real* phaseBuffer,
    QLComplex* dotres,
    const QLComplex * __restrict__ buffer, UINT uiDim, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        QLComplex dotresOne = _zeroc;
        for (UINT i = 0; i < uiDim; ++i)
        {
            absBuffer[idx * uiDim + i] = _cuCabsf(buffer[idx * uiDim + i]);
            phaseBuffer[idx * uiDim + i] = __cuCargf(buffer[idx * uiDim + i]);

            dotresOne = _cuCaddf(dotresOne, _cuCmulf(_cuConjf(buffer[idx * uiDim + i]), buffer[uiMax * uiDim + i]));
        }
        dotres[idx] = dotresOne;

        if (0 == idx)
        {
            for (UINT i = 0; i < uiDim; ++i)
            {
                absBuffer[uiMax * uiDim + i] = _cuCabsf(buffer[uiMax * uiDim + i]);
                phaseBuffer[uiMax * uiDim + i] = __cuCargf(buffer[uiMax * uiDim + i]);
            }
        }
    }
}



__global__ void _QL_LAUNCH_BOUND
_kernelKMeans2DReadFileXY(
    const Real* __restrict__ all,
    Real* x,
    Real* y,
    UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (idx < uiMax)
    {
        x[idx] = all[2 * idx];
        y[idx] = all[2 * idx + 1];
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelKMeans2DInitialCenter(
    Real* x,
    Real* y,
    Real minX, Real maxX, Real minY, Real maxY,
    UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (idx < uiMax)
    {
        x[idx] = minX + _deviceRandomF(threadIdx.x) * (maxX - minX);
        y[idx] = minY + _deviceRandomF(threadIdx.x) * (maxY - minY);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelKMeans2DPointToAbsPhase(
    const Real* __restrict__ x,
    const Real* __restrict__ y,
    Real minX, Real maxX, Real minY, Real maxY,
    Real* angleY, Real* angleZ, UINT uiMax)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (idx < uiMax)
    {
        angleY[idx] = (x[idx] - minX) / (maxX - minX) * PI;
        angleZ[idx] = -(y[idx] - minY) / (maxY - minY) * PI2;
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


void QLQuantumKmeans::LoadFile(const CCString& sReferenceCSV, Real** targetAbs, Real** targetPhase, UINT& uiDim, UINT& uiCount)
{
    UINT uiBlock = 0;
    UINT uiThread = 0;
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
    __DECOMPOSE(elementCount, uiBlock, uiThread);
    _kernelQKMADVectorForm << <uiBlock, uiThread >> > (*targetAbs, *targetPhase, pDeviceVectors, uiDim * uiCount);
    checkCudaErrors(cudaFree(pDeviceVectors));
}

/**
* load the file sReferenceCSV
* 
* repeat for "testRepeat" times
* 
* amplitude amplification is applied for aaRepeat times
* aaRepeat >= 1
* 
*/
void QLQuantumKmeans::TestCircuitBuildState(const CCString& sReferenceCSV, const CCString& sAmplitudeSave, const CCString& sMeasureRate, UINT vectorCount, UINT testRepeat, UINT aaRepeat)
{
    UINT uiBlock = 0;
    UINT uiThread = 0;
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

    __DECOMPOSE(testRepeat, uiBlock, uiThread);
    _kernelQKMADPickVector << <uiBlock, uiThread >> > (pDeviceWorkSpaceAbsU, pDeviceWorkSpacePhaseU, pDeviceAbs, pDevicePhase, uiVectorDim, testRepeat, uiReferenceVectorCount);

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

    LONGLONG veclen = 1LL << (uiVectorPower + uiVectorCountPower + 2);
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
        }
        else
        {
            idxs.AddItem(static_cast<BYTE>(2));
        }
    }
    idxs.AddItem(static_cast<BYTE>(0));
    idxs.AddItem(static_cast<BYTE>(0));
    TArray<BYTE> subspaceQubits;
    subspaceQubits.Append(ByteSequnce, uiVectorPower);

    TArray<QLComplex> finalAbs;
    TArray<QLComplex> measureRate;

    appPushLogDate(FALSE);
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(uiVectorPower + uiVectorCountPower + 2, evn);
    for (UINT i = 0; i < testRepeat; ++i)
    {
        __DECOMPOSE(vectorCount, uiBlock, uiThread);
        _kernelQKMADPickVectorAndCheckDot << <uiBlock, uiThread >> > (
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
        
        QLGate totalCircuit;
        totalCircuit.AddQubits(static_cast<BYTE>(uiVectorCountPower + uiVectorPower));

        totalCircuit.AppendGate(ampGate1, ampGate);
        totalCircuit.AppendGate(ampGate2, ampDaggerGate);

        // 3 - Amplitude amplification
        QLGate afterAA = AmplitudeAmplification(totalCircuit, subspaceQubits, 0, aaRepeat);

        // 4 - controlled collapse
        QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
        for (UINT j = 0; j < uiVectorPower; ++j)
        {
            TArray<BYTE> bitsctr;
            bitsctr.AddItem(static_cast<BYTE>(j));
            afterAA.AppendGate(ctrCollapse, bitsctr);
        }

        // 5 - Simulate and get the measure rate
        TArray<SBasicOperation> ops = afterAA.GetOperation(afterAA.m_lstQubits);
        SIZE_T opssize = ops.Num();

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
        //finalstate.Print();
        finalstate.ReShape(uiVectorDim, vectorCount);

        // 6 - Compare the results
        Real norm = expectedRes.Norm2();
        expectedRes = expectedRes / norm;
        finalstate.Sub(expectedRes);
        Real diff = finalstate.Norm2();
        if (diff > F(0.001))
        {
            appGeneral(_T("difference is too large: %f\n"), diff);
        }

        // 7 - Record the results
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
    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
    appPopLogDate();

    if (!sAmplitudeSave.IsEmpty())
    {
        QLMatrix finalAmplitudes = QLMatrix::CopyCreate(vectorCount, testRepeat, finalAbs.GetData());
        SaveCSVR(finalAmplitudes, sAmplitudeSave);
    }

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

TArray<QLComplex> QLQuantumKmeans::TestCircuitBuildState(UINT uiVectorDim, UINT vectorCount, UINT testRepeat, UINT aaRepeat)
{
    UINT uiBlock = 0;
    UINT uiThread = 0;
    TArray<QLComplex> ret;

    UINT uiVectorPower = MostSignificantPowerTwo(uiVectorDim);
    UINT uiVectorCountPower = MostSignificantPowerTwo(vectorCount);

    QLComplex* pDeviceComplexBuffer = NULL;
    Real* pDeviceWorkSpaceAbs = NULL;
    Real* pDeviceWorkSpacePhase = NULL;
    QLComplex* pDeviceDotRes = NULL;
    Real* pDeviceWorkSpaceY = NULL;
    Real* pDeviceWorkSpaceZ = NULL;

    checkCudaErrors(cudaMalloc((void**)&pDeviceComplexBuffer, sizeof(QLComplex) * uiVectorDim * (vectorCount + 1)));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceAbs, sizeof(Real) * uiVectorDim * (vectorCount + 1)));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpacePhase, sizeof(Real) * uiVectorDim * (vectorCount + 1)));
    checkCudaErrors(cudaMalloc((void**)&pDeviceDotRes, sizeof(QLComplex) * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceY, sizeof(Real) * uiVectorDim * vectorCount));
    checkCudaErrors(cudaMalloc((void**)&pDeviceWorkSpaceZ, sizeof(Real) * uiVectorDim * vectorCount));

    Real* pHostWorkSpaceY = reinterpret_cast<Real*>(malloc(sizeof(Real) * uiVectorDim * vectorCount));
    Real* pHostWorkSpaceZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * uiVectorDim * vectorCount));

    LONGLONG veclen = 1LL << (uiVectorPower + uiVectorCountPower + 2);
    //LONGLONG veclen2 = 1LL << (uiVectorPower + uiVectorCountPower);, this is to test original successful rate
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
        }
        else
        {
            idxs.AddItem(static_cast<BYTE>(2));
        }
    }
    idxs.AddItem(static_cast<BYTE>(0));
    idxs.AddItem(static_cast<BYTE>(0));
    TArray<BYTE> subspaceQubits;
    subspaceQubits.Append(ByteSequnce, uiVectorPower);
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(uiVectorPower + uiVectorCountPower + 2, evn);
    Qureg vec2 = createQureg(uiVectorPower + uiVectorCountPower, evn);

    for (UINT i = 0; i < testRepeat; ++i)
    {
        //==============================
        // 1 - generate random vectors
        __DECOMPOSE((vectorCount + 1) * uiVectorDim, uiBlock, uiThread);
        _kernelQKMADRandomVector_Step1 << <uiBlock, uiThread >> > (
            pDeviceComplexBuffer, 
            (vectorCount + 1) * uiVectorDim);

        __DECOMPOSE(vectorCount + 1, uiBlock, uiThread);
        _kernelQKMADRandomVector_Step2 << <uiBlock, uiThread >> > (
            pDeviceComplexBuffer, 
            uiVectorDim, 
            vectorCount + 1);

        __DECOMPOSE(vectorCount, uiBlock, uiThread);
        _kernelQKMADRandomVector_Step3 << <uiBlock, uiThread >> > (
            pDeviceWorkSpaceAbs, pDeviceWorkSpacePhase, pDeviceDotRes,
            pDeviceComplexBuffer, uiVectorDim, vectorCount);

        QLComplex* pHostDotRes = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * vectorCount));
        checkCudaErrors(cudaMemcpy(pHostDotRes, pDeviceDotRes, sizeof(QLComplex) * vectorCount, cudaMemcpyDeviceToHost));
        QLMatrix expectedRes(vectorCount, 1, pHostDotRes);

        //QLComplex* pHostShow = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * uiVectorDim * (vectorCount + 1)));
        //checkCudaErrors(cudaMemcpy(pHostShow, pDeviceComplexBuffer, sizeof(QLComplex) * uiVectorDim * (vectorCount + 1), cudaMemcpyDeviceToHost));
        //QLMatrix showmtr(vectorCount + 1, uiVectorDim, pHostShow);
        //showmtr.Transpose();
        //showmtr.Print("v");

        //==================================
        // 1 - Build the amplitude encoder for V
        CalculateDegrees(pDeviceWorkSpaceAbs, pDeviceWorkSpacePhase, vectorCount, uiVectorPower, pDeviceWorkSpaceY, pDeviceWorkSpaceZ);
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceY, pDeviceWorkSpaceY, sizeof(Real) * uiVectorDim * vectorCount, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceZ, pDeviceWorkSpaceZ, sizeof(Real) * uiVectorDim * vectorCount, cudaMemcpyDeviceToHost));
        QLGate ampGate1 = ExchangeToYZGate(uiVectorCountPower, uiVectorPower, pHostWorkSpaceY, pHostWorkSpaceZ, FALSE);

        // 2 - Build the amplitude encoder for U
        CalculateDegrees(pDeviceWorkSpaceAbs + uiVectorDim * vectorCount, pDeviceWorkSpacePhase + uiVectorDim * vectorCount, 1, uiVectorPower, pDeviceWorkSpaceY, pDeviceWorkSpaceZ);
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceY, pDeviceWorkSpaceY, sizeof(Real) * uiVectorDim, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(pHostWorkSpaceZ, pDeviceWorkSpaceZ, sizeof(Real) * uiVectorDim, cudaMemcpyDeviceToHost));
        QLGate ampGate2 = ExchangeToYZGate(uiVectorPower, pHostWorkSpaceY, pHostWorkSpaceZ, FALSE);
        ampGate2.Dagger();

        QLGate totalCircuit;
        totalCircuit.AddQubits(static_cast<BYTE>(uiVectorCountPower + uiVectorPower));

        totalCircuit.AppendGate(ampGate1, ampGate);
        totalCircuit.AppendGate(ampGate2, ampDaggerGate);

        // 3 - Amplitude amplification
        QLGate afterAA = AmplitudeAmplification(totalCircuit, subspaceQubits, 0, aaRepeat);

        // 4 - controlled collapse
        QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
        for (UINT j = 0; j < uiVectorPower; ++j)
        {
            TArray<BYTE> bitsctr;
            bitsctr.AddItem(static_cast<BYTE>(j));
            afterAA.AppendGate(ctrCollapse, bitsctr);
        }

        // 5 - Simulate and get the measure rate
        TArray<SBasicOperation> ops = afterAA.GetOperation(afterAA.m_lstQubits);
        SIZE_T opssize = ops.Num();

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
        ret.AddItem(_make_cuComplex(probToBuild, F(0.0)));
        syncQuESTEnv(evn);
        copyStateFromGPU(vec);
        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            res[line2].x = static_cast<Real>(vec.stateVec.real[line2]);
            res[line2].y = static_cast<Real>(vec.stateVec.imag[line2]);
        }
        QLMatrix finalstate = ShowStateVectorDetail(res, idxs, FALSE);
        finalstate.ElementAbs();
        //finalstate.Print();
        finalstate.ReShape(uiVectorDim, vectorCount);

        // 6 - Compare the results
        // This has been tested, therefore not this
        //expectedRes.ElementAbs();
        //Real norm = expectedRes.Norm2();
        //expectedRes = expectedRes / norm;
        //QLMatrix diffmtr = finalstate - expectedRes;
        //Real diff = diffmtr.Norm2();
        //if (diff > F(0.001))
        //{
        //    appGeneral(_T("difference is too large: %f\n"), diff);
        //    expectedRes.Print("expect");
        //    finalstate.Print("finalstate");
        //}

        //=======================================
        //find out build rate before AA
        //This has been tested, therefore not this
        //for (UINT j = 0; j < uiVectorPower; ++j)
        //{
        //    TArray<BYTE> bitsctr;
        //    bitsctr.AddItem(static_cast<BYTE>(j));
        //    totalCircuit.AppendGate(ctrCollapse, bitsctr);
        //}
        //syncQuESTEnv(evn);
        //vec2.stateVec.real[0] = F(1.0);
        //vec2.stateVec.imag[0] = F(0.0);
        //for (LONGLONG line2 = 1; line2 < veclen2; ++line2)
        //{
        //    vec2.stateVec.real[line2] = F(0.0);
        //    vec2.stateVec.imag[line2] = F(0.0);
        //}
        //copyStateToGPU(vec2);
        //Real probToBuildBeforeAA = F(1.0);
        //TArray<SBasicOperation> ops2 = totalCircuit.GetOperation(totalCircuit.m_lstQubits);
        //SIZE_T opssize2 = ops2.Num();
        //for (SIZE_T j = 0; j < opssize2; ++j)
        //{
        //    probToBuildBeforeAA = probToBuildBeforeAA * QLGate::PerformBasicOperation(vec2, ops2[static_cast<INT>(j)]);
        //}

        // 7 - Record the results
        //appGeneral(_T("dim %d, count %d, repeat %d, aa %d, error %f, build rate %f, build rate before AA %f\n"), 
        //    uiVectorDim, vectorCount, i, aaRepeat, diff, probToBuild, probToBuildBeforeAA);
        //appGeneral(_T("dim %d, count %d, repeat %d, aa %d, error %f, build rate %f\n"), 
        //    uiVectorDim, vectorCount, i, aaRepeat, diff, probToBuild);
        appGeneral(_T("dim %d, count %d, repeat %d, aa %d, build rate %f\n"),
            uiVectorDim, vectorCount, i, aaRepeat, probToBuild);
    }

    destroyQureg(vec2, evn);
    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    checkCudaErrors(cudaFree(pDeviceComplexBuffer));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceAbs));
    checkCudaErrors(cudaFree(pDeviceWorkSpacePhase));
    checkCudaErrors(cudaFree(pDeviceDotRes));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceY));
    checkCudaErrors(cudaFree(pDeviceWorkSpaceZ));

    appSafeFree(pHostWorkSpaceY);
    appSafeFree(pHostWorkSpaceZ);
    appSafeFree(res);

    return ret;
}

/**
* 
*/
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
    TArray<BYTE> qubitsSubspace;
    qubitsSubspace.Append(ByteSequnce, uiVectorPower);

    QLGate afterAA = AmplitudeAmplification(totalCircuit, qubitsSubspace, 0, 1);

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
    appGeneral(_T("norm = %f\n"), norm * norm / 4);


    //================= Compare Classical result =====================
    QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
    for (UINT j = 0; j < uiVectorPower; ++j)
    {
        TArray<BYTE> bitsctr;
        bitsctr.AddItem(static_cast<BYTE>(j));
        afterAA.AppendGate(ctrCollapse, bitsctr);
    }
    //LONGLONG veclen = 1LL << (afterAA.GetQubitCount());
    LONGLONG veclen = 1LL << (totalCircuit.GetQubitCount());
    QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * veclen));
    //TArray<SBasicOperation> ops = afterAA.GetOperation(afterAA.m_lstQubits);
    for (UINT j = 0; j < uiVectorPower; ++j)
    {
        totalCircuit.AppendGate(ctrCollapse, static_cast<BYTE>(j));
    }
    TArray<SBasicOperation> ops = totalCircuit.GetOperation(afterAA.m_lstQubits);
    SIZE_T opssize = ops.Num();

    QuESTEnv evn = createQuESTEnv();
    //Qureg vec = createQureg(afterAA.GetQubitCount(), evn);
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
    //for (BYTE i = 0; i < afterAA.GetQubitCount(); ++i)
    for (BYTE i = 0; i < totalCircuit.GetQubitCount(); ++i)
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
    QLMatrix finalstate = ShowStateVectorDetail(res, idxs, FALSE);
    finalstate.ElementAbs();
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

/**
* 
*/
void QLQuantumKmeans::Kmeans2D(
    const CCString& sPoints, 
    const CCString& sSaveK, 
    const CCString& sSaveCenter, 
    const CCString& sRepeat,
    BYTE k, UINT iteration, UINT uiMinHit)
{
    UINT uiBlock = 0;
    UINT uiThread = 0;
    CCString kfileName = _T("");
    CCString cfileName = _T("");
    CCString rfileName = _T("");
    UINT w, h;
    TArray<Real> orignalPointsArray = ReadCSVAR(sPoints, w, h);

    Real* pOrignalX = NULL;
    Real* pOrignalY = NULL;
    Real* pWorkingSpaceOrignal = NULL;
    //UINT* pWorkingSpaceCount = NULL;
    checkCudaErrors(cudaMalloc((void**)&pWorkingSpaceOrignal, sizeof(Real) * h * 2));
    //checkCudaErrors(cudaMalloc((void**)&pWorkingSpaceCount, sizeof(UINT) * h));
    checkCudaErrors(cudaMalloc((void**)&pOrignalX, sizeof(Real) * h));
    checkCudaErrors(cudaMalloc((void**)&pOrignalY, sizeof(Real) * h));

    checkCudaErrors(cudaMemcpy(pWorkingSpaceOrignal, orignalPointsArray.GetData(), sizeof(Real) * h * 2, cudaMemcpyHostToDevice));
    __DECOMPOSE(h, uiBlock, uiThread);
    _kernelKMeans2DReadFileXY << <uiBlock, uiThread >> > (pWorkingSpaceOrignal, pOrignalX, pOrignalY, h);

    checkCudaErrors(cudaMemcpy(pWorkingSpaceOrignal, pOrignalX, sizeof(Real) * h, cudaMemcpyDeviceToDevice));
    const Real minX = ReduceMin(pWorkingSpaceOrignal, h);
    checkCudaErrors(cudaMemcpy(pWorkingSpaceOrignal, pOrignalX, sizeof(Real) * h, cudaMemcpyDeviceToDevice));
    const Real maxX = ReduceMax(pWorkingSpaceOrignal, h);

    checkCudaErrors(cudaMemcpy(pWorkingSpaceOrignal, pOrignalY, sizeof(Real) * h, cudaMemcpyDeviceToDevice));
    const Real minY = ReduceMin(pWorkingSpaceOrignal, h);
    checkCudaErrors(cudaMemcpy(pWorkingSpaceOrignal, pOrignalY, sizeof(Real) * h, cudaMemcpyDeviceToDevice));
    const Real maxY = ReduceMax(pWorkingSpaceOrignal, h);

    Real* pAngleY = NULL;
    Real* pAngleZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&pAngleY, sizeof(Real) * h));
    checkCudaErrors(cudaMalloc((void**)&pAngleZ, sizeof(Real) * h));
    _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pOrignalX, pOrignalY, minX, maxX, minY, maxY, pAngleY, pAngleZ, h);
    Real* pHostAngleY = reinterpret_cast<Real*>(malloc(sizeof(Real) * h));
    Real* pHostAngleZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * h));

    checkCudaErrors(cudaMemcpy(pHostAngleY, pAngleY, sizeof(Real) * h, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleZ, pAngleZ, sizeof(Real) * h, cudaMemcpyDeviceToHost));

    Real* pCenterX = NULL;
    Real* pCenterY = NULL;
    Real* pCenterAngleY = NULL;
    Real* pCenterAngleZ = NULL;
    //BYTE* pDeviceKBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&pCenterX, sizeof(Real) * k));
    checkCudaErrors(cudaMalloc((void**)&pCenterY, sizeof(Real) * k));
    checkCudaErrors(cudaMalloc((void**)&pCenterAngleY, sizeof(Real) * k));
    checkCudaErrors(cudaMalloc((void**)&pCenterAngleZ, sizeof(Real) * k));
    //checkCudaErrors(cudaMalloc((void**)&pDeviceKBuffer, sizeof(BYTE) * h));
    BYTE* pHostKBuffer = reinterpret_cast<BYTE*>(malloc(sizeof(BYTE) * h));
    Real* pHostCenterAngleY = reinterpret_cast<Real*>(malloc(sizeof(Real) * k));
    Real* pHostCenterAngleZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * k));
    Real* pHostCenterXY = reinterpret_cast<Real*>(malloc(sizeof(Real) * k * 2));

    Real* pHostCenterXBuffer = reinterpret_cast<Real*>(malloc(sizeof(Real) * k));
    Real* pHostCenterYBuffer = reinterpret_cast<Real*>(malloc(sizeof(Real) * k));
    INT* pHostCenterCountBuffer = reinterpret_cast<INT*>(malloc(sizeof(INT) * k));

    __DECOMPOSE(k, uiBlock, uiThread);
    _kernelKMeans2DInitialCenter << <uiBlock, uiThread >> > (pCenterX, pCenterY, minX, maxX, minY, maxY, k);
    checkCudaErrors(cudaMemcpy(pHostCenterXY, pCenterX, sizeof(Real) * k, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostCenterXY + k, pCenterY, sizeof(Real) * k, cudaMemcpyDeviceToHost));
    _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pCenterX, pCenterY, minX, maxX, minY, maxY, pCenterAngleY, pCenterAngleZ, k);
    checkCudaErrors(cudaMemcpy(pHostCenterAngleY, pCenterAngleY, sizeof(Real) * k, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostCenterAngleZ, pCenterAngleZ, sizeof(Real) * k, cudaMemcpyDeviceToHost));

    UINT kPower = MostSignificantPowerTwo(k);
    UINT veclen = 1UL << (kPower + 1);
    QLGate cc(EBasicOperation::EBO_CC, 0);
    QLGate had(EBasicOperation::EBO_H);
    TArray<BYTE> measurebits;
    for (BYTE i = 0; i < kPower; ++i)
    {
        measurebits.AddItem(i + 1);
    }
    TArray<BYTE> frybits;
    for (BYTE i = 0; i < kPower + 1; ++i)
    {
        frybits.AddItem(static_cast<BYTE>(kPower - i));
    }

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(kPower + 1, evn);

    INT* hitCount = reinterpret_cast<INT*>(malloc(sizeof(INT) * k));
    TArray<UINT> repeatsForIteration;

    //Iterations start
    for (UINT itr = 0; itr < iteration; ++itr)
    {
        //export host center
        appGeneral(_T("start iteration: %d\n"), itr);

        cfileName.Format(_T("%s_%d.csv"), sSaveCenter.c_str(), itr);
        SaveCSVAR(pHostCenterXY, k, 2, cfileName);

        //appGeneral(_T("x=%f, y=%f"), pHostCenterXY[0], pHostCenterXY[k]);
        UINT repeat = 0;
        memset(pHostCenterXBuffer, 0, sizeof(Real) * k);
        memset(pHostCenterYBuffer, 0, sizeof(Real) * k);
        memset(pHostCenterCountBuffer, 0, sizeof(INT) * k);

        //calculate K for each vector
        for (UINT uiV = 0; uiV < h; ++uiV)
        {
            QLGate wholeCircuit;
            wholeCircuit.AddQubits(static_cast<BYTE>(kPower + 1));
            for (BYTE hadIdx = 0; hadIdx < kPower; ++hadIdx)
            {
                wholeCircuit.AppendGate(had, hadIdx + 1);
            }
            QLGate amplitudeVi = FRyz(pHostCenterAngleY, pHostCenterAngleZ, kPower + 1);
            wholeCircuit.AppendGate(amplitudeVi, frybits);

            //ry(theta) rz(-phi)
            //and dagger, which is rz(phi), ry(-theta)
            QLGate rz(EBasicOperation::EBO_RZ, -pHostAngleZ[uiV]);
            QLGate ry(EBasicOperation::EBO_RY, pHostAngleY[uiV]);
            QLGate amp2;
            amp2.AddQubits(1);
            amp2.AppendGate(ry, 0);
            amp2.AppendGate(rz, 0);
            amp2.Dagger();
            
            wholeCircuit.AppendGate(amp2, 0);
            if (uiMinHit <= 1)
            {
                wholeCircuit.AppendGate(cc, 0);
            }

            syncQuESTEnv(evn);
            memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
            memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
            vec.stateVec.real[0] = F(1.0);
            vec.stateVec.imag[0] = F(0.0);
            copyStateToGPU(vec);

            TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
            SIZE_T opssize = ops.Num();
            for (SIZE_T i = 0; i < opssize; ++i)
            {
                QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
            }
            syncQuESTEnv(evn);

            if (uiMinHit > 1)
            {
                copyStateFromGPU(vec);
                while (0 != measure(vec, 0))
                {
                    copyStateToGPU(vec);
                    ++repeat;
                }
            }

#pragma region test state-vector

            //copyStateFromGPU(vec);
            //for k=4
            //QLMatrix show = ShowStateVectorDetail(vec, 3, 0, 2, 2);
            //show.ElementAbs();
            //show.Print();

#pragma endregion

            if (uiMinHit > 1)
            {
                copyStateFromGPU(vec);
            }
            
            UINT measureRes = 0;
            for (INT j = 0; j < measurebits.Num(); ++j)
            {
                INT out = measure(vec, measurebits[j]);
                if (1 == out)
                {
                    measureRes = measureRes | (1U << j);
                }
            }

            if (uiMinHit > 1)
            {
                memset(hitCount, 0, sizeof(INT) * k);
                hitCount[measureRes] = 1;
                UBOOL bHasHit = FALSE;
                while (!bHasHit)
                {
                    copyStateToGPU(vec);
                    ++repeat;
                    measureRes = 0;
                    for (INT j = 0; j < measurebits.Num(); ++j)
                    {
                        INT out = measure(vec, measurebits[j]);
                        if (1 == out)
                        {
                            measureRes = measureRes | (1U << j);
                        }
                    }
                    hitCount[measureRes] = hitCount[measureRes] + 1;
                    if (hitCount[measureRes] >= static_cast<INT>(uiMinHit))
                    {
                        bHasHit = TRUE;
                        pHostKBuffer[uiV] = static_cast<BYTE>(measureRes);
                    }
                }
            }
            else
            {
                pHostKBuffer[uiV] = static_cast<BYTE>(measureRes);
            }

            pHostCenterXBuffer[measureRes] = pHostCenterXBuffer[measureRes] + orignalPointsArray[2 * uiV];
            pHostCenterYBuffer[measureRes] = pHostCenterYBuffer[measureRes] + orignalPointsArray[2 * uiV + 1];
            pHostCenterCountBuffer[measureRes] = pHostCenterCountBuffer[measureRes] + 1;
        }

        repeatsForIteration.AddItem(repeat);

        //calculate center again
        //checkCudaErrors(cudaMemcpy(pDeviceKBuffer, pHostKBuffer, sizeof(BYTE) * h, cudaMemcpyHostToDevice));
        for (BYTE centerK = 0; centerK < k; ++centerK)
        {
            //Real sumx = ConditionalSum(pOrignalX, 1, 0, pDeviceKBuffer, centerK, h, pWorkingSpaceOrignal);
            //Real sumy = ConditionalSum(pOrignalY, 1, 0, pDeviceKBuffer, centerK, h, pWorkingSpaceOrignal);
            //UINT count = ConditionalCount(pDeviceKBuffer, centerK, h, pWorkingSpaceCount);
            //pHostCenterXY[centerK] = sumx / count;
            //pHostCenterXY[k + centerK] = sumy / count;

            pHostCenterXY[centerK] = pHostCenterXBuffer[centerK] / pHostCenterCountBuffer[centerK];
            pHostCenterXY[k + centerK] = pHostCenterYBuffer[centerK] / pHostCenterCountBuffer[centerK];

            //appGeneral(_T("calc res = %f, %f  vs %f, %f\n"), pHostCenterXY[centerK], pHostCenterXY[k + centerK], pHostCenterXBuffer[centerK], pHostCenterYBuffer[centerK]);
        }
        
        //export K file
        kfileName.Format(_T("%s_%d.csv"), sSaveK.c_str(), itr);
        SaveCSVAB(pHostKBuffer, 1, h, kfileName);

        //calculate angles for new center
        __DECOMPOSE(k, uiBlock, uiThread);
        checkCudaErrors(cudaMemcpy(pCenterX, pHostCenterXY, sizeof(Real) * k, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(pCenterY, pHostCenterXY + k, sizeof(Real) * k, cudaMemcpyHostToDevice));
        _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pCenterX, pCenterY, minX, maxX, minY, maxY, pCenterAngleY, pCenterAngleZ, k);
        checkCudaErrors(cudaMemcpy(pHostCenterAngleY, pCenterAngleY, sizeof(Real) * k, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(pHostCenterAngleZ, pCenterAngleZ, sizeof(Real) * k, cudaMemcpyDeviceToHost));
    }

    //final export host center
    cfileName.Format(_T("%s_final.csv"), sSaveCenter.c_str());
    SaveCSVAR(pHostCenterXY, 2, k, cfileName);

    if (uiMinHit > 1)
    {
        rfileName.Format(_T("%s.csv"), sRepeat.c_str());
        SaveCSVAUI(repeatsForIteration.GetData(), 1, iteration, rfileName);
    }

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    checkCudaErrors(cudaFree(pAngleY));
    checkCudaErrors(cudaFree(pAngleZ));
    appSafeFree(pHostAngleY);
    appSafeFree(pHostAngleZ);
    appSafeFree(hitCount);

    appSafeFree(pHostCenterXBuffer);
    appSafeFree(pHostCenterYBuffer);
    appSafeFree(pHostCenterCountBuffer);

    checkCudaErrors(cudaFree(pCenterX));
    checkCudaErrors(cudaFree(pCenterY));
    checkCudaErrors(cudaFree(pCenterAngleY));
    checkCudaErrors(cudaFree(pCenterAngleZ));
    //checkCudaErrors(cudaFree(pDeviceKBuffer));
    appSafeFree(pHostKBuffer);
    appSafeFree(pHostCenterAngleY);
    appSafeFree(pHostCenterAngleZ);
    appSafeFree(pHostCenterXY);

    checkCudaErrors(cudaFree(pOrignalX));
    checkCudaErrors(cudaFree(pOrignalY));
    checkCudaErrors(cudaFree(pWorkingSpaceOrignal));
    //checkCudaErrors(cudaFree(pWorkingSpaceCount));
}

void QLQuantumKmeans::KNN2D(const CCString& sTrainingPoints, const CCString& sTestPoints,
    const CCString& sLoadK, const CCString& sSaveK, const CCString& sRepeat,
    UINT kHit, UINT uiMaxCluster)
{
    UINT uiBlock = 0;
    UINT uiThread = 0;
    UINT w1, h1;
    UINT w2, h2;
    TArray<Real> orignalPointsArrayTraining = ReadCSVAR(sTrainingPoints, w1, h1);
    TArray<Real> orignalPointsArrayTesting = ReadCSVAR(sTestPoints, w2, h2);
    TArray<INT> orignalPointsCluster = ReadCSVAI(sLoadK, w1, h1);

    Real* pTrainingX = NULL;
    Real* pTrainingY = NULL;
    Real* pTestingX = NULL;
    Real* pTestingY = NULL;
    Real* pWorkingSpaceTraining = NULL;
    Real* pWorkingSpaceTesting = NULL;

    checkCudaErrors(cudaMalloc((void**)&pWorkingSpaceTraining, sizeof(Real) * h1 * 2));
    checkCudaErrors(cudaMalloc((void**)&pWorkingSpaceTesting, sizeof(Real) * h2 * 2));
    checkCudaErrors(cudaMalloc((void**)&pTrainingX, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pTrainingY, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pTestingX, sizeof(Real) * h2));
    checkCudaErrors(cudaMalloc((void**)&pTestingY, sizeof(Real) * h2));

    checkCudaErrors(cudaMemcpy(pWorkingSpaceTraining, orignalPointsArrayTraining.GetData(), sizeof(Real) * h1 * 2, cudaMemcpyHostToDevice));
    __DECOMPOSE(h1, uiBlock, uiThread);
    _kernelKMeans2DReadFileXY << <uiBlock, uiThread >> > (pWorkingSpaceTraining, pTrainingX, pTrainingY, h1);
    checkCudaErrors(cudaMemcpy(pWorkingSpaceTesting, orignalPointsArrayTesting.GetData(), sizeof(Real) * h2 * 2, cudaMemcpyHostToDevice));
    __DECOMPOSE(h2, uiBlock, uiThread);
    _kernelKMeans2DReadFileXY << <uiBlock, uiThread >> > (pWorkingSpaceTesting, pTestingX, pTestingY, h2);

    checkCudaErrors(cudaMemcpy(pWorkingSpaceTraining, pTrainingX, sizeof(Real) * h1, cudaMemcpyDeviceToDevice));
    Real minX = ReduceMin(pWorkingSpaceTraining, h1);
    checkCudaErrors(cudaMemcpy(pWorkingSpaceTraining, pTrainingX, sizeof(Real) * h1, cudaMemcpyDeviceToDevice));
    Real maxX = ReduceMax(pWorkingSpaceTraining, h1);

    checkCudaErrors(cudaMemcpy(pWorkingSpaceTraining, pTrainingY, sizeof(Real) * h1, cudaMemcpyDeviceToDevice));
    Real minY = ReduceMin(pWorkingSpaceTraining, h1);
    checkCudaErrors(cudaMemcpy(pWorkingSpaceTraining, pTrainingY, sizeof(Real) * h1, cudaMemcpyDeviceToDevice));
    Real maxY = ReduceMax(pWorkingSpaceTraining, h1);

    checkCudaErrors(cudaMemcpy(pWorkingSpaceTesting, pTestingX, sizeof(Real) * h2, cudaMemcpyDeviceToDevice));
    const Real minX2 = ReduceMin(pWorkingSpaceTesting, h2);
    minX = minX2 < minX ? minX2 : minX;
    checkCudaErrors(cudaMemcpy(pWorkingSpaceTesting, pTestingX, sizeof(Real) * h2, cudaMemcpyDeviceToDevice));
    const Real maxX2 = ReduceMax(pWorkingSpaceTesting, h2);
    maxX = maxX2 > maxX ? maxX2 : maxX;

    checkCudaErrors(cudaMemcpy(pWorkingSpaceTesting, pTestingY, sizeof(Real) * h2, cudaMemcpyDeviceToDevice));
    const Real minY2 = ReduceMin(pWorkingSpaceTesting, h2);
    minY = minY2 < minY ? minY2 : minY;
    checkCudaErrors(cudaMemcpy(pWorkingSpaceTesting, pTestingY, sizeof(Real) * h2, cudaMemcpyDeviceToDevice));
    const Real maxY2 = ReduceMax(pWorkingSpaceTesting, h2);
    maxY = maxY2 > maxY ? maxY2 : maxY;

    Real* pAngleYTraining = NULL;
    Real* pAngleZTraining = NULL;
    Real* pAngleYTesting = NULL;
    Real* pAngleZTesting = NULL;
    checkCudaErrors(cudaMalloc((void**)&pAngleYTraining, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pAngleZTraining, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pAngleYTesting, sizeof(Real) * h2));
    checkCudaErrors(cudaMalloc((void**)&pAngleZTesting, sizeof(Real) * h2));

    __DECOMPOSE(h1, uiBlock, uiThread);
    _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pTrainingX, pTrainingY, 
        minX, maxX, minY, maxY, pAngleYTraining, pAngleZTraining, h1);
    __DECOMPOSE(h2, uiBlock, uiThread);
    _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pTestingX, pTestingY, 
        minX, maxX, minY, maxY, pAngleYTesting, pAngleZTesting, h2);

    Real* pHostAngleYTraining = reinterpret_cast<Real*>(malloc(sizeof(Real) * h1));
    Real* pHostAngleZTraining = reinterpret_cast<Real*>(malloc(sizeof(Real) * h1));
    Real* pHostAngleYTesting = reinterpret_cast<Real*>(malloc(sizeof(Real) * h2));
    Real* pHostAngleZTesting = reinterpret_cast<Real*>(malloc(sizeof(Real) * h2));

    checkCudaErrors(cudaMemcpy(pHostAngleYTraining, pAngleYTraining, sizeof(Real) * h1, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleZTraining, pAngleZTraining, sizeof(Real) * h1, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleYTesting, pAngleYTesting, sizeof(Real) * h2, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleZTesting, pAngleZTesting, sizeof(Real) * h2, cudaMemcpyDeviceToHost));

    UINT kPower = MostSignificantPowerTwo(h1);
    UINT veclen = 1UL << (kPower + 1);
    //QLGate cc(EBasicOperation::EBO_CC, 0);
    QLGate had(EBasicOperation::EBO_H);
    TArray<BYTE> measurebits;
    for (BYTE i = 0; i < kPower; ++i)
    {
        measurebits.AddItem(i + 1);
    }
    TArray<BYTE> frybits;
    for (BYTE i = 0; i < kPower + 1; ++i)
    {
        frybits.AddItem(static_cast<BYTE>(kPower - i));
    }

    //INT* pointHitCount = reinterpret_cast<INT*>(malloc(sizeof(INT) * h1));
    INT* hitCount = reinterpret_cast<INT*>(malloc(sizeof(INT) * uiMaxCluster));

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(kPower + 1, evn);

    TArray<INT> clusterAssignment;
    TArray<UINT> repeats;
    //Iterations start
    for (UINT uiV = 0; uiV < h2; ++uiV)
    {
        QLGate wholeCircuit;
        wholeCircuit.AddQubits(static_cast<BYTE>(kPower + 1));
        for (BYTE hadIdx = 0; hadIdx < kPower; ++hadIdx)
        {
            wholeCircuit.AppendGate(had, hadIdx + 1);
        }
        QLGate amplitudeVi = FRyz(pHostAngleYTraining, pHostAngleZTraining, kPower + 1);
        wholeCircuit.AppendGate(amplitudeVi, frybits);

        //ry(theta) rz(-phi)
        //and dagger, which is rz(phi), ry(-theta)
        QLGate rz(EBasicOperation::EBO_RZ, -pHostAngleZTesting[uiV]);
        QLGate ry(EBasicOperation::EBO_RY, pHostAngleYTesting[uiV]);
        QLGate amp2;
        amp2.AddQubits(1);
        amp2.AppendGate(ry, 0);
        amp2.AppendGate(rz, 0);
        amp2.Dagger();

        wholeCircuit.AppendGate(amp2, 0);
        //wholeCircuit.AppendGate(cc, 0);

        syncQuESTEnv(evn);
        memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
        memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
        vec.stateVec.real[0] = F(1.0);
        vec.stateVec.imag[0] = F(0.0);
        copyStateToGPU(vec);

        TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
        SIZE_T opssize = ops.Num();
        for (SIZE_T i = 0; i < opssize; ++i)
        {
            QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        }
        syncQuESTEnv(evn);

        copyStateFromGPU(vec);

#pragma region test state-vector

        //if (0 == uiV)
        //{
        //    TArray<BYTE> getout;
        //    getout.AddItem(0);
        //    for (BYTE toshowbit = 0; toshowbit < kPower; ++toshowbit)
        //    {
        //        getout.AddItem(2);
        //    }
        //    QLMatrix show = ShowStateVectorDetail(vec, getout);
        //    show.Mul(_make_cuComplex(F(23.34282757393504), F(0.0)));
        //    show.Print();
        //}

#pragma endregion

        //memset(pointHitCount, 0, sizeof(INT) * h1);
        memset(hitCount, 0, sizeof(INT) * uiMaxCluster);

        UBOOL bHasHit = FALSE;
        UINT repeat = 0;
        while (!bHasHit)
        {
            ++repeat;
            if (0 == measure(vec, 0))
            {
                UINT measureRes = 0;
                for (INT j = 0; j < measurebits.Num(); ++j)
                {
                    INT out = measure(vec, measurebits[j]);
                    if (1 == out)
                    {
                        measureRes = measureRes | (1U << j);
                    }
                }

                UINT thisHit = orignalPointsCluster[measureRes];
                hitCount[thisHit] = hitCount[thisHit] + 1;
                if (hitCount[thisHit] >= static_cast<INT>(kHit))
                {
                    bHasHit = TRUE;
                    clusterAssignment.AddItem(thisHit);
                    repeats.AddItem(repeat);
                }
            }

            if (!bHasHit)
            {
                copyStateToGPU(vec);
            }
        }

        appGeneral(_T("finish %d / %d ...\n"), uiV, h2);
    }

    //final export host center
    SaveCSVAI(clusterAssignment.GetData(), 1, h2, sSaveK);

    SaveCSVAUI(repeats.GetData(), 1, h2, sRepeat);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    checkCudaErrors(cudaFree(pAngleYTraining));
    checkCudaErrors(cudaFree(pAngleZTraining));
    checkCudaErrors(cudaFree(pAngleYTesting));
    checkCudaErrors(cudaFree(pAngleZTesting));

    appSafeFree(pHostAngleYTraining);
    appSafeFree(pHostAngleZTraining);
    appSafeFree(pHostAngleYTesting);
    appSafeFree(pHostAngleZTesting);

    //appSafeFree(pointHitCount);
    appSafeFree(hitCount);

    checkCudaErrors(cudaFree(pTrainingX));
    checkCudaErrors(cudaFree(pTrainingY));
    checkCudaErrors(cudaFree(pTestingX));
    checkCudaErrors(cudaFree(pTestingY));
    checkCudaErrors(cudaFree(pWorkingSpaceTraining));
    checkCudaErrors(cudaFree(pWorkingSpaceTesting));
}

void QLQuantumKmeans::KNN3D(const CCString& sTrainingPoints, const CCString& sTestPoints,
    const CCString& sLoadK, const CCString& sSaveK, const CCString& sRepeat,
    UINT kHit, UINT uiMaxCluster)
{
    //UINT uiBlock = 0;
    //UINT uiThread = 0;
    UINT w1, h1;
    UINT w2, h2;
    TArray<Real> orignalPointsArrayTraining = ReadCSVAR(sTrainingPoints, w1, h1);
    TArray<Real> orignalPointsArrayTesting = ReadCSVAR(sTestPoints, w2, h2);
    TArray<INT> orignalPointsCluster = ReadCSVAI(sLoadK, w1, h1);

    Real* pTrainingPoints = NULL;
    Real* pTrainingY = NULL;
    Real* pTestingPoints = NULL;
    Real* pTestingY = NULL;
    Real* pHostTrainingY = reinterpret_cast<Real*>(malloc(sizeof(Real) * h1 * 4));
    Real* pHostTestY = reinterpret_cast<Real*>(malloc(sizeof(Real) * h2 * 4));

    checkCudaErrors(cudaMalloc((void**)&pTrainingPoints, sizeof(Real) * h1 * 4));
    checkCudaErrors(cudaMalloc((void**)&pTrainingY, sizeof(Real) * h1 * 4));
    checkCudaErrors(cudaMalloc((void**)&pTestingPoints, sizeof(Real) * h2 * 4));
    checkCudaErrors(cudaMalloc((void**)&pTestingY, sizeof(Real) * h2 * 4));

    checkCudaErrors(cudaMemcpy(pTrainingPoints, orignalPointsArrayTraining.GetData(), sizeof(Real) * h1 * 4, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(pTestingPoints, orignalPointsArrayTesting.GetData(), sizeof(Real) * h2 * 4, cudaMemcpyHostToDevice));

    CalculateDegreesReal(pTrainingPoints, h1, 2, pTrainingY);
    CalculateDegreesRealForEach(pTestingPoints, h2, 2, pTestingY);
    checkCudaErrors(cudaMemcpy(pHostTrainingY, pTrainingY, sizeof(Real) * h1 * 4, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostTestY, pTestingY, sizeof(Real) * h2 * 4, cudaMemcpyDeviceToHost));

    UINT kPower = MostSignificantPowerTwo(h1);
    UINT veclen = 1UL << (kPower + 4);
    TArray<BYTE> measurebits;
    for (BYTE i = 0; i < kPower; ++i)
    {
        measurebits.AddItem(i + 2);
    }

    //INT* pointHitCount = reinterpret_cast<INT*>(malloc(sizeof(INT) * h1));
    INT* hitCount = reinterpret_cast<INT*>(malloc(sizeof(INT) * uiMaxCluster));

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(kPower + 4, evn);

    TArray<INT> clusterAssignment;
    TArray<UINT> repeats;
    //Iterations start
    for (UINT uiV = 0; uiV < h2; ++uiV)
    {
        QLGate wholeCircuit;
        wholeCircuit.AddQubits(static_cast<BYTE>(kPower + 2));
        QLGate amplitudeVi = ExchangeToYGate(kPower, 2, pHostTrainingY);
        wholeCircuit.AppendGate(amplitudeVi, amplitudeVi.m_lstQubits);

        //ry(theta) rz(-phi)
        //and dagger, which is rz(phi), ry(-theta)
        QLGate amplitudeU = ExchangeToYGate(2, pHostTestY + uiV * 4);
        amplitudeU.Dagger();

        wholeCircuit.AppendGate(amplitudeU, 0, 1);
        TArray<BYTE> subspace;
        subspace.AddItem(0);
        subspace.AddItem(1);

        QLGate afterAA = AmplitudeAmplification(wholeCircuit, subspace, 0, 1);

        syncQuESTEnv(evn);
        memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
        memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
        vec.stateVec.real[0] = F(1.0);
        vec.stateVec.imag[0] = F(0.0);
        copyStateToGPU(vec);

        TArray<SBasicOperation> ops = afterAA.GetOperation(afterAA.m_lstQubits);
        SIZE_T opssize = ops.Num();
        Real fbuildRate = F(1.0);
        for (SIZE_T i = 0; i < opssize; ++i)
        {
            fbuildRate = fbuildRate * QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        }
        syncQuESTEnv(evn);

        copyStateFromGPU(vec);

#pragma region test state-vector

        //if (0 == uiV)
        //{
        //    TArray<BYTE> getout;
        //    getout.AddItem(0);
        //    for (BYTE toshowbit = 0; toshowbit < kPower; ++toshowbit)
        //    {
        //        getout.AddItem(2);
        //    }
        //    QLMatrix show = ShowStateVectorDetail(vec, getout);
        //    show.Mul(_make_cuComplex(F(23.34282757393504), F(0.0)));
        //    show.Print();
        //}

#pragma endregion

        UBOOL bHasHit = FALSE;
        UINT repeat = 0;
        memset(hitCount, 0, sizeof(INT) * uiMaxCluster);
        while (!bHasHit)
        {
            ++repeat;
            if (0 == measure(vec, 0))
            {
                UINT measureRes = 0;
                for (INT j = 0; j < measurebits.Num(); ++j)
                {
                    INT out = measure(vec, measurebits[j]);
                    if (1 == out)
                    {
                        measureRes = measureRes | (1U << j);
                    }
                }

                //pointHitCount[measureRes] = pointHitCount[measureRes] + 1;
                //if (1 == pointHit || pointHitCount[measureRes] == pointHit)
                //{
                    UINT thisHit = orignalPointsCluster[measureRes];
                    hitCount[thisHit] = hitCount[thisHit] + 1;
                    if (hitCount[thisHit] >= static_cast<INT>(kHit))
                    {
                        bHasHit = TRUE;
                        clusterAssignment.AddItem(thisHit);
                        repeats.AddItem(repeat);
                    }
                //}
            }

            if (!bHasHit)
            {
                copyStateToGPU(vec);
            }
        }

        appGeneral(_T("finish %d / %d ...\n"), uiV, h2);
    }

    //final export host center
    SaveCSVAI(clusterAssignment.GetData(), 1, h2, sSaveK);

    SaveCSVAUI(repeats.GetData(), 1, h2, sRepeat);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    appSafeFree(pHostTrainingY);
    appSafeFree(pHostTestY);

    //appSafeFree(pointHitCount);
    appSafeFree(hitCount);

    checkCudaErrors(cudaFree(pTrainingY));
    checkCudaErrors(cudaFree(pTestingY));
    checkCudaErrors(cudaFree(pTrainingPoints));
    checkCudaErrors(cudaFree(pTestingPoints));
}


void QLQuantumKmeans::KNN2DAnsatz(const CCString& sAnsatz, const CCString& sTestPoints, const CCString& sSaveK, const CCString& sRepeat, UINT kHit)
{
    //UINT uiBlock = 0;
    //UINT uiThread = 0;
    UINT w, h;
    TArray<Real> ansatzParam = ReadCSVAR(sAnsatz, w, h);
    TArray<QLComplex> orignalPointsArray = ReadCSVA(sTestPoints, w, h);

    CTwoLocal ansatz(9, 1, ESingleLayer::RYRZ, ELinkLayer::CZ, ELinkStyle::SCA);
    ansatz.SetParameters(ansatzParam);
    QLGate ansatzGate = ansatz.BuildStateWithParam();
    TArray<UINT> repeats;

    UINT veclen = 512;
    
    TArray<INT> clusterAssignment;
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(9, evn);

    for (UINT uiV = 0; uiV < h; ++uiV)
    {
        UINT hitCount[4] = { 0,0,0,0 };
        QLGate vectorAE = AmplitudeEncodeOneVector(orignalPointsArray.GetData() + 2 * uiV, 1, FALSE);
        QLGate wholeCircuit;
        wholeCircuit.AddQubits(9);
        wholeCircuit.AppendGate(ansatzGate, ansatzGate.m_lstQubits);
        vectorAE.Dagger();
        wholeCircuit.AppendGate(vectorAE, 0);

        syncQuESTEnv(evn);
        memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
        memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
        vec.stateVec.real[0] = F(1.0);
        vec.stateVec.imag[0] = F(0.0);
        copyStateToGPU(vec);

        TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
        SIZE_T opssize = ops.Num();
        for (SIZE_T i = 0; i < opssize; ++i)
        {
            QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        }
        syncQuESTEnv(evn);

        copyStateFromGPU(vec);

        UBOOL bHasHit = FALSE;
        UINT repeat = 0;
        while (!bHasHit)
        {
            ++repeat;
            if (0 == measure(vec, 0))
            {
                //UINT measureRes = 0;
                UINT thisHit = static_cast<UINT>(measure(vec, 7) + (measure(vec, 8) << 1));
                hitCount[thisHit] = hitCount[thisHit] + 1;
                if (hitCount[thisHit] >= static_cast<INT>(kHit))
                {
                    bHasHit = TRUE;
                    clusterAssignment.AddItem(thisHit);
                    repeats.AddItem(repeat);
                }
            }

            if (!bHasHit)
            {
                copyStateToGPU(vec);
            }
        }

        appGeneral(_T("finish %d / %d ...\n"), uiV, h);
    }

    //final export host center
    SaveCSVAI(clusterAssignment.GetData(), 1, h, sSaveK);
    SaveCSVAUI(repeats.GetData(), 1, h, sRepeat);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
}

void QLQuantumKmeans::KNNAnsatz(
    const CCString& sAnsatz, 
    const CCString& sTestPoints, 
    const CCString& sScore,
    BYTE ansatzQubits, 
    BYTE byMeasureQubits,
    UBOOL bAdaptive, UINT uiAnsatzLevel, UINT uiRepeat)
{
    //UINT uiBlock = 0;
    //UINT uiThread = 0;
    UINT w, h, ha;
    TArray<Real> ansatzParam = ReadCSVAR(sAnsatz, w, ha);
    TArray<QLComplex> orignalPointsArray = ReadCSVA(sTestPoints, w, h);
    BYTE vectorPower = static_cast<BYTE>(MostSignificantPowerTwo(w));

    QLComplex* pDeviceVectors = NULL;
    Real* pDeviceAbs = NULL;
    Real* pDevicePhase = NULL;
    Real* pDeviceY = NULL;
    Real* pDeviceZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceVectors, sizeof(QLComplex) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDeviceAbs, sizeof(Real) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDevicePhase, sizeof(Real) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDeviceY, sizeof(Real) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDeviceZ, sizeof(Real) * w * h));
    checkCudaErrors(cudaMemcpy(pDeviceVectors, orignalPointsArray.GetData(), sizeof(QLComplex) * w * h, cudaMemcpyHostToDevice));
    CalculateDegreesForEach(pDeviceVectors, pDeviceAbs, pDevicePhase, h, vectorPower, pDeviceY, pDeviceZ);
    Real* pHostY = reinterpret_cast<Real*>(malloc(sizeof(Real) * w * h));
    Real* pHostZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * w * h));
    checkCudaErrors(cudaMemcpy(pHostY, pDeviceY, sizeof(Real) * w * h, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostZ, pDeviceZ, sizeof(Real) * w * h, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(pDeviceVectors));
    checkCudaErrors(cudaFree(pDeviceAbs));
    checkCudaErrors(cudaFree(pDevicePhase));
    checkCudaErrors(cudaFree(pDeviceY));
    checkCudaErrors(cudaFree(pDeviceZ));

    QLGate ansatzGate;
    if (bAdaptive)
    {
        uiAnsatzLevel = (ha / (ansatzQubits << 1U)) - 1;
        CTwoLocalAdaptive ansatz(ansatzQubits, ESingleLayer::RYRZ, ELinkLayer::CZ, ELinkStyle::Circular);
        ansatz.SetMaxLayer(1000);
        for (UINT i = 0; i < uiAnsatzLevel; ++i)
        {
            ansatz.IncreaseAdaptive();
        }
        ansatz.SetParameters(ansatzParam);
        ansatzGate = ansatz.BuildStateWithParam();
    }
    else
    {
        CTwoLocal ansatz(ansatzQubits, uiAnsatzLevel, ESingleLayer::RYRZ, ELinkLayer::CZ, ELinkStyle::SCA);
        ansatz.SetParameters(ansatzParam);
        ansatzGate = ansatz.BuildStateWithParam();
    }

    UINT uiPossibleRes = 1 << byMeasureQubits;
    Real* scores = reinterpret_cast<Real*>(malloc(sizeof(Real) * h * uiPossibleRes));

    UINT veclen = 1UL << static_cast<UINT>(ansatzQubits);
    QLGate cc(EBasicOperation::EBO_CC, 0);
    
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(ansatzQubits, evn);
    QLGate wholeCircuitBefore;
    wholeCircuitBefore.AddQubits(ansatzQubits);
    wholeCircuitBefore.AppendGate(ansatzGate, ansatzGate.m_lstQubits);
    syncQuESTEnv(evn);
    memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
    memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
    vec.stateVec.real[0] = F(1.0);
    vec.stateVec.imag[0] = F(0.0);
    copyStateToGPU(vec);

    TArray<SBasicOperation> opsBefore = wholeCircuitBefore.GetOperation(wholeCircuitBefore.m_lstQubits);
    SIZE_T opssizeBefore = opsBefore.Num();
    for (SIZE_T i = 0; i < opssizeBefore; ++i)
    {
        QLGate::PerformBasicOperation(vec, opsBefore[static_cast<INT>(i)]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    Real* res = reinterpret_cast<Real*>(appAlloca(sizeof(Real) * uiPossibleRes));
    INT* rescount = reinterpret_cast<INT*>(appAlloca(sizeof(INT) * uiPossibleRes));
    TArray<INT> qubitsToSee;
    for (BYTE bym = 0; bym < byMeasureQubits; ++bym)
    {
        qubitsToSee.AddItem(static_cast<INT>(ansatzQubits - 1 - bym));
    }

    //Real* realpart = reinterpret_cast<Real*>(malloc(sizeof(Real) * veclen));
    //Real* imagpart = reinterpret_cast<Real*>(malloc(sizeof(Real) * veclen));
    //memcpy(realpart, vec.stateVec.real, sizeof(Real) * veclen);
    //memcpy(imagpart, vec.stateVec.imag, sizeof(Real) * veclen);
    for (UINT uiV = 0; uiV < h; ++uiV)
    {
        //QLGate vectorAE = AmplitudeEncodeOneVector(orignalPointsArray.GetData() + w * uiV, vectorPower, FALSE);
        QLGate vectorAE = ExchangeToYZGate(vectorPower, pHostY + w * uiV, pHostZ + w * uiV, FALSE);
        QLGate wholeCircuit;
        wholeCircuit.AddQubits(ansatzQubits);
        vectorAE.Dagger();
        wholeCircuit.AppendGate(vectorAE, vectorAE.m_lstQubits);
        if (uiV == (h - 1))
        {

        }
        for (BYTE qm = 0; qm < vectorPower; ++qm)
        {
            wholeCircuit.AppendGate(cc, qm);
        }

        if (0 != uiV)
        {
            //memcpy(vec.stateVec.real, realpart, sizeof(Real) * veclen);
            //memcpy(vec.stateVec.imag, imagpart, sizeof(Real) * veclen);
            syncQuESTEnv(evn);
            copyStateToGPU(vec);
        }

        TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
        SIZE_T opssize = ops.Num();
        for (SIZE_T i = 0; i < opssize; ++i)
        {
            QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        }
        //syncQuESTEnv(evn);
        calcProbOfAllOutcomes(res, vec, qubitsToSee.GetData(), byMeasureQubits);

        //copyStateFromGPU(vec);
        //UINT uiHitRes = 0;
        memset(rescount, 0, sizeof(INT) * uiPossibleRes);
        for (UINT hit = 0; hit < uiRepeat; ++hit)
        {
            Real fHit = RandomF();
            Real fRight = res[0];
            for (UINT ires = 0; ires < uiPossibleRes; ++ires)
            {
                if (fHit < fRight)
                {
                    ++rescount[ires];
                    break;
                }
                
                if (ires == (uiPossibleRes - 2))
                {
                    ++rescount[uiPossibleRes - 1];
                    break;
                }
                fRight = fRight + res[ires + 1];
            }
        }
        CCString sRes;
        const Real fDenorm = static_cast<Real>(uiRepeat);
        for (UINT ires = 0; ires < uiPossibleRes; ++ires)
        {
            Real fScore = rescount[ires] / fDenorm;
            sRes = sRes + appToString(fScore) + _T(", ");
            scores[uiV * uiPossibleRes + ires] = fScore;
        }
        appGeneral(_T("finish %d / %d , hit = %s ...\n"), uiV, h, sRes.c_str());
    }

    appSafeFree(pHostY);
    appSafeFree(pHostZ);
    //appSafeFree(realpart);
    //appSafeFree(imagpart);

    //final export host center
    SaveCSVAR(scores, uiPossibleRes, h, sScore);
    appSafeFree(scores);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
}

void QLQuantumKmeans::KNNAE(const CCString& sTrainingPoints, const CCString& sTestPoints, const CCString& sScore, BYTE byMeasureQubits, UINT uiRepeat)
{
    //UINT uiBlock = 0;
    //UINT uiThread = 0;
    UINT w, h;
    TArray<QLComplex> trainpoints = ReadCSVA(sTrainingPoints, w, h);
    QLGate aegate = AmplitudeEncodeVectors(trainpoints.GetData(), MostSignificantPowerTwo(h), MostSignificantPowerTwo(w), FALSE);
    BYTE allBytes = static_cast<BYTE>(aegate.m_lstQubits.Num());
    TArray<QLComplex> orignalPointsArray = ReadCSVA(sTestPoints, w, h);
    BYTE vectorPower = static_cast<BYTE>(MostSignificantPowerTwo(w));

    QLComplex* pDeviceVectors = NULL;
    Real* pDeviceAbs = NULL;
    Real* pDevicePhase = NULL;
    Real* pDeviceY = NULL;
    Real* pDeviceZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceVectors, sizeof(QLComplex) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDeviceAbs, sizeof(QLComplex) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDevicePhase, sizeof(QLComplex) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDeviceY, sizeof(QLComplex) * w * h));
    checkCudaErrors(cudaMalloc((void**)&pDeviceZ, sizeof(QLComplex) * w * h));
    checkCudaErrors(cudaMemcpy(pDeviceVectors, orignalPointsArray.GetData(), sizeof(QLComplex) * w * h, cudaMemcpyHostToDevice));
    CalculateDegreesForEach(pDeviceVectors, pDeviceAbs, pDevicePhase, h, vectorPower, pDeviceY, pDeviceZ);
    Real* pHostY = reinterpret_cast<Real*>(malloc(sizeof(Real) * w * h));
    Real* pHostZ = reinterpret_cast<Real*>(malloc(sizeof(Real) * w * h));
    checkCudaErrors(cudaMemcpy(pHostY, pDeviceY, sizeof(Real) * w * h, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostZ, pDeviceZ, sizeof(Real) * w * h, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(pDeviceVectors));
    checkCudaErrors(cudaFree(pDeviceAbs));
    checkCudaErrors(cudaFree(pDevicePhase));
    checkCudaErrors(cudaFree(pDeviceY));
    checkCudaErrors(cudaFree(pDeviceZ));

    UINT uiPossibleRes = 1 << byMeasureQubits;
    Real* scores = reinterpret_cast<Real*>(malloc(sizeof(Real) * h * uiPossibleRes));

    UINT veclen = 1UL << static_cast<UINT>(allBytes);
    QLGate cc(EBasicOperation::EBO_CC, 0);

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(allBytes, evn);

    
    QLGate wholeCircuit;
    wholeCircuit.AddQubits(allBytes);
    wholeCircuit.AppendGate(aegate, aegate.m_lstQubits);

    syncQuESTEnv(evn);
    memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
    memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
    vec.stateVec.real[0] = F(1.0);
    vec.stateVec.imag[0] = F(0.0);
    copyStateToGPU(vec);
    TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
    SIZE_T opssize = ops.Num();
    for (SIZE_T i = 0; i < opssize; ++i)
    {
        QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
    }
    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    Real* realpart = reinterpret_cast<Real*>(malloc(sizeof(Real) * veclen));
    Real* imagpart = reinterpret_cast<Real*>(malloc(sizeof(Real) * veclen));
    memcpy(realpart, vec.stateVec.real, sizeof(Real) * veclen);
    memcpy(imagpart, vec.stateVec.imag, sizeof(Real) * veclen);

    Real* res = reinterpret_cast<Real*>(appAlloca(sizeof(Real) * uiPossibleRes));
    INT* rescount = reinterpret_cast<INT*>(appAlloca(sizeof(INT) * uiPossibleRes));
    TArray<INT> qubitsToSee;
    for (BYTE bym = 0; bym < byMeasureQubits; ++bym)
    {
        qubitsToSee.AddItem(static_cast<INT>(allBytes - 1 - bym));
    }

    for (UINT uiV = 0; uiV < h; ++uiV)
    {
        QLGate newcircuit;
        newcircuit.AddQubits(allBytes);
        //QLGate vectorAE = AmplitudeEncodeOneVector(orignalPointsArray.GetData() + w * uiV, vectorPower, FALSE);
        QLGate vectorAE = ExchangeToYZGate(vectorPower, pHostY + w * uiV, pHostZ + w * uiV, FALSE);
        vectorAE.Dagger();
        newcircuit.AppendGate(vectorAE, vectorAE.m_lstQubits);
        for (BYTE qm = 0; qm < vectorPower; ++qm)
        {
            newcircuit.AppendGate(cc, qm);
        }
        
        if (0 != uiV)
        {
            memcpy(vec.stateVec.real, realpart, sizeof(Real) * veclen);
            memcpy(vec.stateVec.imag, imagpart, sizeof(Real) * veclen);
            syncQuESTEnv(evn);
            copyStateToGPU(vec);
        }

        TArray<SBasicOperation> ops2 = newcircuit.GetOperation(newcircuit.m_lstQubits);
        SIZE_T opssize2 = ops2.Num();
        for (SIZE_T i = 0; i < opssize2; ++i)
        {
            QLGate::PerformBasicOperation(vec, ops2[static_cast<INT>(i)]);
        }

        calcProbOfAllOutcomes(res, vec, qubitsToSee.GetData(), byMeasureQubits);

        //copyStateFromGPU(vec);
        //UINT uiHitRes = 0;
        memset(rescount, 0, sizeof(INT) * uiPossibleRes);
        for (UINT hit = 0; hit < uiRepeat; ++hit)
        {
            Real fHit = RandomF();
            Real fRight = res[0];
            for (UINT ires = 0; ires < uiPossibleRes; ++ires)
            {
                if (fHit < fRight)
                {
                    ++rescount[ires];
                    break;
                }

                if (ires == (uiPossibleRes - 2))
                {
                    ++rescount[uiPossibleRes - 1];
                    break;
                }
                fRight = fRight + res[ires + 1];
            }
        }
        CCString sRes;
        const Real fDenorm = static_cast<Real>(uiRepeat);
        for (UINT ires = 0; ires < uiPossibleRes; ++ires)
        {
            Real fScore = rescount[ires] / fDenorm;
            sRes = sRes + appToString(fScore) + _T(", ");
            scores[uiV * uiPossibleRes + ires] = fScore;
        }
        appGeneral(_T("finish %d / %d , hit = %s ...\n"), uiV, h, sRes.c_str());
    }

    //final export host center
    SaveCSVAR(scores, uiPossibleRes, h, sScore);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
}

void QLQuantumKmeans::QAnomaly2D(const CCString& sReferenceCSV, const CCString& sPointCSV, const CCString& sBuildRate,
    Real minX, Real maxX, Real minY, Real maxY)
{
    UINT uiBlock = 0;
    UINT uiThread = 0;
    UINT w1, h1;
    UINT w2, h2;
    TArray<Real> orignalPointsArrayTraining = ReadCSVAR(sReferenceCSV, w1, h1);
    TArray<Real> orignalPointsArrayTesting = ReadCSVAR(sPointCSV, w2, h2);

    Real* pTrainingX = NULL;
    Real* pTrainingY = NULL;
    Real* pTestingX = NULL;
    Real* pTestingY = NULL;
    Real* pWorkingSpaceTraining = NULL;
    Real* pWorkingSpaceTesting = NULL;

    checkCudaErrors(cudaMalloc((void**)&pWorkingSpaceTraining, sizeof(Real) * h1 * 2));
    checkCudaErrors(cudaMalloc((void**)&pWorkingSpaceTesting, sizeof(Real) * h2 * 2));
    checkCudaErrors(cudaMalloc((void**)&pTrainingX, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pTrainingY, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pTestingX, sizeof(Real) * h2));
    checkCudaErrors(cudaMalloc((void**)&pTestingY, sizeof(Real) * h2));

    checkCudaErrors(cudaMemcpy(pWorkingSpaceTraining, orignalPointsArrayTraining.GetData(), sizeof(Real) * h1 * 2, cudaMemcpyHostToDevice));
    __DECOMPOSE(h1, uiBlock, uiThread);
    _kernelKMeans2DReadFileXY << <uiBlock, uiThread >> > (pWorkingSpaceTraining, pTrainingX, pTrainingY, h1);
    checkCudaErrors(cudaMemcpy(pWorkingSpaceTesting, orignalPointsArrayTesting.GetData(), sizeof(Real) * h2 * 2, cudaMemcpyHostToDevice));
    __DECOMPOSE(h2, uiBlock, uiThread);
    _kernelKMeans2DReadFileXY << <uiBlock, uiThread >> > (pWorkingSpaceTesting, pTestingX, pTestingY, h2);

    Real* pAngleYTraining = NULL;
    Real* pAngleZTraining = NULL;
    Real* pAngleYTesting = NULL;
    Real* pAngleZTesting = NULL;
    checkCudaErrors(cudaMalloc((void**)&pAngleYTraining, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pAngleZTraining, sizeof(Real) * h1));
    checkCudaErrors(cudaMalloc((void**)&pAngleYTesting, sizeof(Real) * h2));
    checkCudaErrors(cudaMalloc((void**)&pAngleZTesting, sizeof(Real) * h2));

    __DECOMPOSE(h1, uiBlock, uiThread);
    _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pTrainingX, pTrainingY,
        minX, maxX, minY, maxY, pAngleYTraining, pAngleZTraining, h1);
    __DECOMPOSE(h2, uiBlock, uiThread);
    _kernelKMeans2DPointToAbsPhase << <uiBlock, uiThread >> > (pTestingX, pTestingY,
        minX, maxX, minY, maxY, pAngleYTesting, pAngleZTesting, h2);

    Real* pHostAngleYTraining = reinterpret_cast<Real*>(malloc(sizeof(Real) * h1));
    Real* pHostAngleZTraining = reinterpret_cast<Real*>(malloc(sizeof(Real) * h1));
    Real* pHostAngleYTesting = reinterpret_cast<Real*>(malloc(sizeof(Real) * h2));
    Real* pHostAngleZTesting = reinterpret_cast<Real*>(malloc(sizeof(Real) * h2));

    checkCudaErrors(cudaMemcpy(pHostAngleYTraining, pAngleYTraining, sizeof(Real) * h1, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleZTraining, pAngleZTraining, sizeof(Real) * h1, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleYTesting, pAngleYTesting, sizeof(Real) * h2, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostAngleZTesting, pAngleZTesting, sizeof(Real) * h2, cudaMemcpyDeviceToHost));

    TArray<Real> fBuildRates;
    UINT kPower = MostSignificantPowerTwo(h1);
    UINT veclen = 1UL << (kPower + 1);
    QLGate cc(EBasicOperation::EBO_CC, 0);
    QLGate had(EBasicOperation::EBO_H);

    TArray<BYTE> frybits;
    for (BYTE i = 0; i < kPower + 1; ++i)
    {
        frybits.AddItem(static_cast<BYTE>(kPower - i));
    }

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(kPower + 1, evn);

    for (UINT uiV = 0; uiV < h2; ++uiV)
    {
        QLGate wholeCircuit;
        wholeCircuit.AddQubits(static_cast<BYTE>(kPower + 1));
        for (BYTE hadIdx = 0; hadIdx < kPower; ++hadIdx)
        {
            wholeCircuit.AppendGate(had, hadIdx + 1);
        }
        QLGate amplitudeVi = FRyz(pHostAngleYTraining, pHostAngleZTraining, kPower + 1);
        wholeCircuit.AppendGate(amplitudeVi, frybits);

        //ry(theta) rz(-phi)
        //and dagger, which is rz(phi), ry(-theta)
        QLGate rz(EBasicOperation::EBO_RZ, -pHostAngleZTesting[uiV]);
        QLGate ry(EBasicOperation::EBO_RY, pHostAngleYTesting[uiV]);
        QLGate amp2;
        amp2.AddQubits(1);
        amp2.AppendGate(ry, 0);
        amp2.AppendGate(rz, 0);
        amp2.Dagger();

        wholeCircuit.AppendGate(amp2, 0);
        wholeCircuit.AppendGate(cc, 0);

        syncQuESTEnv(evn);
        memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
        memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
        vec.stateVec.real[0] = F(1.0);
        vec.stateVec.imag[0] = F(0.0);
        copyStateToGPU(vec);

        TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
        SIZE_T opssize = ops.Num();
        Real fBuildRate = F(1.0);
        for (SIZE_T i = 0; i < opssize; ++i)
        {
            fBuildRate = fBuildRate * QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        }
        fBuildRates.AddItem(fBuildRate);
        syncQuESTEnv(evn);

        copyStateFromGPU(vec);

#pragma region test state-vector

        //if (0 == uiV)
        //{
        //    TArray<BYTE> getout;
        //    getout.AddItem(0);
        //    for (BYTE toshowbit = 0; toshowbit < kPower; ++toshowbit)
        //    {
        //        getout.AddItem(2);
        //    }
        //    QLMatrix show = ShowStateVectorDetail(vec, getout);
        //    show.Mul(_make_cuComplex(F(23.34282757393504), F(0.0)));
        //    show.Print();
        //}

#pragma endregion

        appGeneral(_T("br = %2.12f, finish %d / %d ...\n"), fBuildRate, uiV, h2);
    }

    SaveCSVAR(fBuildRates.GetData(), 1, h2, sBuildRate);

    //final export host center
    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    checkCudaErrors(cudaFree(pAngleYTraining));
    checkCudaErrors(cudaFree(pAngleZTraining));
    checkCudaErrors(cudaFree(pAngleYTesting));
    checkCudaErrors(cudaFree(pAngleZTesting));

    appSafeFree(pHostAngleYTraining);
    appSafeFree(pHostAngleZTraining);
    appSafeFree(pHostAngleYTesting);
    appSafeFree(pHostAngleZTesting);

    checkCudaErrors(cudaFree(pTrainingX));
    checkCudaErrors(cudaFree(pTrainingY));
    checkCudaErrors(cudaFree(pTestingX));
    checkCudaErrors(cudaFree(pTestingY));
    checkCudaErrors(cudaFree(pWorkingSpaceTraining));
    checkCudaErrors(cudaFree(pWorkingSpaceTesting));
}

void QLQuantumKmeans::QAnomaly3D(const CCString& sTrainingPoints, const CCString& sTestPoints, const CCString& sSaveScore)
{
    UINT w1, h1;
    UINT w2, h2;
    TArray<Real> orignalPointsArrayTraining = ReadCSVAR(sTrainingPoints, w1, h1);
    TArray<Real> orignalPointsArrayTesting = ReadCSVAR(sTestPoints, w2, h2);

    Real* pTrainingPoints = NULL;
    Real* pTrainingY = NULL;
    Real* pTestingPoints = NULL;
    Real* pTestingY = NULL;
    Real* pHostTrainingY = reinterpret_cast<Real*>(malloc(sizeof(Real) * h1 * 4));
    Real* pHostTestY = reinterpret_cast<Real*>(malloc(sizeof(Real) * h2 * 4));

    checkCudaErrors(cudaMalloc((void**)&pTrainingPoints, sizeof(Real) * h1 * 4));
    checkCudaErrors(cudaMalloc((void**)&pTrainingY, sizeof(Real) * h1 * 4));
    checkCudaErrors(cudaMalloc((void**)&pTestingPoints, sizeof(Real) * h2 * 4));
    checkCudaErrors(cudaMalloc((void**)&pTestingY, sizeof(Real) * h2 * 4));

    checkCudaErrors(cudaMemcpy(pTrainingPoints, orignalPointsArrayTraining.GetData(), sizeof(Real) * h1 * 4, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(pTestingPoints, orignalPointsArrayTesting.GetData(), sizeof(Real) * h2 * 4, cudaMemcpyHostToDevice));

    CalculateDegreesReal(pTrainingPoints, h1, 2, pTrainingY);
    CalculateDegreesRealForEach(pTestingPoints, h2, 2, pTestingY);
    checkCudaErrors(cudaMemcpy(pHostTrainingY, pTrainingY, sizeof(Real) * h1 * 4, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(pHostTestY, pTestingY, sizeof(Real) * h2 * 4, cudaMemcpyDeviceToHost));

    UINT kPower = MostSignificantPowerTwo(h1);
    //2 for vector power, 2 for amplitude of CSwap, 1 for CSWap measure
    UINT veclen = 1UL << (kPower + 2);
    
    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(kPower + 2, evn);

    //QLGate cswap = CreateControlledSwap(1);
    //QLGate h(EBasicOperation::EBO_H);
    QLGate cc(EBasicOperation::EBO_CC, 0);
    TArray<Real> fOverlaps;

    //Iterations start
    for (UINT uiV = 0; uiV < h2; ++uiV)
    {
        QLGate wholeCircuit;
        wholeCircuit.AddQubits(static_cast<BYTE>(kPower + 2));
        QLGate amplitudeVi = ExchangeToYGate(kPower, 2, pHostTrainingY);
        wholeCircuit.AppendGate(amplitudeVi, amplitudeVi.m_lstQubits);

        //ry(theta) rz(-phi)
        //and dagger, which is rz(phi), ry(-theta)
        QLGate amplitudeU = ExchangeToYGate(2, pHostTestY + uiV * 4);
        amplitudeU.Dagger();

        wholeCircuit.AppendGate(amplitudeU, 0, 1);

        wholeCircuit.AppendGate(cc, 0);
        wholeCircuit.AppendGate(cc, 1);

        //wholeCircuit.AppendGate(h, kPower + 4);
        //wholeCircuit.AppendGate(cswap, kPower + 4, 0, kPower + 2);
        //wholeCircuit.AppendGate(cswap, kPower + 4, 1, kPower + 3);
        //wholeCircuit.AppendGate(h, kPower + 4);

        syncQuESTEnv(evn);
        memset(vec.stateVec.real, 0, sizeof(Real) * veclen);
        memset(vec.stateVec.imag, 0, sizeof(Real) * veclen);
        vec.stateVec.real[0] = F(1.0);
        vec.stateVec.imag[0] = F(0.0);
        copyStateToGPU(vec);

        TArray<SBasicOperation> ops = wholeCircuit.GetOperation(wholeCircuit.m_lstQubits);
        SIZE_T opssize = ops.Num();
        Real fProb = F(1.0);
        for (SIZE_T i = 0; i < opssize; ++i)
        {
            fProb = fProb * QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
        }

        syncQuESTEnv(evn);
        
        //copyStateFromGPU(vec);

#pragma region test state-vector

        //if (0 == uiV)
        //{
        //    TArray<BYTE> getout;
        //    getout.AddItem(0);
        //    for (BYTE toshowbit = 0; toshowbit < kPower; ++toshowbit)
        //    {
        //        getout.AddItem(2);
        //    }
        //    QLMatrix show = ShowStateVectorDetail(vec, getout);
        //    show.Mul(_make_cuComplex(F(23.34282757393504), F(0.0)));
        //    show.Print();
        //}

#pragma endregion

        //INT tomeasure[1] = { kPower + 4 };
        //Real measureprob[2];
        //calcProbOfAllOutcomes(measureprob, vec, tomeasure, 1);
        //Real fZeroProb = F(2.0) * measureprob[0] - 1;
        fOverlaps.AddItem(fProb);

        appGeneral(_T("finish br:%f,  %d / %d ...\n"), fProb, uiV, h2);
    }

    //final export host center
    SaveCSVAR(fOverlaps.GetData(), 1, h2, sSaveScore);

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    appSafeFree(pHostTrainingY);
    appSafeFree(pHostTestY);

    checkCudaErrors(cudaFree(pTrainingY));
    checkCudaErrors(cudaFree(pTestingY));
    checkCudaErrors(cudaFree(pTrainingPoints));
    checkCudaErrors(cudaFree(pTestingPoints));
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================