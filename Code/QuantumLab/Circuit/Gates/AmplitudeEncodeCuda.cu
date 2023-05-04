//=============================================================================
// FILENAME : AmplitudeEncodeCuda.cu
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/03/2023 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* uiMax = len(v) * #v
*/
__global__ void _QL_LAUNCH_BOUND
_kernelAE_VtoAbsAndPhase(const QLComplex * __restrict__ v, Real* absBuffer, Real* phaseBuffer, UINT uiMax)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        absBuffer[idx] = _cuCabsf(v[idx]);
        phaseBuffer[idx] = __cuCargf(v[idx]);
    }
}

/**
* when level = n
* calculate 1 << (lengthPower - 1 - n) angles
* if n != (lengthPower - 1), combine the abs
* 
* level = 0:
* c s c s c s c s c s c s c s c s
* cs x cs x cs x cs x cs x cs x cs x cs x
* to put dddddddd xxxxxxxx
* 
* level = 1:
* c x s x c x s x c x s c x s x
* cs xxx cs xxx cs xxx cs xxx
* to put xxxxxxxx dddd xxxx
* 
* level = 2:
* c xxx s xxx c xxx s xxx
*/
__global__ void _QL_LAUNCH_BOUND
_kernelAE_Rotations(Real* YLst, Real* ZLst, Real* absBuffer, Real* phaseBuffer, UINT lengthPower, UINT vectorCount, UINT uiLevel)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    const UINT uiVLen = 1U << lengthPower;

    if (lengthPower == uiLevel && NULL != ZLst && NULL != phaseBuffer)
    {
        //only record the total phase
        ZLst[vectorCount * (uiVLen - 1) + idx] = phaseBuffer[uiVLen * idx];
        return;
    }

    const UINT numberOfDegrees = 1U << (lengthPower - uiLevel - 1);
    const UINT uiIdOfV = idx / numberOfDegrees;

    if (uiIdOfV < vectorCount)
    {
        const UINT uiIdInV = idx % numberOfDegrees;

        const UINT uiLeft = uiIdOfV * uiVLen + uiIdInV * (1U << (uiLevel + 1));
        const UINT uiRight = uiLeft + (1U << uiLevel);
        UINT putStart = (1U << lengthPower) - (1U << (lengthPower - uiLevel));
        putStart = vectorCount * putStart + uiIdOfV * numberOfDegrees + uiIdInV;

        const Real cs = absBuffer[uiLeft];
        const Real sn = absBuffer[uiRight];
        YLst[putStart] = _atan2(sn, cs) * F(2.0);

        if (1 != numberOfDegrees)
        {
            absBuffer[uiLeft] = absBuffer[uiLeft] * absBuffer[uiLeft] + absBuffer[uiRight] * absBuffer[uiRight];
            if (absBuffer[uiLeft] > _QL_FLT_EPSILON)
            {
                absBuffer[uiLeft] = _sqrt(absBuffer[uiLeft]);
            }
            else
            {
                absBuffer[uiLeft] = F(0.0);
            }
        }

        if (NULL != ZLst && NULL != phaseBuffer)
        {
            ZLst[putStart] = phaseBuffer[uiLeft] - phaseBuffer[uiRight];
            phaseBuffer[uiLeft] = F(0.5) * (phaseBuffer[uiLeft] + phaseBuffer[uiRight]);
        }
    }
}

/**
* (factor, r1 + r2 I, r3 + r3 I, ...)
*/
__global__ void _QL_LAUNCH_BOUND
_kernelAE_Normalize(const Real* __restrict__ deviceData, QLComplex* deviceNormalized, UINT realVectorLength, UINT complexVectorLength, UINT uiMax, Real factor)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (idx < uiMax)
    {
        Real sum = factor * factor;
        const UINT cpxIdStart = idx * complexVectorLength;
        const UINT realIdStart = idx * realVectorLength;
        deviceNormalized[cpxIdStart] = _make_cuComplex(factor, F(0.0));
        for (UINT i = 0; i < realVectorLength / 2; ++i)
        {
            const UINT idStart = realIdStart + 2 * i;
            const Real re = deviceData[idStart];
            const Real im = deviceData[idStart + 1];
            deviceNormalized[cpxIdStart + i + 1] = _make_cuComplex(re, im);
            sum = sum + re * re + im * im;
        }
        sum = _sqrt(sum);
        for (UINT i = 0; i < complexVectorLength; ++i)
        {
            deviceNormalized[cpxIdStart + i] = cuCdivf_cr(deviceNormalized[cpxIdStart + i], sum);
        }
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelAE_CompareDistance(const Real* __restrict__ deviceData, const QLComplex* __restrict__ deviceNormalized, UINT realVectorLength, UINT complexVectorLength, Real* result, UINT uiMax)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (idx < uiMax)
    {
        const UINT idStartAll = (idx << 1U);
        UINT idStart = idStartAll * realVectorLength;
        Real distance = F(0.0);
        for (UINT i = 0; i < realVectorLength; ++i)
        {
            distance += (deviceData[idStart + i] - deviceData[idStart + i + realVectorLength]) * (deviceData[idStart + i] - deviceData[idStart + i + realVectorLength]);
        }
        result[idStartAll] = _sqrt(distance);

        idStart = idStartAll * complexVectorLength;
        QLComplex sum = _zeroc;
        for (UINT i = 0; i < complexVectorLength; ++i)
        {
            sum = _cuCaddf(sum, _cuCmulf(_cuConjf(deviceNormalized[idStart + i]), deviceNormalized[idStart + i + complexVectorLength]));
        }
        result[idStartAll + 1] = _cuCabsf(sum);
    }
}

#pragma endregion

void QLAPI NormalizeRealToComplex(const Real* deviceData, QLComplex* deviceNormalized, UINT realVectorLength, UINT complexVectorLength, UINT uiVectorCount, Real factor)
{
    UINT iBlock1 = uiVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(uiVectorCount, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = uiVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(uiVectorCount, iBlock1) : uiVectorCount;
    _kernelAE_Normalize << <iBlock1, iThread1 >> > (deviceData, deviceNormalized, realVectorLength, complexVectorLength, uiVectorCount, factor);
}

void QLAPI CalculateDegrees(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ)
{
    UINT vectorLength = 1U << vectorPower;
    UINT uiLen = vectorCount * vectorLength;
    UINT iBlock1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, iBlock1) : uiLen;
    _kernelAE_VtoAbsAndPhase << <iBlock1, iThread1 >> > (deviceV, absBuffer, phaseBuffer, uiLen);

    for (UINT i = 0; i <= vectorPower; ++i)
    {
        if (i != vectorPower)
        {
            UINT uiLen2 = vectorCount * (1U << (vectorPower - i - 1));
            UINT iBlock2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, _QL_LAUNCH_MAX_THREAD) : 1;
            UINT iThread2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, iBlock2) : uiLen2;
            _kernelAE_Rotations << <iBlock2, iThread2 >> > (deviceY, deviceZ, absBuffer, phaseBuffer, vectorPower, vectorCount, i);
        }
        else
        {
            UINT uiLen2 = vectorCount;
            UINT iBlock2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, _QL_LAUNCH_MAX_THREAD) : 1;
            UINT iThread2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, iBlock2) : uiLen2;
            _kernelAE_Rotations << <iBlock2, iThread2 >> > (deviceY, deviceZ, absBuffer, phaseBuffer, vectorPower, vectorCount, i);
        }
    }
}

void QLAPI CalculateDegrees(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ, Real* hostY, Real* hostZ)
{
    CalculateDegrees(deviceV, absBuffer, phaseBuffer, vectorCount, vectorPower, deviceY, deviceZ);
    UINT uiLen = vectorCount * (1U << vectorPower);

    checkCudaErrors(cudaMemcpy(hostY, deviceY, sizeof(Real) * uiLen, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(hostZ, deviceZ, sizeof(Real) * uiLen, cudaMemcpyDeviceToHost));
}

QLGate QLAPI ExchangeToYZGate(UINT vectorPower, Real* hostY, Real* hostZ, UBOOL bHasPhase)
{
    /**
    * the degrees is 
    * 00000000 1111 22 3 4
    
    * yz(00000000) at 0123
    * yz(1111) at 012
    * yz(22) at 01
    * we need yz(3) at 0
    * 
    * and then global phase at 0
    */

    QLGate retGate;
    retGate.AddQubits(static_cast<BYTE>(vectorPower));
    retGate.m_sName = _T("AE");

    UINT vectorLength = (1U << vectorPower);
    UINT idx = vectorLength - 2;
    for (UINT i = 0; i < vectorPower; ++i)
    {
        TArray<BYTE> bits;
        for (UINT j = 0; j < i + 1; ++j)
        {
            //we need
            //3, 32, 321, 3210 
            bits.AddItem(static_cast<BYTE>(vectorPower - 1 - j));
        }
        if (1 == bits.Num())
        {
            QLGate ry(EBasicOperation::EBO_RY, hostY[idx]);
            QLGate rz(EBasicOperation::EBO_RZ, -hostZ[idx]);
            retGate.AppendGate(ry, bits);
            retGate.AppendGate(rz, bits);
            idx = idx - 2;
        }
        else
        {
            QLGate fryz = FRyz(hostY + idx, hostZ + idx, static_cast<UINT>(bits.Num()));
            idx = idx - (1U << (i + 1));
            retGate.AppendGate(fryz, bits);
        }
    }

    if (bHasPhase)
    {
        QLGate ph(EBasicOperation::EBO_Phase, hostZ[vectorLength - 1]);
        TArray <BYTE> bit0;
        bit0.AddItem(static_cast<BYTE>(vectorPower - 1));
        retGate.AppendGate(ph, bit0);
    }
    return retGate;
}

/**
* Apply the phases to the vectors
*/
static void inline hostAE_CalculatePhases(const Real* phases, Real* degreeZ, UINT lengthPower)
{
    UINT phaseCount = 1U << lengthPower;
    Real a0 = F(0.0);
    for (INT i = 0; i < phaseCount; ++i)
    {
        a0 += phases[i];
    }
    a0 = a0 / phaseCount;
    degreeZ[phaseCount - 1] = a0;

    UINT idxStart = phaseCount - 2;
    for (INT i = 0; i < lengthPower; ++i)
    {
        UINT degreeCount = 1U << i;
        Real factor = F(1.0) / (1U << (lengthPower - i));
        UINT sumLength = 1U << (lengthPower - i - 1);
        for (UINT j = 0; j < degreeCount; ++j)
        {
            Real plus = F(0.0);
            for (UINT k = 2 * j * sumLength; k < (2 * j + 1) * sumLength; ++k)
            {
                plus = plus + phases[k];
            }
            Real minus = F(0.0);
            for (UINT k = (2 * j + 1) * sumLength; k < (2 * j + 2) * sumLength; ++k)
            {
                minus = minus + phases[k];
            }
            degreeZ[idxStart + j] = factor * (plus - minus) * F(-2.0);
        }
        idxStart = idxStart - (degreeCount << 1U);
    }
}

QLGate QLAPI ExchangeToYZGate(UINT vectorCountPower, UINT vectorPower, Real* hostY, Real* hostZ, UBOOL bHasPhase)
{
    QLGate retGate;
    retGate.AddQubits(static_cast<BYTE>(vectorCountPower + vectorPower));
    retGate.m_sName = _T("QRAM");
    QLGate h(EBasicOperation::EBO_H);
    for (UINT i = 0; i < vectorCountPower; ++i)
    {
        TArray<BYTE> hadmard;
        hadmard.AddItem(static_cast<BYTE>(i + vectorPower));
        retGate.AppendGate(h, hadmard);
        retGate.m_lstAdditionalQubitsData.AddItem(i + vectorPower);
    }

    UINT vectorCount = (1U << vectorCountPower);
    UINT vectorLength = (1U << vectorPower);

    for (UINT i = 0; i < vectorPower; ++i)
    {
        TArray<BYTE> bits;
        for (UINT j = 0; j < vectorCountPower; ++j)
        {
            //For example, if we have 16 vectors (length-8), we need
            //5432, 54321, 543210, ...
            //here we need 3, 43, and 543, where vectorCountPower = 3
            bits.AddItem(static_cast<BYTE>(vectorCountPower + vectorPower - j - 1));
        }
        for (UINT j = 0; j <= i; ++j)
        {
            //put 2,1,0
            bits.AddItem(static_cast<BYTE>(vectorPower - 1 - j));
        }

        UINT startInV = vectorCount * (vectorLength - (1U << (i + 1)));

        QLGate fryz = FRyz(hostY + startInV, hostZ + startInV, static_cast<UINT>(bits.Num()));

        retGate.AppendGate(fryz, bits);
    }

    if (bHasPhase)
    {
        UINT uiIdxStart = vectorCount * (vectorLength - 1);
        hostAE_CalculatePhases(hostZ + uiIdxStart, hostY + uiIdxStart, vectorCountPower);
        UINT idx = vectorCount - 2;
        for (UINT i = 0; i < vectorCountPower; ++i)
        {
            TArray<BYTE> bits;
            for (UINT j = 0; j < i + 1; ++j)
            {
                //we need
                //3, 32, 321, 3210 
                bits.AddItem(static_cast<BYTE>(vectorCountPower - 1 - j + vectorPower));
            }
            if (1 == bits.Num())
            {
                QLGate rz(EBasicOperation::EBO_RZ, hostY[uiIdxStart + idx]);
                retGate.AppendGate(rz, bits);
                idx = idx - 2;
            }
            else
            {
                QLGate frz = FRz(hostY + uiIdxStart + idx, static_cast<UINT>(bits.Num()));
                idx = idx - (1U << (i + 1));
                retGate.AppendGate(frz, bits);
            }
        }

        QLGate ph(EBasicOperation::EBO_Phase, hostY[vectorCount + uiIdxStart - 1]);
        TArray <BYTE> bit0;
        bit0.AddItem(static_cast<BYTE>(vectorCountPower + vectorPower - 1));
        retGate.AppendGate(ph, bit0);
    }

    return retGate;
}

QLGate QLAPI AmplitudeEncodeOneVector(const QLComplex* hostv, UINT vectorPower, UBOOL bHasPhase)
{
    UINT vLength = 1U << vectorPower;

    QLComplex* v = NULL;
    Real* absBuffer = NULL;
    Real* phaseBuffer = NULL;
    Real* YBuffer = NULL;
    Real* ZBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&v, sizeof(QLComplex) * vLength));
    checkCudaErrors(cudaMalloc((void**)&absBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMalloc((void**)&phaseBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMalloc((void**)&YBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMalloc((void**)&ZBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMemcpy(v, hostv, sizeof(QLComplex) * vLength, cudaMemcpyHostToDevice));

    Real* YHostBuffer = (Real*)malloc(sizeof(Real) * vLength);
    Real* ZHostBuffer = (Real*)malloc(sizeof(Real) * vLength);

    CalculateDegrees(v, absBuffer, phaseBuffer, 1, vectorPower, YBuffer, ZBuffer, YHostBuffer, ZHostBuffer);

    checkCudaErrors(cudaFree(v));
    checkCudaErrors(cudaFree(absBuffer));
    checkCudaErrors(cudaFree(phaseBuffer));
    checkCudaErrors(cudaFree(YBuffer));
    checkCudaErrors(cudaFree(ZBuffer));

    QLGate ret = ExchangeToYZGate(vectorPower, YHostBuffer, ZHostBuffer, bHasPhase);

    free(YHostBuffer);
    free(ZHostBuffer);

    return ret;
}

QLGate QLAPI AmplitudeEncodeVectors(const QLComplex* hostv, UINT vectorCountPower, UINT vectorPower, UBOOL bHasPhase)
{
    UINT vLength = 1U << vectorPower;
    UINT vCount = 1U << vectorCountPower;
    UINT totalLength = vLength * vCount;

    QLComplex* v = NULL;
    Real* absBuffer = NULL;
    Real* phaseBuffer = NULL;
    Real* YBuffer = NULL;
    Real* ZBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&v, sizeof(QLComplex) * totalLength));
    checkCudaErrors(cudaMalloc((void**)&absBuffer, sizeof(Real) * totalLength));
    checkCudaErrors(cudaMalloc((void**)&phaseBuffer, sizeof(Real) * totalLength));
    checkCudaErrors(cudaMalloc((void**)&YBuffer, sizeof(Real) * totalLength));
    checkCudaErrors(cudaMalloc((void**)&ZBuffer, sizeof(Real) * totalLength));
    checkCudaErrors(cudaMemcpy(v, hostv, sizeof(QLComplex) * totalLength, cudaMemcpyHostToDevice));

    Real* YHostBuffer = (Real*)malloc(sizeof(Real) * totalLength);
    Real* ZHostBuffer = (Real*)malloc(sizeof(Real) * totalLength);

    CalculateDegrees(v, absBuffer, phaseBuffer, vCount, vectorPower, YBuffer, ZBuffer, YHostBuffer, ZHostBuffer);

    checkCudaErrors(cudaFree(v));
    checkCudaErrors(cudaFree(absBuffer));
    checkCudaErrors(cudaFree(phaseBuffer));
    checkCudaErrors(cudaFree(YBuffer));
    checkCudaErrors(cudaFree(ZBuffer));

    QLGate ret = ExchangeToYZGate(vectorCountPower, vectorPower, YHostBuffer, ZHostBuffer, bHasPhase);

    free(YHostBuffer);
    free(ZHostBuffer);

    return ret;
}



/**
* This function receive host data which has k+1 vectors, with length = uiVectorLength
* 
* uiVectorLength = 6, 14, 30, 62, 126j
* factors        = 1.5, 2.5, 4.0, 5.0, 6.0
* 
* 
*/
void QLAPI TestProbalityToBuildAmplitude(UINT uiVectorLength, UINT uiKPower, Real fFactor, Real& probToBuild, Real& probToObtainLargest, UBOOL bDebugPrint)
{
    UINT uiK = 1U << uiKPower;

    QLMatrix randomR(uiVectorLength, uiK + 1);
    randomR.RandomOneReal();
    TArray<Real> realList = randomR.ToVectorRe();

    UINT uiComplexLength = 1 + uiVectorLength / 2;

    Real* tempRealBuffer = NULL;
    QLComplex* tempComplexBuffer = NULL;
    Real* tempDeviceYBuffer = NULL;
    Real* tempDeviceZBuffer = NULL;
    Real* tempDeviceAbsBuffer = NULL;
    Real* tempDevicePhaseBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&tempRealBuffer, sizeof(Real) * uiVectorLength * (uiK + 1)));
    checkCudaErrors(cudaMalloc((void**)&tempComplexBuffer, sizeof(QLComplex) * uiComplexLength * (uiK + 1)));
    checkCudaErrors(cudaMalloc((void**)&tempDeviceYBuffer, sizeof(Real) * uiComplexLength * uiK));
    checkCudaErrors(cudaMalloc((void**)&tempDeviceZBuffer, sizeof(Real) * uiComplexLength * uiK));
    checkCudaErrors(cudaMalloc((void**)&tempDeviceAbsBuffer, sizeof(Real) * uiComplexLength * uiK));
    checkCudaErrors(cudaMalloc((void**)&tempDevicePhaseBuffer, sizeof(Real) * uiComplexLength * uiK));

    checkCudaErrors(cudaMemcpy(tempRealBuffer, realList.GetData(), sizeof(Real) * uiVectorLength * (uiK + 1), cudaMemcpyHostToDevice));
    QLComplex* tempComplexHostBuffer = (QLComplex*)malloc(sizeof(QLComplex) * uiComplexLength * (uiK + 1));
    Real* hostY = (Real*)malloc(sizeof(Real) * uiComplexLength * (uiK + 1));
    Real* hostZ = (Real*)malloc(sizeof(Real) * uiComplexLength * (uiK + 1));
    UINT uiVectorPower = MostSignificantPowerTwo(uiComplexLength);

    NormalizeRealToComplex(tempRealBuffer, tempComplexBuffer, uiVectorLength, uiComplexLength, (uiK + 1), fFactor);
    checkCudaErrors(cudaMemcpy(tempComplexHostBuffer, tempComplexBuffer, sizeof(QLComplex) * uiComplexLength * (uiK + 1), cudaMemcpyDeviceToHost));

    TArray<BYTE> bits;
    for (UINT i = 0; i < uiVectorPower; ++i)
    {
        bits.AddItem(static_cast<BYTE>(i));
    }

    TArray<BYTE> idxs;
    for (UINT i = 0; i < uiVectorPower + uiKPower; ++i)
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

    CalculateDegrees(tempComplexBuffer, tempDeviceAbsBuffer, tempDevicePhaseBuffer, uiK, uiVectorPower, tempDeviceYBuffer, tempDeviceZBuffer, hostY, hostZ);
    CalculateDegrees(tempComplexBuffer + uiComplexLength * uiK, tempDeviceAbsBuffer, tempDevicePhaseBuffer, 1, uiVectorPower, tempDeviceYBuffer, tempDeviceZBuffer, hostY + uiK * uiComplexLength, hostZ + uiK * uiComplexLength);

    TArray<Real> allDots;
    Real normalizeFactor = F(1.0);
    Real maxAmplitude = F(-1.0);
    UINT idxOfAmplitude = 0;
    for (UINT j = 0; j < uiK; ++j)
    {
        QLComplex sum = _zeroc;
        for (UINT k = 0; k < uiComplexLength; ++k)
        {
            sum = _cuCaddf(sum, _cuCmulf(_cuConjf(tempComplexHostBuffer[uiComplexLength * j + k]), tempComplexHostBuffer[uiComplexLength * uiK + k]));
        }
        Real fsumSq = __cuCabsSqf(sum);
        allDots.AddItem(sqrt(fsumSq));
        normalizeFactor = normalizeFactor + fsumSq;

        if (fsumSq > maxAmplitude)
        {
            idxOfAmplitude = j;
            maxAmplitude = fsumSq;
        }
    }

    //this is wrong
    QLGate initialState = ExchangeToYZGate(uiKPower, uiVectorPower, hostY, hostZ, FALSE);
    QLGate daggerState = ExchangeToYZGate(uiVectorPower, hostY + uiK * uiComplexLength, hostZ + uiK * uiComplexLength, FALSE);
    daggerState.Dagger();

    initialState.AppendGate(daggerState, bits);

    QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
    for (UINT j = 0; j < uiVectorPower; ++j)
    {
        TArray<BYTE> bitsctr;
        bitsctr.AddItem(static_cast<BYTE>(j));
        initialState.AppendGate(ctrCollapse, bitsctr);
    }

#pragma region simulate

    //test the build of the state
    TArray<SBasicOperation> ops = initialState.GetOperation(initialState.m_lstQubits);
    SIZE_T opssize = ops.Num();

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(initialState.GetQubitCount(), evn);

    LONGLONG veclen = 1LL << initialState.GetQubitCount();
    QLComplex* res = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * veclen));

    syncQuESTEnv(evn);
    vec.stateVec.real[0] = F(1.0);
    vec.stateVec.imag[0] = F(0.0);
    for (LONGLONG line2 = 1; line2 < veclen; ++line2)
    {
        vec.stateVec.real[line2] = F(0.0);
        vec.stateVec.imag[line2] = F(0.0);
    }
    copyStateToGPU(vec);

    probToBuild = F(1.0);
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

    QLMatrix finalstate = ShowStateVectorDetail(res, idxs, FALSE);

    //test the largest state, which is idxOfAmplitude
    probToObtainLargest = F(1.0);
    for (UINT i = 0; i < uiKPower; ++i)
    {
        probToObtainLargest *= static_cast<Real>(collapseToOutcome(vec, i + uiVectorPower, (idxOfAmplitude & (1U << i)) ? 1 : 0));
    }
    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

#pragma endregion

    if (bDebugPrint)
    {
        appPushLogDate(FALSE);
        for (UINT j = 0; j < uiK; ++j)
        {
            Real amplitude = _cuCabsf(finalstate.Get(0, j));
            appGeneral(_T("%d: %f - %f\n"), j, amplitude, allDots[j] / amplitude / sqrt(normalizeFactor));
        }
        appGeneral(_T("\n"));
        appPopLogDate();
        appGeneral(_T("Probability: expected %f, real %f, probability to measure the largest: %d - %f\n"), normalizeFactor / uiK, probToBuild, idxOfAmplitude, sqrt(probToObtainLargest));
    }

    checkCudaErrors(cudaFree(tempRealBuffer));
    checkCudaErrors(cudaFree(tempComplexBuffer));
    checkCudaErrors(cudaFree(tempDeviceYBuffer));
    checkCudaErrors(cudaFree(tempDeviceZBuffer));
    checkCudaErrors(cudaFree(tempDeviceAbsBuffer));
    checkCudaErrors(cudaFree(tempDevicePhaseBuffer));

    free(tempComplexHostBuffer);
    free(hostY);
    free(hostZ);
}

/**
* test the distance
*/
void QLAPI TestDistance(const CCString& sSaveFileName, UINT uiVectorLength, Real fFactor, UINT testPairCount)
{
    UINT uiPaireAll = testPairCount << 1U;
    QLMatrix randomR(uiPaireAll, uiVectorLength);
    randomR.RandomOneReal();
    TArray<Real> realList = randomR.ToVectorRe();

    UINT uiComplexLength = 1 + uiVectorLength / 2;

    Real* tempRealBuffer = NULL;
    QLComplex* tempComplexBuffer = NULL;
    Real* resBuffer = NULL;
    Real* resHostBuffer = reinterpret_cast<Real*>(malloc(sizeof(Real) * uiPaireAll));

    checkCudaErrors(cudaMalloc((void**)&tempRealBuffer, sizeof(Real) * uiVectorLength * uiPaireAll));
    checkCudaErrors(cudaMalloc((void**)&tempComplexBuffer, sizeof(QLComplex) * uiComplexLength * uiPaireAll));
    checkCudaErrors(cudaMalloc((void**)&resBuffer, sizeof(Real) * uiPaireAll));

    checkCudaErrors(cudaMemcpy(tempRealBuffer, realList.GetData(), sizeof(Real) * uiVectorLength * uiPaireAll, cudaMemcpyHostToDevice));
    NormalizeRealToComplex(tempRealBuffer, tempComplexBuffer, uiVectorLength, uiComplexLength, uiPaireAll, fFactor);

    //QLComplex* matrixelement = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * uiComplexLength * uiPaireAll));
    //checkCudaErrors(cudaMemcpy(matrixelement, tempComplexBuffer, sizeof(QLComplex) * uiComplexLength * uiPaireAll, cudaMemcpyDeviceToHost));
    //QLMatrix norm(uiPaireAll, uiComplexLength, matrixelement);

    //randomR.Transpose();
    //randomR.Print(_T("r"));
    //norm.Transpose();
    //norm.Print(_T("c"));

    UINT iBlock1 = testPairCount > _QL_LAUNCH_MAX_THREAD ? Ceil(testPairCount, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = testPairCount > _QL_LAUNCH_MAX_THREAD ? Ceil(testPairCount, iBlock1) : testPairCount;
    _kernelAE_CompareDistance << <iBlock1, iThread1 >> > (tempRealBuffer, tempComplexBuffer, uiVectorLength, uiComplexLength, resBuffer, testPairCount);

    checkCudaErrors(cudaMemcpy(resHostBuffer, resBuffer, sizeof(Real) * uiPaireAll, cudaMemcpyDeviceToHost));
    SaveCSVAR(resHostBuffer, 2, testPairCount, sSaveFileName);

    checkCudaErrors(cudaFree(tempRealBuffer));
    checkCudaErrors(cudaFree(tempComplexBuffer));
    checkCudaErrors(cudaFree(resBuffer));

    free(resHostBuffer);
}

/**
* do not use this
*/
void QLAPI TestProbalityToBuildAmplitude(const CCString& sFileName, UINT uiVectorLength, UINT uiK, Real fFactor, UINT uiTestCount)
{
    UINT w, h;
    TArray<Real> data = ReadCSVAR(sFileName, w, h);
    appGeneral(_T("data in %s : %d x %d\n"), sFileName, w, h);

    assert(w == uiVectorLength);

    Real* realNumbers = (Real*)malloc(sizeof(Real) * uiVectorLength * uiK);
    UINT uiComplexLength = 1 + uiVectorLength / 2;

    Real* tempRealBuffer = NULL;
    QLComplex* tempComplexBuffer = NULL;
    Real* tempDeviceYBuffer = NULL;
    Real* tempDeviceZBuffer = NULL;
    Real* tempDeviceAbsBuffer = NULL;
    Real* tempDevicePhaseBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&tempRealBuffer, sizeof(Real) * w * h));
    checkCudaErrors(cudaMalloc((void**)&tempComplexBuffer, sizeof(QLComplex) * uiComplexLength * h));
    checkCudaErrors(cudaMalloc((void**)&tempDeviceYBuffer, sizeof(Real) * uiComplexLength * uiK));
    checkCudaErrors(cudaMalloc((void**)&tempDeviceZBuffer, sizeof(Real) * uiComplexLength * uiK));
    checkCudaErrors(cudaMalloc((void**)&tempDeviceAbsBuffer, sizeof(Real) * uiComplexLength * uiK));
    checkCudaErrors(cudaMalloc((void**)&tempDevicePhaseBuffer, sizeof(Real) * uiComplexLength * uiK));

    checkCudaErrors(cudaMemcpy(tempRealBuffer, data.GetData(), sizeof(Real) * w * h, cudaMemcpyHostToDevice));
    QLComplex* tempComplexHostBuffer = (QLComplex*)malloc(sizeof(QLComplex) * uiComplexLength * h);
    Real* hostY = (Real*)malloc(sizeof(Real) * uiComplexLength * (uiK + 1));
    Real* hostZ = (Real*)malloc(sizeof(Real) * uiComplexLength * (uiK + 1));
    UINT uiVectorPower = MostSignificantPowerTwo(uiComplexLength);
    UINT uiVectorCountPower = MostSignificantPowerTwo(uiK);

    NormalizeRealToComplex(tempRealBuffer, tempComplexBuffer, w, uiComplexLength, h, fFactor);
    checkCudaErrors(cudaMemcpy(tempComplexHostBuffer, tempComplexBuffer, sizeof(QLComplex) * uiComplexLength * h, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);

    TArray<BYTE> bits;
    for (UINT i = 0; i < uiVectorPower; ++i)
    {
        bits.AddItem(static_cast<BYTE>(i));
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

    for (UINT i = 0; i < uiTestCount; ++i)
    {
        TArray<UINT> toChoose;

        UINT uiIdxStart = i * uiComplexLength * (uiK + 1);

        //for (UINT k = 0; k < uiComplexLength; ++k)
        //{
        //    appGeneral(_T("%f + %f I, "), tempComplexHostBuffer[uiIdxStart + uiComplexLength * 0 + k].x, tempComplexHostBuffer[uiIdxStart + uiComplexLength * 0 + k].y);
        //}
        //appGeneral(_T("\n"));

        //for (UINT k = 0; k < uiComplexLength; ++k)
        //{
        //    appGeneral(_T("%f + %f I, "), tempComplexHostBuffer[uiIdxStart + uiComplexLength * 1 + k].x, tempComplexHostBuffer[uiIdxStart + uiComplexLength * 1 + k].y);
        //}
        //appGeneral(_T("\n"));

        //for (UINT k = 0; k < uiComplexLength; ++k)
        //{
        //    appGeneral(_T("%f + %f I, "), tempComplexHostBuffer[uiIdxStart + uiComplexLength * uiK + k].x, tempComplexHostBuffer[uiIdxStart + uiComplexLength * uiK + k].y);
        //}
        //appGeneral(_T("\n"));

        CalculateDegrees(tempComplexBuffer + uiIdxStart, tempDeviceAbsBuffer, tempDevicePhaseBuffer, uiK, uiVectorPower, tempDeviceYBuffer, tempDeviceZBuffer, hostY, hostZ);
        CalculateDegrees(tempComplexBuffer + uiIdxStart + uiComplexLength * uiK, tempDeviceAbsBuffer, tempDevicePhaseBuffer, 1, uiVectorPower, tempDeviceYBuffer, tempDeviceZBuffer, hostY + uiK * uiComplexLength, hostZ + uiK * uiComplexLength);

        for (UINT j = 0; j < uiK; ++j)
        {
            QLComplex sum = _zeroc;
            for (UINT k = 0; k < uiComplexLength; ++k)
            {
                sum = _cuCaddf(sum, _cuCmulf(_cuConjf(tempComplexHostBuffer[uiIdxStart + uiComplexLength * j + k]), tempComplexHostBuffer[uiIdxStart + uiComplexLength * uiK + k]));
            }
            appGeneral(_T("%f "), _cuCabsf(sum));
        }
        appGeneral(_T("\n"));

        //this is wrong
        QLGate initialState = ExchangeToYZGate(uiVectorCountPower, uiVectorPower, hostY, hostZ, FALSE);
        QLGate daggerState = ExchangeToYZGate(uiVectorPower, hostY + uiK * uiComplexLength, hostZ + uiK * uiComplexLength, FALSE);
        daggerState.Dagger();

        initialState.AppendGate(daggerState, bits);
        
        QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
        for (UINT j = 0; j < uiVectorPower; ++j)
        {
            TArray<BYTE> bitsctr;
            bitsctr.AddItem(static_cast<BYTE>(j));
            initialState.AppendGate(ctrCollapse, bitsctr);
        }

        //test whether the state is good
        QLSimulatorParametersVector param;
        param.m_byQubitCount = static_cast<BYTE>(initialState.m_lstQubits.Num());
        param.m_MasterGate = initialState;
        param.m_bPrint = FALSE;
        QLSimulatorOutputVector out;
        QLSimulatorVector sim;
        sim.Simulate(&param, &out);

        QLMatrix finalstate = ShowStateVectorDetail(out.m_OutputMatrix.HostBuffer(), idxs, FALSE);
        finalstate.Print("res");

        for (UINT j = 0; j < uiK; ++j)
        {
            appGeneral(_T("%f "), _cuCabsf(finalstate.Get(0, j)));
        }
        appGeneral(_T("\n"));

        appGeneral(_T("%f "), out.m_fProbability);
    }

    appPopLogDate();

    checkCudaErrors(cudaFree(tempRealBuffer));
    checkCudaErrors(cudaFree(tempComplexBuffer));
    checkCudaErrors(cudaFree(tempDeviceYBuffer));
    checkCudaErrors(cudaFree(tempDeviceZBuffer));
    checkCudaErrors(cudaFree(tempDeviceAbsBuffer));
    checkCudaErrors(cudaFree(tempDevicePhaseBuffer));

    free(tempComplexHostBuffer);
    free(hostY);
    free(hostZ);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================