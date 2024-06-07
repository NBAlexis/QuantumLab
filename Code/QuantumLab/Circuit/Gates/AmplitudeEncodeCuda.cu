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

__global__ void _QL_LAUNCH_BOUND
_kernelAE_VtoAbs(const QLComplex* __restrict__ v, Real* absBuffer, UINT uiMax)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax)
    {
        absBuffer[idx] = v[idx].x;
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
        UINT putStart = uiVLen - (1U << (lengthPower - uiLevel));
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

__global__ void _QL_LAUNCH_BOUND
_kernelAE_RotationsReal(Real* YLst, Real* absBuffer, UINT lengthPower, UINT vectorCount, UINT uiLevel)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    const UINT uiVLen = 1U << lengthPower;

    if (lengthPower == uiLevel)
    {
        //only record the total phase
        return;
    }

    const UINT numberOfDegrees = 1U << (lengthPower - uiLevel - 1);
    const UINT uiIdOfV = idx / numberOfDegrees;

    if (uiIdOfV < vectorCount)
    {
        const UINT uiIdInV = idx % numberOfDegrees;

        const UINT uiLeft = uiIdOfV * uiVLen + uiIdInV * (1U << (uiLevel + 1));
        const UINT uiRight = uiLeft + (1U << uiLevel);
        UINT putStart = uiVLen - (1U << (lengthPower - uiLevel));
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
    }
}

/**
* similar as _kernelAE_Rotations
* But, _kernelAE_Rotations is to build |phi_i>|i>
* it put the degrees as
* 
* d1 d1 d1 d1 D1 D1 D1 D1, d2 d2 D2 D2, d3 D3, x
* 
* where d1 is d1 of v1, D1 is d1 of v2.
* 
* This function put it as:
* v1:
* d1 d1 d1 d1, d2 d2, d3, x
* v2:
* D1 D1 D1 D1, D2 D2, D3, x
* ...
*/
__global__ void _QL_LAUNCH_BOUND
_kernelAE_RotationsForEachVector(Real* YLst, Real* ZLst, Real* absBuffer, Real* phaseBuffer, UINT lengthPower, UINT vectorCount, UINT uiLevel)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    const UINT uiVLen = 1U << lengthPower;

    if (lengthPower == uiLevel && NULL != ZLst && NULL != phaseBuffer)
    {
        //only record the total phase
        ZLst[idx * uiVLen + (uiVLen - 1)] = phaseBuffer[uiVLen * idx];
        return;
    }

    const UINT numberOfDegrees = 1U << (lengthPower - uiLevel - 1);
    const UINT uiIdOfV = idx / numberOfDegrees;

    if (uiIdOfV < vectorCount)
    {
        const UINT uiIdInV = idx % numberOfDegrees;

        const UINT uiLeft = uiIdOfV * uiVLen + uiIdInV * (1U << (uiLevel + 1));
        const UINT uiRight = uiLeft + (1U << uiLevel);
        UINT putStart = uiVLen - (1U << (lengthPower - uiLevel));
        putStart = uiIdOfV * uiVLen + putStart + uiIdInV;

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

__global__ void _QL_LAUNCH_BOUND
_kernelAE_RotationsRealForEachVector(Real* YLst, Real* absBuffer, UINT lengthPower, UINT vectorCount, UINT uiLevel)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
    const UINT uiVLen = 1U << lengthPower;

    if (lengthPower == uiLevel)
    {
        //only record the total phase
        return;
    }

    const UINT numberOfDegrees = 1U << (lengthPower - uiLevel - 1);
    const UINT uiIdOfV = idx / numberOfDegrees;

    if (uiIdOfV < vectorCount)
    {
        const UINT uiIdInV = idx % numberOfDegrees;

        const UINT uiLeft = uiIdOfV * uiVLen + uiIdInV * (1U << (uiLevel + 1));
        const UINT uiRight = uiLeft + (1U << uiLevel);
        UINT putStart = uiVLen - (1U << (lengthPower - uiLevel));
        //putStart = vectorCount * putStart + uiIdOfV * numberOfDegrees + uiIdInV;
        putStart = uiIdOfV * uiVLen + putStart + uiIdInV;

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

void QLAPI CalculateDegrees(Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ)
{
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

void QLAPI CalculateDegreesForEach(Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ)
{
    for (UINT i = 0; i <= vectorPower; ++i)
    {
        if (i != vectorPower)
        {
            UINT uiLen2 = vectorCount * (1U << (vectorPower - i - 1));
            UINT iBlock2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, _QL_LAUNCH_MAX_THREAD) : 1;
            UINT iThread2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, iBlock2) : uiLen2;
            _kernelAE_RotationsForEachVector << <iBlock2, iThread2 >> > (deviceY, deviceZ, absBuffer, phaseBuffer, vectorPower, vectorCount, i);
        }
        else
        {
            UINT uiLen2 = vectorCount;
            UINT iBlock2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, _QL_LAUNCH_MAX_THREAD) : 1;
            UINT iThread2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, iBlock2) : uiLen2;
            _kernelAE_RotationsForEachVector << <iBlock2, iThread2 >> > (deviceY, deviceZ, absBuffer, phaseBuffer, vectorPower, vectorCount, i);
        }
    }
}

void QLAPI CalculateDegreesReal(Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY)
{
    for (UINT i = 0; i < vectorPower; ++i)
    {
        UINT uiLen2 = vectorCount * (1U << (vectorPower - i - 1));
        UINT iBlock2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, iBlock2) : uiLen2;
        _kernelAE_RotationsReal << <iBlock2, iThread2 >> > (deviceY, absBuffer, vectorPower, vectorCount, i);
    }
}

void QLAPI CalculateDegreesRealForEach(Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY)
{
    for (UINT i = 0; i < vectorPower; ++i)
    {
        UINT uiLen2 = vectorCount * (1U << (vectorPower - i - 1));
        UINT iBlock2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, _QL_LAUNCH_MAX_THREAD) : 1;
        UINT iThread2 = uiLen2 > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen2, iBlock2) : uiLen2;
        _kernelAE_RotationsRealForEachVector << <iBlock2, iThread2 >> > (deviceY, absBuffer, vectorPower, vectorCount, i);
    }
}

void QLAPI CalculateDegrees(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ)
{
    UINT vectorLength = 1U << vectorPower;
    UINT uiLen = vectorCount * vectorLength;
    UINT iBlock1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, iBlock1) : uiLen;
    _kernelAE_VtoAbsAndPhase << <iBlock1, iThread1 >> > (deviceV, absBuffer, phaseBuffer, uiLen);
    CalculateDegrees(absBuffer, phaseBuffer, vectorCount, vectorPower, deviceY, deviceZ);
}

void QLAPI CalculateDegreesForEach(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ)
{
    UINT vectorLength = 1U << vectorPower;
    UINT uiLen = vectorCount * vectorLength;
    UINT iBlock1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, iBlock1) : uiLen;
    _kernelAE_VtoAbsAndPhase << <iBlock1, iThread1 >> > (deviceV, absBuffer, phaseBuffer, uiLen);
    CalculateDegreesForEach(absBuffer, phaseBuffer, vectorCount, vectorPower, deviceY, deviceZ);
}

void QLAPI CalculateDegreesReal(const QLComplex* deviceV, Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY)
{
    UINT vectorLength = 1U << vectorPower;
    UINT uiLen = vectorCount * vectorLength;
    UINT iBlock1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = uiLen > _QL_LAUNCH_MAX_THREAD ? Ceil(uiLen, iBlock1) : uiLen;
    _kernelAE_VtoAbs << <iBlock1, iThread1 >> > (deviceV, absBuffer, uiLen);
    CalculateDegreesReal(absBuffer, vectorCount, vectorPower, deviceY);
}

void QLAPI CalculateDegrees(const QLComplex* deviceV, Real* absBuffer, Real* phaseBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* deviceZ, Real* hostY, Real* hostZ)
{
    CalculateDegrees(deviceV, absBuffer, phaseBuffer, vectorCount, vectorPower, deviceY, deviceZ);
    UINT uiLen = vectorCount * (1U << vectorPower);

    checkCudaErrors(cudaMemcpy(hostY, deviceY, sizeof(Real) * uiLen, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(hostZ, deviceZ, sizeof(Real) * uiLen, cudaMemcpyDeviceToHost));
}

void QLAPI CalculateDegreesReal(const QLComplex* deviceV, Real* absBuffer, UINT vectorCount, UINT vectorPower, Real* deviceY, Real* hostY)
{
    CalculateDegreesReal(deviceV, absBuffer, vectorCount, vectorPower, deviceY);
    UINT uiLen = vectorCount * (1U << vectorPower);

    checkCudaErrors(cudaMemcpy(hostY, deviceY, sizeof(Real) * uiLen, cudaMemcpyDeviceToHost));
}

/**
* the degrees is
* 00000000 1111 22 3 4

* yz(00000000) at 0123
* yz(1111) at 012
* yz(22) at 01
* we need yz(3) at 0
*
* and then global phase at 0
* 
* |i> use higher qubits, and phi use lower qubits
* For example, 8 phi_i with len(phi) = 4, there are 5 qubits, |i> is 2,3,4 and |phi> is 0,1
*/
QLGate QLAPI ExchangeToYZGate(UINT vectorPower, Real* hostY, Real* hostZ, UBOOL bHasPhase)
{
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
* If the vector is real vector
*/
QLGate QLAPI ExchangeToYGate(UINT vectorPower, Real* hostY)
{
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
            retGate.AppendGate(ry, bits);
            idx = idx - 2;
        }
        else
        {
            QLGate fry = FRy(hostY + idx, static_cast<UINT>(bits.Num()));
            idx = idx - (1U << (i + 1));
            retGate.AppendGate(fry, bits);
        }
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
    for (UINT i = 0; i < phaseCount; ++i)
    {
        a0 += phases[i];
    }
    a0 = a0 / phaseCount;
    degreeZ[phaseCount - 1] = a0;

    UINT idxStart = phaseCount - 2;
    for (UINT i = 0; i < lengthPower; ++i)
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
        retGate.m_lstAdditionalQubitsData.AddItem(static_cast<BYTE>(i + vectorPower));
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

QLGate QLAPI ExchangeToYGate(UINT vectorCountPower, UINT vectorPower, Real* hostY)
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
        retGate.m_lstAdditionalQubitsData.AddItem(static_cast<BYTE>(i + vectorPower));
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

        QLGate fry = FRy(hostY + startInV, static_cast<UINT>(bits.Num()));

        retGate.AppendGate(fry, bits);
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

QLGate QLAPI AmplitudeEncodeOneVectorReal(const QLComplex* hostv, UINT vectorPower)
{
    UINT vLength = 1U << vectorPower;

    QLComplex* v = NULL;
    Real* absBuffer = NULL;
    Real* YBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&v, sizeof(QLComplex) * vLength));
    checkCudaErrors(cudaMalloc((void**)&absBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMalloc((void**)&YBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMemcpy(v, hostv, sizeof(QLComplex) * vLength, cudaMemcpyHostToDevice));

    Real* YHostBuffer = (Real*)malloc(sizeof(Real) * vLength);

    CalculateDegreesReal(v, absBuffer, 1, vectorPower, YBuffer, YHostBuffer);

    checkCudaErrors(cudaFree(v));
    checkCudaErrors(cudaFree(absBuffer));
    checkCudaErrors(cudaFree(YBuffer));

    QLGate ret = ExchangeToYGate(vectorPower, YHostBuffer);

    free(YHostBuffer);

    return ret;
}

QLGate QLAPI AmplitudeEncodeOneVectorReal(const Real* hostv, UINT vectorPower)
{
    UINT vLength = 1U << vectorPower;

    Real* absBuffer = NULL;
    Real* YBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&absBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMalloc((void**)&YBuffer, sizeof(Real) * vLength));
    checkCudaErrors(cudaMemcpy(absBuffer, hostv, sizeof(Real) * vLength, cudaMemcpyHostToDevice));

    Real* YHostBuffer = (Real*)malloc(sizeof(Real) * vLength);

    CalculateDegreesReal(absBuffer, 1, vectorPower, YBuffer);
    checkCudaErrors(cudaMemcpy(YHostBuffer, YBuffer, sizeof(Real) * vLength, cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(absBuffer));
    checkCudaErrors(cudaFree(YBuffer));

    QLGate ret = ExchangeToYGate(vectorPower, YHostBuffer);

    free(YHostBuffer);

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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================