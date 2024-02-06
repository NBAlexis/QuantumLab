//=============================================================================
// FILENAME : FRGate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [30/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

#pragma region kernels

__device__ __inline__ UINT _deviceGrayCode(UINT n)
{
    return n ^ (n >> 1);
}

__device__  __inline__ UINT _deviceCounting1Bit(UINT n)
{
    UINT c;
    for (c = 0; n; ++c)
    {
        n &= n - 1;
    }
    return c;
}

__device__  __inline__ UINT _deviceBitWiseInnerProduct(UINT a, UINT b)
{
    return _deviceCounting1Bit(a & b);
}

//Why is it slower than CPU!!!
__global__ void _QL_LAUNCH_BOUND
_kernelSplitAngles(
    const Real* __restrict__ angles,
    Real* res,
    UINT i,
    UINT uiLength)
{
    const UINT j = (threadIdx.x + blockIdx.x * blockDim.x);
    if (j < uiLength)
    {
        UINT uiSign = _deviceBitWiseInnerProduct(j, _deviceGrayCode(i));
        res[j] = (uiSign & 1) ? (-angles[j]) : (angles[j]);
    }
}

//__global__ void _QL_LAUNCH_BOUND
//_kernelSplitAngles(
//    const Real* __restrict__ angles,
//    Real* res,
//    UINT uiLength)
//{
//    const UINT i = (threadIdx.x + blockIdx.x * blockDim.x);
//    Real theta = F(0.0);
//    for (UINT j = 0; j < uiLength; ++j)
//    {
//        UINT uiSign = _deviceBitWiseInnerProduct(j, _deviceGrayCode(i));
//        theta = theta + ((uiSign & 1) ? (-angles[j]) : (angles[j]));
//    }
//    res[i] = theta / uiLength;
//}

#pragma endregion

CSpliteAngleBufferHelper*  _splitAnglePointer = NULL;

CSpliteAngleBufferHelper::CSpliteAngleBufferHelper()
    : m_bHasBuffer(FALSE)
    , m_pDeviceBufferConst(NULL)
    , m_pDeviceBuffer(NULL)
    , m_pHostRes(NULL)
{
    _splitAnglePointer = this;
    InitialBuffer();
}

CSpliteAngleBufferHelper::~CSpliteAngleBufferHelper()
{
    if (NULL != m_pDeviceBufferConst)
    {
        checkCudaErrors(cudaFree(m_pDeviceBufferConst));
        m_pDeviceBufferConst = NULL;
    }
    if (NULL != m_pDeviceBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceBuffer));
        m_pDeviceBuffer = NULL;
    }
    if (NULL != m_pHostRes)
    {
        appSafeFree(m_pHostRes);
    }
}

void CSpliteAngleBufferHelper::InitialBuffer()
{
    if (m_bHasBuffer)
    {
        return;
    }
    m_bHasBuffer = TRUE;
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceBufferConst, sizeof(Real) *
        (1UL << static_cast<UINT>(_kMaxSupport))
    ));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceBuffer, sizeof(Real) *
        (1UL << static_cast<UINT>(_kMaxSupport))
    ));
    //m_pHostRes = reinterpret_cast<Real*>(malloc(sizeof(Real) *
    //    (1UL << static_cast<UINT>(_kMaxSupport))
    //));
}

TArray<Real> CSpliteAngleBufferHelper::SpliteAngles(const Real* angles, UINT length)
{
    InitialBuffer();

    TArray<Real> retList;
    checkCudaErrors(cudaMemcpy(m_pDeviceBufferConst, angles, sizeof(Real) * length, cudaMemcpyHostToDevice));
    UINT uiBlock = 0;
    UINT uiThread = 0;
    __DECOMPOSE(length, uiBlock, uiThread);

    for (UINT i = 0; i < length; ++i)
    {
        _kernelSplitAngles << <uiBlock, uiThread >> > (m_pDeviceBufferConst, m_pDeviceBuffer, i, length);
        retList.AddItem(ReduceSum(m_pDeviceBuffer, length) / length);
    }

    return retList;
}

//TArray<Real> CSpliteAngleBufferHelper::SpliteAngles(const Real* angles, UINT length)
//{
//    InitialBuffer();
//
//    TArray<Real> retList;
//    checkCudaErrors(cudaMemcpy(m_pDeviceBufferConst, angles, sizeof(Real) * length, cudaMemcpyHostToDevice));
//    UINT uiBlock = 0;
//    UINT uiThread = 0;
//    __DECOMPOSE(length, uiBlock, uiThread);
//    _kernelSplitAngles << <uiBlock, uiThread >> > (m_pDeviceBufferConst, m_pDeviceBuffer, length);
//    checkCudaErrors(cudaMemcpy(m_pHostRes, m_pDeviceBuffer, sizeof(Real) * length, cudaMemcpyDeviceToHost));
//    retList.Append(m_pHostRes, length);
//
//    return retList;
//}

TArray<Real> QLAPI SpliteAngles(const Real* angles, UINT length)
{
    if (length > 2048 && NULL != _splitAnglePointer)
    {
        return _splitAnglePointer->SpliteAngles(angles, length);
    }

    TArray<Real> retList;
    for (UINT i = 0; i < length; ++i)
    {
        Real theta = F(0.0);
        for (UINT j = 0; j < length; ++j)
        {
            UINT uiSign = BitWiseInnerProduct(j, GrayCode(i));
            theta = theta + ((uiSign & 1) ? (-angles[j]) : (angles[j]));
        }
        retList.AddItem(theta / length);
    }
    return retList;
}

TArray<Real> QLAPI SpliteAngles(const TArray<Real>& angles, UINT length)
{
    return SpliteAngles(angles.GetData(), length);
}

QLGate QLAPI FRy(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_RY, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRy");

    TArray <Real> theta = SpliteAngles(angles, degreeNumber);
    TArray <BYTE> target;
    target.AddItem(static_cast<BYTE>(numberOfQubits - 1));
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, theta[i]);
        ret.AppendGate(ry, target);

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        ctrIdx = numberOfQubits - 2 - ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
        cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }

    return ret;
}

QLGate QLAPI FRy(const Real* angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_RY, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRy");

    TArray <Real> theta = SpliteAngles(angles, degreeNumber);
    TArray <BYTE> target;
    target.AddItem(static_cast<BYTE>(numberOfQubits - 1));
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, theta[i]);
        ret.AppendGate(ry, target);

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        ctrIdx = numberOfQubits - 2 - ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
        cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }

    return ret;
}

QLGate QLAPI FRz(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_RZ, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRz");

    TArray <Real> theta = SpliteAngles(angles, degreeNumber);
    TArray <BYTE> target;
    target.AddItem(static_cast<BYTE>(numberOfQubits - 1));
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate rz(EBasicOperation::EBO_RZ, theta[i]);
        ret.AppendGate(rz, target);

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        ctrIdx = numberOfQubits - 2 - ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
        cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }

    return ret;
}

QLGate QLAPI FRz(const Real* angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_RZ, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRz");

    TArray <Real> theta = SpliteAngles(angles, degreeNumber);
    TArray <BYTE> target;
    target.AddItem(static_cast<BYTE>(numberOfQubits - 1));
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate rz(EBasicOperation::EBO_RZ, theta[i]);
        ret.AppendGate(rz, target);

        UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
        ctrIdx = numberOfQubits - 2 - ctrIdx;

        TArray <BYTE> cnotbits;
        cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
        cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
        ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
    }

    return ret;
}

QLGate QLAPI FRyz(const TArray<Real>& anglesY, const TArray<Real>& anglesZ, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(anglesY.Num()) < degreeNumber || static_cast<UINT>(anglesZ.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        QLGate rzdagger(EBasicOperation::EBO_RZ, -anglesZ[0]);
        //rzdagger.Dagger();
        QLGate ry(EBasicOperation::EBO_RY, anglesY[0]);

        QLGate ret;
        ret.AddQubits(1);
        ret.m_sName = _T("FRyz");
        TArray<BYTE> lstQubits;
        lstQubits.AddItem(0);

        ret.AppendGate(ry, lstQubits);
        ret.AppendGate(rzdagger, lstQubits);
        return ret;
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRyz");

    TArray <Real> thetaY = SpliteAngles(anglesY, degreeNumber);
    TArray <Real> thetaZ = SpliteAngles(anglesZ, degreeNumber);

    TArray <BYTE> target;
    target.AddItem(static_cast<BYTE>(numberOfQubits - 1));
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, thetaY[i]);
        ret.AppendGate(ry, target);

        if (i != degreeNumber - 1)
        {
            UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;
            TArray <BYTE> cnotbits;
            cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
            cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }
    }

    for (UINT i = 0; i < degreeNumber; ++i)
    {
        UINT j = degreeNumber - 1 - i;
        if (0 != i)
        {
            UINT ctrIdx = GrayCodeDifferent(j, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;

            TArray <BYTE> cnotbits;
            cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
            cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }

        QLGate rz(EBasicOperation::EBO_RZ, -thetaZ[j]);
        ret.AppendGate(rz, target);
    }

    return ret;
}

QLGate QLAPI FRyz(const Real* anglesY, const Real* anglesZ, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);

    if (1 == numberOfQubits)
    {
        QLGate rzdagger(EBasicOperation::EBO_RZ, -anglesZ[0]);
        //rzdagger.Dagger();
        QLGate ry(EBasicOperation::EBO_RY, anglesY[0]);

        QLGate ret;
        ret.AddQubits(1);
        ret.m_sName = _T("FRyz");
        TArray<BYTE> lstQubits;
        lstQubits.AddItem(0);

        ret.AppendGate(ry, lstQubits);
        ret.AppendGate(rzdagger, lstQubits);
        return ret;
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRyz");

    TArray <Real> thetaY = SpliteAngles(anglesY, degreeNumber);
    TArray <Real> thetaZ = SpliteAngles(anglesZ, degreeNumber);

    TArray <BYTE> target;
    target.AddItem(static_cast<BYTE>(numberOfQubits - 1));
    for (UINT i = 0; i < degreeNumber; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, thetaY[i]);
        ret.AppendGate(ry, target);

        if (i != degreeNumber - 1)
        {
            UINT ctrIdx = GrayCodeDifferent(i, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;
            TArray <BYTE> cnotbits;
            cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
            cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }
    }

    for (UINT i = 0; i < degreeNumber; ++i)
    {
        UINT j = degreeNumber - 1 - i;
        if (0 != i)
        {
            UINT ctrIdx = GrayCodeDifferent(j, degreeNumber);
            ctrIdx = numberOfQubits - 2 - ctrIdx;

            TArray <BYTE> cnotbits;
            cnotbits.AddItem(static_cast<BYTE>(ctrIdx));
            cnotbits.AddItem(static_cast<BYTE>(numberOfQubits - 1));
            ret.AppendGate(QLGate(EBasicOperation::EBO_CX), cnotbits);
        }

        QLGate rz(EBasicOperation::EBO_RZ, -thetaZ[j]);
        ret.AppendGate(rz, target);
    }

    return ret;
}

QLGate QLAPI FRp(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_P, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRp");

    QLGate x(EBasicOperation::EBO_X);
    TArray<BYTE> allBits;
    for (UINT i = 0; i < numberOfQubits - 1; ++i)
    {
        TArray<BYTE> targetbit;
        targetbit.AddItem(static_cast<BYTE>(i));
        ret.AppendGate(x, targetbit);

        allBits.AddItem(static_cast<BYTE>(i));
    }
    allBits.AddItem(static_cast<BYTE>(numberOfQubits - 1));

    QLGate cnp = CreateCnP(static_cast<BYTE>(numberOfQubits - 1), angles[0]);
    ret.AppendGate(cnp, allBits);

    for (UINT i = 1; i < degreeNumber; ++i)
    {
        UINT diff = i ^ (i - 1);
        for (UINT j = 0; j < numberOfQubits - 1; ++j)
        {
            if (diff & (1 << j))
            {
                TArray<BYTE> cnotbit;
                cnotbit.AddItem(static_cast<BYTE>(j));
                ret.AppendGate(x, cnotbit);
            }
        }

        QLGate cnp2 = CreateCnP(static_cast<BYTE>(numberOfQubits - 1), angles[i]);
        ret.AppendGate(cnp2, allBits);
    }
    return ret;
}

QLGate QLAPI FRPh(const TArray<Real>& angles, UINT numberOfQubits)
{
    assert(numberOfQubits >= 1);
    UINT degreeNumber = 1U << (numberOfQubits - 1);
    if (static_cast<UINT>(angles.Num()) < degreeNumber)
    {
        appCrucial(_T("degree number wrong!\n"));
        return QLGate();
    }

    if (1 == numberOfQubits)
    {
        return QLGate(EBasicOperation::EBO_Phase, angles[0]);
    }

    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(numberOfQubits));
    ret.m_sName = _T("FRPh");

    QLGate x(EBasicOperation::EBO_X);
    TArray<BYTE> allBits;
    for (UINT i = 0; i < numberOfQubits - 1; ++i)
    {
        TArray<BYTE> targetbit;
        targetbit.AddItem(static_cast<BYTE>(i));
        ret.AppendGate(x, targetbit);

        allBits.AddItem(static_cast<BYTE>(i));
    }
    allBits.AddItem(static_cast<BYTE>(numberOfQubits - 1));

    QLGate cnp = CreateCnPh(static_cast<BYTE>(numberOfQubits - 1), angles[0]);
    ret.AppendGate(cnp, allBits);

    for (UINT i = 1; i < degreeNumber; ++i)
    {
        UINT diff = i ^ (i - 1);
        for (UINT j = 0; j < numberOfQubits - 1; ++j)
        {
            if (diff & (1 << j))
            {
                TArray<BYTE> cnotbit;
                cnotbit.AddItem(static_cast<BYTE>(j));
                ret.AppendGate(x, cnotbit);
            }
        }

        QLGate cnp2 = CreateCnPh(static_cast<BYTE>(numberOfQubits - 1), angles[i]);
        ret.AppendGate(cnp2, allBits);
    }
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================