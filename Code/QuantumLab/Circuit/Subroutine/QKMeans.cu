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
_kernelQKMeansSpliteK(BYTE* kBuffer, UINT uiMax, BYTE kToSplit, BYTE* splitTo, BYTE countOfSplitTo)
{
    UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < uiMax && kBuffer[idx] == kToSplit)
    {
        kBuffer[idx] = splitTo[__r->_deviceRandomI(threadIdx.x, countOfSplitTo)];
    }
}

#pragma endregion

/**
* 
def ApplyVectorList(vv, circuit: QuantumCircuit, qubits: list):
    retlsty1, retlsty2, retlstz1, retlstz2 = Calculate4DegreeList(vv)
    qubits1 = qubits[2:len(qubits)]
    [5,4,3,2,1]
    [5,4,3,2,1,0]
    qubits2 = qubits[1:len(qubits)]
    for q in qubits1:
        circuit.h(q)
    qubits1.reverse()
    qubits2.reverse()
    FRy(np.array(retlsty1) * 2, circuit, qubits[1], qubits1)
    FRy(np.array(retlsty2) * 2, circuit, qubits[0], qubits2)
    FRz(np.array(retlstz1) * 2, circuit, qubits[0], qubits2)
    FRz(np.array(retlstz2) * 2, circuit, qubits[1], qubits1)

*/


QLQuantumKmeans::QLQuantumKmeans(BYTE maxK)
    : m_byMaxK(maxK)
    , m_byQubit(0)
    , m_byVectorCount(maxK + 1)
    , m_uiRepeat(1)
    , m_uiN(0)
    , m_uiBlock(1)
    , m_uiThread(1)
    , m_uiBlockN(1)
    , m_uiThreadN(1)
    , m_uiBlockC(1)
    , m_uiThreadC(1)

    , m_pDeviceVBuffer(NULL)
    , m_pDeviceCVBuffer(NULL)
    , m_pDeviceKCounts(NULL)
    , m_pHostKCounts(NULL)
    , m_pDeviceTempKBuffer(NULL)

    , m_pDeviceY1Buffer(NULL)
    , m_pDeviceY2Buffer(NULL)
    , m_pDeviceZ1Buffer(NULL)
    , m_pDeviceZ2Buffer(NULL)
    , m_pHostY1Buffer(NULL)
    , m_pHostY2Buffer(NULL)
    , m_pHostZ1Buffer(NULL)
    , m_pHostZ2Buffer(NULL)

    , m_pDeviceData(NULL)
    , m_pHostDataY1Buffer(NULL)
    , m_pHostDataY2Buffer(NULL)
    , m_pHostDataZ1Buffer(NULL)
    , m_pHostDataZ2Buffer(NULL)
    , m_pHostKValues(NULL)
    , m_pDeviceKValues(NULL)

    , m_pDeviceRealWorkingBuffer(NULL)
    , m_pDeviceUIntWorkingBuffer(NULL)

    , m_pHostCenters(NULL)
    , m_uiStep(0)
{
    m_byQubit = 2 + static_cast<BYTE>(MostSignificantPowerTwo(maxK));

    m_uiBlock = m_byVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(m_byVectorCount, _QL_LAUNCH_MAX_THREAD) : 1;
    m_uiThread = m_byVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(m_byVectorCount, m_uiBlock) : m_byVectorCount;
    m_uiBlockC = m_byVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(m_byMaxK, _QL_LAUNCH_MAX_THREAD) : 1;
    m_uiThreadC = m_byVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(m_byMaxK, m_uiBlockC) : m_byMaxK;

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceVBuffer, sizeof(Real) * 7 * m_byVectorCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceCVBuffer, sizeof(QLComplex) * 4 * m_byVectorCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceKCounts, sizeof(UINT) * m_byMaxK));
    m_pHostKCounts = (UINT*)malloc(sizeof(UINT) * m_byMaxK);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceTempKBuffer, sizeof(BYTE) * m_byMaxK));

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceY1Buffer, sizeof(Real) * m_byVectorCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceZ1Buffer, sizeof(Real) * m_byVectorCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceY2Buffer, sizeof(Real) * 2 * m_byVectorCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceZ2Buffer, sizeof(Real) * 2 * m_byVectorCount));

    m_pHostY1Buffer = (Real*)malloc(sizeof(Real) * m_byVectorCount);
    m_pHostZ1Buffer = (Real*)malloc(sizeof(Real) * m_byVectorCount);
    m_pHostY2Buffer = (Real*)malloc(sizeof(Real) * 2 * m_byVectorCount);
    m_pHostZ2Buffer = (Real*)malloc(sizeof(Real) * 2 * m_byVectorCount);

}

QLQuantumKmeans::~QLQuantumKmeans()
{
    checkCudaErrors(cudaFree(m_pDeviceVBuffer));
    checkCudaErrors(cudaFree(m_pDeviceCVBuffer));
    checkCudaErrors(cudaFree(m_pDeviceY1Buffer));
    checkCudaErrors(cudaFree(m_pDeviceZ1Buffer));
    checkCudaErrors(cudaFree(m_pDeviceY2Buffer));
    checkCudaErrors(cudaFree(m_pDeviceZ2Buffer));

    appSafeFree(m_pHostY1Buffer);
    appSafeFree(m_pHostZ1Buffer);
    appSafeFree(m_pHostY2Buffer);
    appSafeFree(m_pHostZ2Buffer);
}


void QLQuantumKmeans::Prepare(const CCString& fileName, const CCString& sStartCenterFile, UINT n)
{
    m_sStartCenterFile = sStartCenterFile;

    UINT w, h;
    TArray<Real> data = ReadCSVAR(fileName, w, h);

    appGeneral(_T("File loaded %d x %d.\n"), w, h);

    m_uiN = (0 == n ? h : n);
    m_uiBlockN = m_uiN > _QL_LAUNCH_MAX_THREAD ? Ceil(m_uiN, _QL_LAUNCH_MAX_THREAD) : 1;
    m_uiThreadN = m_uiN > _QL_LAUNCH_MAX_THREAD ? Ceil(m_uiN, m_uiBlockN) : m_uiN;

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(Real) * 7 * m_uiN));
    checkCudaErrors(cudaMemcpy(m_pDeviceData, data.GetData(), sizeof(Real) * 7 * m_uiN, cudaMemcpyHostToDevice));

    m_pHostKValues = (BYTE*)malloc(sizeof(BYTE) * m_uiN);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceKValues, sizeof(BYTE) * m_uiN));

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRealWorkingBuffer, sizeof(Real) * m_uiN * 7));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceUIntWorkingBuffer, sizeof(UINT) * m_uiN));

    m_pHostCenters = (Real*)(malloc(sizeof(Real) * 7 * m_byMaxK));

    //============== initial degree buffers ===============
    m_pHostDataY1Buffer = (Real*)malloc(sizeof(Real) * m_uiN);
    m_pHostDataY2Buffer = (Real*)malloc(sizeof(Real) * 2 * m_uiN);
    m_pHostDataZ1Buffer = (Real*)malloc(sizeof(Real) * m_uiN);
    m_pHostDataZ2Buffer = (Real*)malloc(sizeof(Real) * 2 * m_uiN);

    QLComplex* pTempComplexBuffer = NULL;
    Real* pTempY1 = NULL;
    Real* pTempY2 = NULL;
    Real* pTempZ1 = NULL;
    Real* pTempZ2 = NULL;
    checkCudaErrors(cudaMalloc((void**)&pTempComplexBuffer, sizeof(QLComplex) * 4 * m_uiN));
    checkCudaErrors(cudaMalloc((void**)&pTempY1, sizeof(Real) * m_uiN));
    checkCudaErrors(cudaMalloc((void**)&pTempY2, sizeof(Real) * 2 * m_uiN));
    checkCudaErrors(cudaMalloc((void**)&pTempZ1, sizeof(Real) * m_uiN));
    checkCudaErrors(cudaMalloc((void**)&pTempZ2, sizeof(Real) * 2 * m_uiN));
    _kernelQKMADVectorToNormalizedVector << <m_uiBlockN, m_uiThreadN >> > (pTempComplexBuffer, m_pDeviceData, m_uiN);
    _kernelQKMADVectorToAngle << <m_uiBlockN, m_uiThreadN >> > (pTempY1, pTempY2, pTempZ1, pTempZ2, pTempComplexBuffer, m_uiN);
    checkCudaErrors(cudaMemcpy(m_pHostDataY1Buffer, pTempY1, sizeof(Real) * m_uiN, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostDataZ1Buffer, pTempZ1, sizeof(Real) * m_uiN, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostDataY2Buffer, pTempY2, sizeof(Real) * 2 * m_uiN, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostDataZ2Buffer, pTempZ2, sizeof(Real) * 2 * m_uiN, cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(pTempComplexBuffer));
    checkCudaErrors(cudaFree(pTempY1));
    checkCudaErrors(cudaFree(pTempY2));
    checkCudaErrors(cudaFree(pTempZ1));
    checkCudaErrors(cudaFree(pTempZ2));
    m_pMeasureCounts = (UINT*)malloc(sizeof(UINT) * m_uiN);

    appGeneral(_T("Buffer prepared.\n"));
}

/**
* host vector buffer is vectors v[i] + w, for example, 64 v[i] + 1 w.
* in this case, vectorCount = 65, qubitCount = 2 + log2(64) = 8
* dimension of host buffer is 7, and are expected to be already initialized
*/
void QLQuantumKmeans::CalculateDegrees()
{
    //checkCudaErrors(cudaMemcpy(m_pDeviceVBuffer, hostVectorBuffer, sizeof(Real) * 7 * m_byVectorCount, cudaMemcpyHostToDevice));

    _kernelQKMADVectorToNormalizedVector << <m_uiBlock, m_uiThread >> > (m_pDeviceCVBuffer, m_pDeviceVBuffer, m_byVectorCount);
    _kernelQKMADVectorToAngle << <m_uiBlock, m_uiThread >> > (m_pDeviceY1Buffer, m_pDeviceY2Buffer, m_pDeviceZ1Buffer, m_pDeviceZ2Buffer, m_pDeviceCVBuffer, m_byVectorCount);

    checkCudaErrors(cudaMemcpy(m_pHostY1Buffer, m_pDeviceY1Buffer, sizeof(Real) * m_byVectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostZ1Buffer, m_pDeviceZ1Buffer, sizeof(Real) * m_byVectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostY2Buffer, m_pDeviceY2Buffer, sizeof(Real) * 2 * m_byVectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostZ2Buffer, m_pDeviceZ2Buffer, sizeof(Real) * 2 * m_byVectorCount, cudaMemcpyDeviceToHost));
}

void QLQuantumKmeans::CalculateDegreesOnlyCenters()
{
    _kernelQKMADVectorToNormalizedVector << <m_uiBlockC, m_uiThreadC >> > (m_pDeviceCVBuffer, m_pDeviceVBuffer, m_byMaxK);
    _kernelQKMADVectorToAngle << <m_uiBlockC, m_uiThreadC >> > (m_pDeviceY1Buffer, m_pDeviceY2Buffer, m_pDeviceZ1Buffer, m_pDeviceZ2Buffer, m_pDeviceCVBuffer, m_byMaxK);

    checkCudaErrors(cudaMemcpy(m_pHostY1Buffer, m_pDeviceY1Buffer, sizeof(Real) * m_byMaxK, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostZ1Buffer, m_pDeviceZ1Buffer, sizeof(Real) * m_byMaxK, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostY2Buffer, m_pDeviceY2Buffer, sizeof(Real) * 2 * m_byMaxK, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostZ2Buffer, m_pDeviceZ2Buffer, sizeof(Real) * 2 * m_byMaxK, cudaMemcpyDeviceToHost));
}

/**
* the fry1 frz1 can be combined
*/
QLGate QLQuantumKmeans::ApplyInitialState() const
{
    QLGate gate;
    gate.AddQubits(static_cast<BYTE>(m_byQubit));
    QLGate h(EBasicOperation::EBO_H);
    for (BYTE i = 2; i < m_byQubit; ++i)
    {
        TArray<BYTE> hadmard;
        hadmard.AddItem(i);
        gate.AppendGate(h, hadmard);
    }

    TArray<BYTE> controller1;
    TArray<BYTE> controller2;
    for (BYTE i = 0; i < m_byQubit; ++i)
    {
        if (0 != i)
        {
            controller1.AddItem(m_byQubit - i);
        }
        controller2.AddItem(m_byQubit - i - 1);
    }

    QLGate fryz1 = FRyz(m_pHostY1Buffer, m_pHostZ1Buffer, m_byQubit - 1);
    QLGate fryz2 = FRyz(m_pHostY2Buffer, m_pHostZ2Buffer, m_byQubit);

    gate.AppendGate(fryz1, controller1);
    gate.AppendGate(fryz2, controller2);

    return gate;
}

void QLQuantumKmeans::ApplyInverseVector(QLGate& gate) const
{
    TArray<BYTE> q1;
    TArray<BYTE> q2;
    q1.AddItem(1);
    q2.AddItem(1);
    q2.AddItem(0);

    QLGate ry(EBasicOperation::EBO_RY, -m_pHostY1Buffer[m_byMaxK]);
    QLGate rz(EBasicOperation::EBO_RZ, m_pHostZ1Buffer[m_byMaxK]);
    Real ryz_y[2] = { m_pHostY2Buffer[2 * m_byMaxK], m_pHostY2Buffer[2 * m_byMaxK + 1]};
    Real ryz_z[2] = { m_pHostZ2Buffer[2 * m_byMaxK], m_pHostZ2Buffer[2 * m_byMaxK + 1]};
    QLGate fryz = FRyz(ryz_y, ryz_z, 2);
    fryz.Dagger();

    gate.AppendGate(rz, q1);
    gate.AppendGate(fryz, q2);
    gate.AppendGate(ry, q1);

    QLGate ctrCollapse(EBasicOperation::EBO_CC, 0);
    TArray<BYTE> q3;
    q3.AddItem(0);

    gate.AppendGate(ctrCollapse, q3);
    gate.AppendGate(ctrCollapse, q1);
}

/**
* hostVectorBuffer is vectors v[i] + w, for example, 64 v[i] + 1 w. 
* in this case, vectorCount = 65, qubitCount = 2 + log2(64) = 8
*/
QLGate QLQuantumKmeans::CompareCircuit()
{
    CalculateDegrees();
    QLGate gate = ApplyInitialState();
    ApplyInverseVector(gate);

    return gate;
}

QLGate QLQuantumKmeans::CompareCircuit(UINT uiIdx)
{
    //only once
    //CalculateDegreesOnlyCenters();

    m_pHostY1Buffer[m_byMaxK] = m_pHostDataY1Buffer[uiIdx];
    m_pHostY2Buffer[2 * m_byMaxK] = m_pHostDataY2Buffer[2 * uiIdx];
    m_pHostY2Buffer[2 * m_byMaxK + 1] = m_pHostDataY2Buffer[2 * uiIdx + 1];

    m_pHostZ1Buffer[m_byMaxK] = m_pHostDataZ1Buffer[uiIdx];
    m_pHostZ2Buffer[2 * m_byMaxK] = m_pHostDataZ2Buffer[2 * uiIdx];
    m_pHostZ2Buffer[2 * m_byMaxK + 1] = m_pHostDataZ2Buffer[2 * uiIdx + 1];

    QLGate gate = ApplyInitialState();
    ApplyInverseVector(gate);

    return gate;
}

BYTE QLQuantumKmeans::Measure(const QLGate& gate, UINT uiRepeat, UINT* count, UINT* measureCount) const
{
    QLSimulatorParametersMeasure param;
    param.m_byQubitCount = static_cast<BYTE>(gate.m_lstQubits.Num());
    param.m_MasterGate = gate;
    param.m_iRepeat = 1;
    if (uiRepeat > 1)
    {
        param.m_iMeasureUntil = static_cast<INT>(uiRepeat);
    }
    param.m_bPrint = FALSE;
    for (BYTE p = 0; p < static_cast<UINT>(gate.m_lstQubits.Num() - 2); ++p)
    {
        param.m_lstMeasureBits.AddItem(2 + p);
    }
    QLSimulatorOutputMeasure res;

    QLSimulatorMeasure sim;
    sim.Simulate(&param, &res);

    if (1 == uiRepeat)
    {
        if (NULL != count)
        {
            memcpy(count, res.m_lstCounts.GetData(), sizeof(UINT) * m_byMaxK);
        }
        if (NULL != measureCount)
        {
            measureCount[0] = 1;
        }
        return static_cast<BYTE>(res.m_lstHist[0]);
    }

    //BYTE byMaxK = 0;
    //UINT uiMax = 0;
    //for (INT i = 0; i < res.m_lstCounts.Num(); ++i)
    //{
    //    if (uiMax < res.m_lstCounts[i])
    //    {
    //        uiMax = res.m_lstCounts[i];
    //        byMaxK = static_cast<BYTE>(i);
    //    }
    //}


    if (NULL != count)
    {
        memcpy(count, res.m_lstCounts.GetData(), sizeof(UINT) * m_byMaxK);
    }

    if (NULL != measureCount)
    {
        measureCount[0] = res.m_uiMeasureUntilCount;
    }

    //return byMaxK;
    return static_cast<BYTE>(res.m_uiMeasureUntilRes);
}

UBOOL QLQuantumKmeans::CalculateCenters(UBOOL bDebug)
{
    UINT iThreadNeeded = Ceil(m_uiN, 2);
    UINT iBlock1 = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, _QL_LAUNCH_MAX_THREAD) : 1;
    UINT iThread1 = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded, iBlock1) : iThreadNeeded;

    checkCudaErrors(cudaMemcpy(m_pDeviceKValues, m_pHostKValues, sizeof(BYTE) * m_uiN, cudaMemcpyHostToDevice));
    for (BYTE byK = 0; byK < m_byMaxK; ++byK)
    {
        _kernelCenterListStep1 << <iBlock1, iThread1 >> > (
            m_pDeviceData,
            m_pDeviceRealWorkingBuffer,
            m_pDeviceUIntWorkingBuffer,
            m_uiN,
            m_pDeviceKValues, byK);
        //checkCudaErrors(cudaDeviceSynchronize());

        const UINT iRequiredDim = (iThreadNeeded + 1) >> 1;
        const UINT iPower = GetReduceDim(iRequiredDim);
        for (UINT i = 0; i <= iPower; ++i)
        {
            UINT iJump = 1 << i;
            const UINT iThreadNeeded2 = 1 << (iPower - i);
            UINT iBlock2 = iThreadNeeded2 > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded2, _QL_LAUNCH_MAX_THREAD) : 1;
            UINT iThread2 = iThreadNeeded2 > _QL_LAUNCH_MAX_THREAD ? Ceil(iThreadNeeded2, iBlock2) : iThreadNeeded2;
            _kernelCenterListStep2 << <iBlock2, iThread2 >> > (m_pDeviceRealWorkingBuffer, m_pDeviceUIntWorkingBuffer, iJump, iThreadNeeded);
            //checkCudaErrors(cudaDeviceSynchronize());
        }

        _kernelCenterListStep3 << <1, 1 >> > (m_pDeviceVBuffer + byK * 7, m_pDeviceKCounts + byK, m_pDeviceRealWorkingBuffer, m_pDeviceUIntWorkingBuffer);
        //checkCudaErrors(cudaDeviceSynchronize());
    }

    //=============================================================================================================
    //check split
    checkCudaErrors(cudaMemcpy(m_pHostKCounts, m_pDeviceKCounts, sizeof(UINT) * m_byMaxK, cudaMemcpyDeviceToHost));
    TArray<BYTE> needToSplit;
    BYTE maxCountK = 0;
    UINT maxCount = 0;
    for (BYTE byK = 0; byK < m_byMaxK; ++byK)
    {
        if (0 == m_pHostKCounts[byK])
        {
            needToSplit.AddItem(byK);
        }
        else if (m_pHostKCounts[byK] > maxCount)
        {
            maxCountK = byK;
            maxCount = m_pHostKCounts[byK];
        }
    }
    if (needToSplit.Num() > 0)
    {
        if (bDebug)
        {
            appGeneral(_T("%d empty k-cluster need to split.\n"), needToSplit.Num());
        }
        needToSplit.AddItem(maxCountK);
        checkCudaErrors(cudaMemcpy(m_pDeviceTempKBuffer, needToSplit.GetData(), sizeof(BYTE) * needToSplit.Num(), cudaMemcpyHostToDevice));
        _kernelQKMeansSpliteK << <m_uiBlockN, m_uiThreadN >> > (m_pDeviceKValues, m_uiN, maxCountK, m_pDeviceTempKBuffer, static_cast<BYTE>(needToSplit.Num()));
        checkCudaErrors(cudaMemcpy(m_pHostKValues, m_pDeviceKValues, sizeof(BYTE) * m_uiN, cudaMemcpyDeviceToHost));
        return TRUE;
    }
    return FALSE;
}

BYTE QLQuantumKmeans::Reclassify(UINT uiIdx, UINT* measureCount)
{
    //checkCudaErrors(cudaMemcpy(&m_pDeviceVBuffer[7 * m_byMaxK], &m_pDeviceData[7 * uiIdx], sizeof(Real) * 7, cudaMemcpyDeviceToDevice));

    QLGate gate = CompareCircuit(uiIdx);
    //QLGate gate = CompareCircuit();

    //appGeneral(_T("y1 %f - %f; y2 (%f, %f) - (%f, %f); z1 %f - %f; z2 (%f, %f) - (%f, %f) \n"),
    //    m_pHostY1Buffer[m_byMaxK],
    //    m_pHostDataY1Buffer[uiIdx],
    //    m_pHostY2Buffer[2 * m_byMaxK],
    //    m_pHostY2Buffer[2 * m_byMaxK + 1],
    //    m_pHostDataY2Buffer[2 * uiIdx],
    //    m_pHostDataY2Buffer[2 * uiIdx + 1],

    //    m_pHostZ1Buffer[m_byMaxK],
    //    m_pHostDataZ1Buffer[uiIdx],
    //    m_pHostZ2Buffer[2 * m_byMaxK],
    //    m_pHostZ2Buffer[2 * m_byMaxK + 1],
    //    m_pHostDataZ2Buffer[2 * uiIdx],
    //    m_pHostDataZ2Buffer[2 * uiIdx + 1]);

    return Measure(gate, m_uiRepeat, NULL, measureCount);
}

UINT QLQuantumKmeans::Reclassify(UBOOL bDebug)
{
    UINT uiChanged = 0;
    for (UINT i = 0; i < m_uiN; ++i)
    {
        static UINT measuredCount[1];
        BYTE byNewK = Reclassify(i, measuredCount);
        if (byNewK != m_pHostKValues[i])
        {
            ++uiChanged;
            m_pHostKValues[i] = byNewK;
        }
        m_pMeasureCounts[i] = measuredCount[0];
    }
    
    UBOOL bSplit = CalculateCenters(bDebug);
    while (bSplit)
    {
        bSplit = CalculateCenters(bDebug);
    }
    CalculateDegreesOnlyCenters();
    
    ++m_uiStep;
    return uiChanged;
}


void QLQuantumKmeans::TestCircuit(const Real* hostVectors)
{
    //========== put vectors to m_pDeviceVBuffer =============
    checkCudaErrors(cudaMemcpy(m_pDeviceVBuffer, hostVectors, sizeof(Real) * 7 * m_byVectorCount, cudaMemcpyHostToDevice));
    QLGate gate = CompareCircuit();

    QLSimulatorParametersVector param;
    param.m_byQubitCount = static_cast<BYTE>(gate.m_lstQubits.Num());
    param.m_MasterGate = gate;
    param.m_bPrint = TRUE;
    QLSimulatorOutputVector out;
    QLSimulatorVector sim;
    sim.Simulate(&param, &out);

    static UINT counts[256];
    static UINT measurecounts[1];
    BYTE mres = Measure(gate, 1000, counts, measurecounts);
    appGeneral(_T("measured: res %d, count %d \n"), mres, measurecounts);
    TArray<QLComplex> state_vector = out.m_OutputMatrix.ToVector();
    for (BYTE byK = 0; byK < m_byMaxK; ++byK)
    {
        QLComplex dot = _make_cuComplex(hostVectors[7 * byK] * hostVectors[7 * m_byMaxK], F(0.0));
        for (BYTE i = 1; i < 4; ++i)
        {
            dot = _cuCaddf(dot, _cuCmulf(
                _make_cuComplex(hostVectors[7 * byK + i], -hostVectors[7 * byK + i + 3]),
                _make_cuComplex(hostVectors[7 * m_byMaxK + i], hostVectors[7 * m_byMaxK + i + 3])));
        }
        Real fMeasured = _cuCabsf(state_vector[4 * byK]);
        Real fExpected = _cuCabsf(dot);
        appGeneral(_T("%d: amplitude %f, expected %f, divide %f, measured: %d, expected/mesure: %f, \n"), byK, fMeasured, fExpected, fExpected / fMeasured, counts[byK], fExpected/ _sqrt(counts[byK]));
    }
}

void QLQuantumKmeans::InitialK(UBOOL bDebug)
{
    //_kernelQKMeansInitialK << <m_uiBlockN, m_uiThreadN >> > (m_pDeviceKValues, m_uiN, m_byMaxK);

    //UBOOL split = CalculateCenters(bDebug);
    //while (split)
    //{
    //    split = CalculateCenters(bDebug);
    //}
    UINT uiStartIndex = 0;
    checkCudaErrors(cudaMemcpy(m_pDeviceVBuffer, m_pDeviceData + 7 * uiStartIndex, sizeof(Real) * 7 * m_byMaxK, cudaMemcpyDeviceToDevice));
    CalculateDegreesOnlyCenters();
}

void QLQuantumKmeans::InitialWithCenterFile()
{
    UINT w, h;
    TArray<Real> data = ReadCSVAR(m_sStartCenterFile, w, h);
    checkCudaErrors(cudaMemcpy(m_pDeviceVBuffer, data.GetData(), sizeof(Real) * 7 * m_byMaxK, cudaMemcpyHostToDevice));
    CalculateDegreesOnlyCenters();
}

void QLQuantumKmeans::KMeans(const CCString& sResultFileName, UINT uiStop, UINT uiRepeat, UINT uiStep, UBOOL bDebug)
{
    m_sSaveNameHead = sResultFileName;
    m_uiStep = uiStep;
    m_uiRepeat = uiRepeat;

    //========== step1: initial =============
    if (m_sStartCenterFile.IsEmpty())
    {
        InitialK(bDebug);
    }
    else
    {
        InitialWithCenterFile();
    }
    
    //========== step2: calculate cetners =============
    //========== step3: reclassify =============
    UINT changed = Reclassify(bDebug);
    ExportDebugInfo();
    if (bDebug)
    {
        appGeneral(_T("Step: %d, Changed: %d\n"), m_uiStep, changed);
    }

    //========== step4: repeat step2 and step3 until finished =============
    while (changed > uiStop)
    {
        changed = Reclassify(bDebug);
        ExportDebugInfo();
        if (bDebug)
        {
            appGeneral(_T("Step: %d, Changed: %d\n"), m_uiStep, changed);
        }
    }

    //========== step5: save k file =============
}

void QLQuantumKmeans::ExportDebugInfo()
{
    SaveCSVAB(m_pHostKValues, 1, m_uiN, m_sSaveNameHead + _T("_k_") + appIntToString(static_cast<INT>(m_uiStep)) + _T(".csv"));

    checkCudaErrors(cudaMemcpy(m_pHostCenters, m_pDeviceVBuffer, sizeof(Real) * 7 * m_byMaxK, cudaMemcpyDeviceToHost));
    SaveCSVAR(m_pHostCenters, 7, m_byMaxK, m_sSaveNameHead + _T("_c_") + appIntToString(static_cast<INT>(m_uiStep)) + _T(".csv"));

    SaveCSVAUI(m_pMeasureCounts, 1, m_uiN, m_sSaveNameHead + _T("_m_") + appIntToString(static_cast<INT>(m_uiStep)) + _T(".csv"));
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================