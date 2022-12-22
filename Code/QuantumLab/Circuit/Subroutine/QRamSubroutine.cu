//=============================================================================
// FILENAME : QRamSubroutine.cu
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


QLQuantumKmeans::QLQuantumKmeans(BYTE maxK, const Real* pData)
    : m_byMaxK(maxK)
    , m_byQubit(0)
    , m_byVectorCount(maxK + 1)
    , m_uiBlock(0)
    , m_uiThread(0)
    , m_pDeviceVBuffer(NULL)
    , m_pDeviceCVBuffer(NULL)
    , m_pDeviceY1Buffer(NULL)
    , m_pDeviceY2Buffer(NULL)
    , m_pDeviceZ1Buffer(NULL)
    , m_pDeviceZ2Buffer(NULL)
    , m_pHostY1Buffer(NULL)
    , m_pHostY2Buffer(NULL)
    , m_pHostZ1Buffer(NULL)
    , m_pHostZ2Buffer(NULL)
{
    m_byQubit = 2 + MostSignificantPowerTwo(maxK);

    m_uiBlock = m_byVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(m_byVectorCount, _QL_LAUNCH_MAX_THREAD) : 1;
    m_uiThread = m_byVectorCount > _QL_LAUNCH_MAX_THREAD ? Ceil(m_byVectorCount, m_uiBlock) : m_byVectorCount;

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceVBuffer, sizeof(Real) * 7 * m_byVectorCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceCVBuffer, sizeof(QLComplex) * 4 * m_byVectorCount));
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

/**
* host vector buffer is vectors v[i] + w, for example, 64 v[i] + 1 w.
* in this case, vectorCount = 65, qubitCount = 2 + log2(64) = 8
* dimension of host buffer is 7, and are expected to be already initialized
*/
void QLQuantumKmeans::CalculateDegrees(const Real* hostVectorBuffer)
{
    checkCudaErrors(cudaMemcpy(m_pDeviceVBuffer, hostVectorBuffer, sizeof(Real) * 7 * m_byVectorCount, cudaMemcpyHostToDevice));

    _kernelQKMADVectorToNormalizedVector << <m_uiBlock, m_uiThread >> > (m_pDeviceCVBuffer, m_pDeviceVBuffer, m_byVectorCount);
    _kernelQKMADVectorToAngle << <m_uiBlock, m_uiThread >> > (m_pDeviceY1Buffer, m_pDeviceY2Buffer, m_pDeviceZ1Buffer, m_pDeviceZ2Buffer, m_pDeviceCVBuffer, m_byVectorCount);

    checkCudaErrors(cudaMemcpy(m_pHostY1Buffer, m_pDeviceY1Buffer, sizeof(Real) * m_byVectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostZ1Buffer, m_pDeviceZ1Buffer, sizeof(Real) * m_byVectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostY2Buffer, m_pDeviceY2Buffer, sizeof(Real) * 2 * m_byVectorCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostZ2Buffer, m_pDeviceZ2Buffer, sizeof(Real) * 2 * m_byVectorCount, cudaMemcpyDeviceToHost));
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
}

/**
* hostVectorBuffer is vectors v[i] + w, for example, 64 v[i] + 1 w. 
* in this case, vectorCount = 65, qubitCount = 2 + log2(64) = 8
*/
QLGate QLQuantumKmeans::CompareCircuit(const Real* hostVectorBuffer)
{
    CalculateDegrees(hostVectorBuffer);
    QLGate gate = ApplyInitialState();
    ApplyInverseVector(gate);

    return gate;
}



#pragma endregion

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================