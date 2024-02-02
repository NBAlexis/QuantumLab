//=============================================================================
// FILENAME : QLRandom.cpp
// 
// DESCRIPTION:
//
//
// REVISION: [dd/mm/yy]
//  [12/09/2022 nbale]
//=============================================================================
#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

__global__ void _QL_LAUNCH_BOUND
_kernalAllocateSeedTable(UINT* pDevicePtr, UINT uiSeed)
{
    UINT uiThread = threadIdx.x;
    QLRandom::_deviceAsignSeeds(pDevicePtr, uiSeed + uiThread, uiThread);
}

__global__ void _QL_LAUNCH_BOUND
_kernalInitialXORWOW(curandState* states, UINT uiSeed)
{
    UINT uiThread = threadIdx.x;
    curand_init(uiSeed, uiThread, 0, &states[uiThread]);
}

__global__ void _QL_LAUNCH_BOUND
_kernalInitialPhilox(curandStatePhilox4_32_10_t* states, UINT uiSeed)
{
    UINT uiThread = threadIdx.x;
    curand_init(uiSeed, uiThread, 0, &states[uiThread]);
}

__global__ void _QL_LAUNCH_BOUND
_kernalInitialMRG(curandStateMRG32k3a* states, UINT uiSeed)
{
    UINT uiThread = threadIdx.x;
    curand_init(uiSeed, uiThread, 0, &states[uiThread]);
}

__global__ void _QL_LAUNCH_BOUND
_kernalInitialSobel32(curandStateSobol32* states, curandDirectionVectors32_t* dirs, UINT uiSeed)
{
    UINT uiThread = threadIdx.x;
    curand_init(dirs[uiThread], uiSeed % 16, &states[uiThread]);
}

__global__ void _QL_LAUNCH_BOUND
_kernalInitialScrambledSobel32(curandStateScrambledSobol32* states, UINT* consts, curandDirectionVectors32_t* dirs, UINT uiSeed)
{
    UINT uiThread = threadIdx.x;
    curand_init(dirs[uiThread], consts[uiThread], uiSeed % __SOBEL_OFFSET_MAX, &states[uiThread]);
}

QLRandom::~QLRandom()
{

    switch (m_eRandomType)
    {
    case ERandom::ER_Schrage:
    {
        checkCudaErrors(cudaFree(m_pDeviceSeedTable));
    }
    break;
    case ERandom::ER_MRG32K3A:
    {
        checkCudaErrors(curandDestroyGenerator(m_HGen));
        checkCudaErrors(cudaFree(m_deviceBuffer));
        checkCudaErrors(cudaFree(m_pDeviceRandStatesMRG));
    }
    break;
    case ERandom::ER_PHILOX4_32_10:
    {
        checkCudaErrors(curandDestroyGenerator(m_HGen));
        checkCudaErrors(cudaFree(m_deviceBuffer));
        checkCudaErrors(cudaFree(m_pDeviceRandStatesPhilox));
    }
    break;
    case ERandom::ER_QUASI_SOBOL32:
    {
        checkCudaErrors(curandDestroyGenerator(m_HGen));
        checkCudaErrors(cudaFree(m_deviceBuffer));
        checkCudaErrors(cudaFree(m_pDeviceRandStatesSobol32));
        checkCudaErrors(cudaFree(m_pDeviceSobolDirVec));
    }
    break;
    case ERandom::ER_SCRAMBLED_SOBOL32:
    {
        checkCudaErrors(curandDestroyGenerator(m_HGen));
        checkCudaErrors(cudaFree(m_deviceBuffer));
        checkCudaErrors(cudaFree(m_pDeviceRandStatesScrambledSobol32));
        checkCudaErrors(cudaFree(m_pDeviceSobolDirVec));
        checkCudaErrors(cudaFree(m_pDeviceSobelConsts));
    }
    break;
    case ERandom::ER_XORWOW:
    default:
    {
        checkCudaErrors(curandDestroyGenerator(m_HGen));
        checkCudaErrors(cudaFree(m_deviceBuffer));
        checkCudaErrors(cudaFree(m_pDeviceRandStatesXORWOW));
    }
    break;
    }
}

//Initial XORWOW only support 512 threads per block
void QLRandom::InitialStatesXORWOW()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRandStatesXORWOW, sizeof(curandState) * m_uiMaxThread));
    _kernalInitialXORWOW << <1, m_uiMaxThread >> > (m_pDeviceRandStatesXORWOW, m_uiHostSeed);
}

//Initial Philox only support 256 threads per block
void QLRandom::InitialStatesPhilox()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRandStatesPhilox, sizeof(curandStatePhilox4_32_10_t) * m_uiMaxThread));
    _kernalInitialPhilox << <1, m_uiMaxThread >> > (m_pDeviceRandStatesPhilox, m_uiHostSeed);
}

//Initial MRG only support 256 threads per block
void QLRandom::InitialStatesMRG()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRandStatesMRG, sizeof(curandStateMRG32k3a) * m_uiMaxThread));
    _kernalInitialMRG << <1, m_uiMaxThread >> > (m_pDeviceRandStatesMRG, m_uiHostSeed);
}

void QLRandom::InitialStatesSobol32()
{
    //support only 20000 dimensions, so using _HC_Volumn instead
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRandStatesSobol32, sizeof(curandStateSobol32) * m_uiMaxThread));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceSobolDirVec, sizeof(curandDirectionVectors32_t) * m_uiMaxThread));

    //int[32]
    curandDirectionVectors32_t* hostVectors32;
    checkCudaErrors(curandGetDirectionVectors32(&hostVectors32, CURAND_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(m_pDeviceSobolDirVec, hostVectors32, sizeof(curandDirectionVectors32_t) * m_uiMaxThread, cudaMemcpyHostToDevice));

    _kernalInitialSobel32 << <1, m_uiMaxThread >> > (m_pDeviceRandStatesSobol32, m_pDeviceSobolDirVec, m_uiHostSeed);
}

void QLRandom::InitialStatesScrambledSobol32()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRandStatesScrambledSobol32, sizeof(curandStateScrambledSobol32) * m_uiMaxThread));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceSobolDirVec, sizeof(curandDirectionVectors32_t) * m_uiMaxThread));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceSobelConsts, sizeof(UINT) * m_uiMaxThread));

    curandDirectionVectors32_t* hostVectors32;
    checkCudaErrors(curandGetDirectionVectors32(&hostVectors32, CURAND_SCRAMBLED_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(m_pDeviceSobolDirVec, hostVectors32, sizeof(curandDirectionVectors32_t) * m_uiMaxThread, cudaMemcpyHostToDevice));

    UINT* hostScrambleConstants32;
    checkCudaErrors(curandGetScrambleConstants32(&hostScrambleConstants32));
    checkCudaErrors(cudaMemcpy(m_pDeviceSobelConsts, hostScrambleConstants32, sizeof(UINT) * m_uiMaxThread, cudaMemcpyHostToDevice));

    _kernalInitialScrambledSobel32 << <1, m_uiMaxThread >> > (m_pDeviceRandStatesScrambledSobol32, m_pDeviceSobelConsts, m_pDeviceSobolDirVec, m_uiHostSeed);
}

void QLRandom::InitialTableSchrage()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceSeedTable, sizeof(UINT) * m_uiMaxThread));
    _kernalAllocateSeedTable << <1, m_uiMaxThread >> > (m_pDeviceSeedTable, m_uiMaxThread);
}

//to do add options to initialize the random
QLRandomInitializer::QLRandomInitializer(ERandom eRandom, UINT uiSeed)
    : m_pRandom(NULL)
    , m_pDeviceRandom(NULL)
{
    if (0 == uiSeed)
    {
        uiSeed = appGetTimeStamp();
    }

    m_pRandom = new QLRandom(uiSeed, _QL_LAUNCH_MAX_THREAD, eRandom);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceRandom, sizeof(QLRandom)));
    checkCudaErrors(cudaMemcpy(m_pDeviceRandom, m_pRandom, sizeof(QLRandom), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpyToSymbol(__r, &m_pDeviceRandom, sizeof(QLRandom*)));

    printf("random initialed\n");

    _hostRandomPointer = m_pRandom;
}

QLRandomInitializer::~QLRandomInitializer()
{
    checkCudaErrors(cudaFree(m_pDeviceRandom));
    appSafeDelete(m_pRandom);
}

__constant__ QLRandom* __r;
//QLRandomInitializer GRandom;
QLRandom* _hostRandomPointer = NULL;

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
