//=============================================================================
// FILENAME : QLRandom.h
// 
// DESCRIPTION:
// This is random number for parallel
//
//
// REVISION: [dd/mm/yy]
//  [12/09/2022 nbale]
//=============================================================================

#ifndef _QLRANDOM_H_
#define _QLRANDOM_H_


#define __SOBEL_OFFSET_MAX (4096)

#define __DefineRandomFuncion(rettype, funcname) __device__ __inline__ static rettype _deviceRandom##funcname(UINT uiThreadIdx) \
{ \
    return __r->_deviceRandom##funcname(uiThreadIdx); \
}


__BEGIN_NAMESPACE

enum class ERandom : UINT
{
    ER_Schrage,

    ER_XORWOW,
    ER_MRG32K3A,
    //ER_MTGP32, //see the document, this has lots of constraints.
    ER_PHILOX4_32_10,
    ER_QUASI_SOBOL32,
    ER_SCRAMBLED_SOBOL32,
};

enum class ERandomSeedType : UINT
{
    ERST_Number,
    ERST_Timestamp,
};

class QLAPI QLRandom
{
public:

    /**
    * There are nine types of random number generators in cuRAND, that fall into two categories.
    *  CURAND_RNG_PSEUDO_XORWOW, CURAND_RNG_PSEUDO_MRG32K3A, CURAND_RNG_PSEUDO_MTGP32, CURAND_RNG_PSEUDO_PHILOX4_32_10 and CURAND_RNG_PSEUDO_MT19937 are pseudorandom number generators.
    * CURAND_RNG_PSEUDO_XORWOW is implemented using the XORWOW algorithm, a member of the xor-shift family of pseudorandom number generators.
    * CURAND_RNG_PSEUDO_MRG32K3A is a member of the Combined Multiple Recursive family of pseudorandom number generators.
    *  CURAND_RNG_PSEUDO_MT19937 and CURAND_RNG_PSEUDO_MTGP32 are members of the Mersenne Twister family of pseudorandom number generators.
    * CURAND_RNG_PSEUDO_MTGP32 has parameters customized for operation on the GPU.
    * CURAND_RNG_PSEUDO_MT19937 has the same parameters as CPU version, but ordering is different.
    * CURNAD_RNG_PSEUDO_MT19937 supports only HOST API and can be used only on architecture sm_35 or higher.
    * CURAND_RNG_PHILOX4_32_10 is a member of Philox family, which is one of the three non-cryptographic Counter Based Random Number Generators presented on SC11 conference by D E Shaw Research.
    *
    *  There are 4 variants of the basic SOBOL¡¯ quasi random number generator. All of the variants generate sequences in up to 20,000 dimensions. CURAND_RNG_QUASI_SOBOL32, CURAND_RNG_QUASI_SCRAMBLED_SOBOL32, CURAND_RNG_QUASI_SOBOL64, and CURAND_RNG_QUASI_SCRAMBLED_SOBOL64 are quasirandom number generator types.
    * CURAND_RNG_QUASI_SOBOL32 is a Sobol¡¯ generator of 32-bit sequences.
    * CURAND_RNG_QUASI_SCRAMBLED_SOBOL32 is a scrambled Sobol¡¯ generator of 32-bit sequences.
    * CURAND_RNG_QUASI_SOBOL64 is a Sobol¡¯ generator of 64-bit sequences.
    * CURAND_RNG_QUASI_SCRAMBLED_SOBOL64 is a scrambled Sobol¡¯ generator of 64-bit sequences.
    */
    QLRandom(UINT uiSeed, UINT uiMaxThread, ERandom er)
        : m_eRandomType(er)
        , m_uiMaxThread(uiMaxThread)
        , m_uiHostSeed(uiSeed)
    {
        switch (er)
        {
        case ERandom::ER_Schrage:
        {
            InitialTableSchrage();
        }
        break;
        case ERandom::ER_MRG32K3A:
        {
            checkCudaErrors(curandCreateGenerator(&m_HGen, CURAND_RNG_PSEUDO_MRG32K3A));
            checkCudaErrors(curandSetPseudoRandomGeneratorSeed(m_HGen, uiSeed));
            checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
            InitialStatesMRG();
        }
        break;
        case ERandom::ER_PHILOX4_32_10:
        {
            checkCudaErrors(curandCreateGenerator(&m_HGen, CURAND_RNG_PSEUDO_PHILOX4_32_10));
            checkCudaErrors(curandSetPseudoRandomGeneratorSeed(m_HGen, uiSeed));
            checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
            InitialStatesPhilox();
        }
        break;
        case ERandom::ER_QUASI_SOBOL32:
        {
            //for sobol, on the host, we use XORWOW
            checkCudaErrors(curandCreateGenerator(&m_HGen, CURAND_RNG_QUASI_SOBOL32));
            checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
            InitialStatesSobol32();
        }
        break;
        case ERandom::ER_SCRAMBLED_SOBOL32:
        {
            //for sobol, on the host, we use XORWOW
            checkCudaErrors(curandCreateGenerator(&m_HGen, CURAND_RNG_QUASI_SCRAMBLED_SOBOL32));
            checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
            InitialStatesScrambledSobol32();
        }
        break;
        case ERandom::ER_XORWOW:
        default:
        {
            checkCudaErrors(curandCreateGenerator(&m_HGen, CURAND_RNG_PSEUDO_XORWOW));
            checkCudaErrors(curandSetPseudoRandomGeneratorSeed(m_HGen, uiSeed));
            checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
            InitialStatesXORWOW();
        }
        break;
        }

        checkCudaErrors(cudaGetLastError());
    }

    ~QLRandom();

    /**
    * Note that this gives [0, 1), and curand_uniform gives (0, 1]
    */
    __device__ __inline__ FLOAT _deviceRandomF(UINT uithreadIdx) const
    {
        switch (m_eRandomType)
        {
        case ERandom::ER_Schrage:
            return static_cast<FLOAT>(AM * _deviceRandomUISchrage(uithreadIdx));
        case ERandom::ER_MRG32K3A:
            return 1.0f - curand_uniform(&(m_pDeviceRandStatesMRG[uithreadIdx]));
        case ERandom::ER_PHILOX4_32_10:
            return 1.0f - curand_uniform(&(m_pDeviceRandStatesPhilox[uithreadIdx]));
        case ERandom::ER_QUASI_SOBOL32:
            return 1.0f - curand_uniform(&(m_pDeviceRandStatesSobol32[(uithreadIdx)]));
        case ERandom::ER_SCRAMBLED_SOBOL32:
            return 1.0f - curand_uniform(&(m_pDeviceRandStatesScrambledSobol32[uithreadIdx]));
        case ERandom::ER_XORWOW:
        default:
            return 1.0f - curand_uniform(&(m_pDeviceRandStatesXORWOW[uithreadIdx]));
        }

        //return 0;
    }

    __device__ __inline__ QLComplex _deviceRandomC(UINT uithreadIdx) const
    {
        const Real f1 = _deviceRandomF(uithreadIdx) * 2.0f - 1.0f;
        const Real f2 = _deviceRandomF(uithreadIdx) * 2.0f - 1.0f;
        return _make_cuComplex(f1, f2);
    }

    /**
    * The standard deviation of it is 1
    */
    __device__ __inline__ Real _deviceRandomGaussF(UINT uithreadIdx) const
    {
        const Real f1 = _deviceRandomF(uithreadIdx);
        const Real f2 = _deviceRandomF(uithreadIdx) * PI2;

        const Real oneMinusf1 = F(1.0) - f1;
        const Real inSqrt = -F(2.0) * _log(oneMinusf1 > F(0.0) ? oneMinusf1 : _QL_FLT_MIN_);
        const Real amplitude = (inSqrt > F(0.0) ? _sqrt(inSqrt) : F(0.0)) * InvSqrt2;
        return _cos(f2) * amplitude;
    }

    __device__ __inline__ QLComplex _deviceRandomGaussC(UINT uithreadIdx) const
    {
        const Real f1 = _deviceRandomF(uithreadIdx);
        const Real f2 = _deviceRandomF(uithreadIdx) * PI;
        const Real oneMinusf1 = F(1.0) - f1;
        const Real inSqrt = -F(2.0) * _log(oneMinusf1 > F(0.0) ? oneMinusf1 : _QL_FLT_MIN_);
        const Real amplitude = (inSqrt > F(0.0) ? _sqrt(inSqrt) : F(0.0)) * InvSqrt2;
        return _make_cuComplex(_cos(f2) * amplitude, _sin(f2) * amplitude);
    }

    __device__ __inline__ UINT _deviceRandomI(UINT uithreadIdx, UINT uiMax) const
    {
        UINT toRet = static_cast<UINT>(uiMax * _deviceRandomF(uithreadIdx));
        if (toRet >= uiMax)
        {
            toRet = uiMax - 1;
        }
        return toRet;
    }

    __device__ __inline__ Real _deviceRandomIF(UINT uithreadIdx, UINT uiMax) const
    {
        return static_cast<Real>(static_cast<UINT>(uiMax * _deviceRandomF(uithreadIdx)));
    }

    __device__ __inline__ QLComplex _deviceRandomZN(UINT uithreadIdx, UINT uiMax) const
    {
        const Real byRandom = static_cast<Real>(_deviceRandomI(uithreadIdx, uiMax));
        const Real arg = PI2 * byRandom / static_cast<Real>(uiMax);
        return _make_cuComplex(_cos(arg), _sin(arg));
    }

    __host__ __inline__ Real GetRandomF()
    {
        if (ERandom::ER_Schrage == m_eRandomType)
        {
            return static_cast<Real>(AM * GetRandomUISchrage());
        }

        curandGenerateUniform(m_HGen, m_deviceBuffer, 1);
        checkCudaErrors(cudaMemcpy(m_hostBuffer, m_deviceBuffer, sizeof(FLOAT), cudaMemcpyDeviceToHost));
        return static_cast<Real>(m_hostBuffer[0]);
    }

    FLOAT* m_deviceBuffer;
    Real m_hostBuffer[1];
    curandGenerator_t m_HGen;
    ERandom m_eRandomType;
    UINT m_uiMaxThread;

protected:

    void InitialStatesXORWOW();
    void InitialStatesPhilox();
    void InitialStatesMRG();
    void InitialStatesSobol32();
    void InitialStatesScrambledSobol32();

    curandState* m_pDeviceRandStatesXORWOW;
    curandStatePhilox4_32_10_t* m_pDeviceRandStatesPhilox;
    curandStateMRG32k3a* m_pDeviceRandStatesMRG;

    curandStateSobol32* m_pDeviceRandStatesSobol32;
    curandDirectionVectors32_t* m_pDeviceSobolDirVec;
    UINT* m_pDeviceSobelConsts;
    curandStateScrambledSobol32* m_pDeviceRandStatesScrambledSobol32;

#pragma region Schrage

public:

    UINT* m_pDeviceSeedTable;

    /**
    * run on device, parally set the table
    */
    __device__ __inline__ static void _deviceAsignSeeds(UINT* devicePtr, UINT uiSeed, UINT uiFatIndex)
    {
        devicePtr[uiFatIndex] = (1664525UL * (uiFatIndex + uiSeed) + 1013904223UL) & 0xffffffff;
    }

protected:

    void InitialTableSchrage();

    __device__ __inline__ UINT _deviceRandomUISchrage(UINT threadIndex) const
    {
        m_pDeviceSeedTable[threadIndex] = ((1664525UL * m_pDeviceSeedTable[threadIndex] + 1013904223UL) & 0xffffffff);
        return m_pDeviceSeedTable[threadIndex];
    }

    __host__ __inline__ UINT GetRandomUISchrage()
    {
        m_uiHostSeed = (1664525UL * m_uiHostSeed + 1013904223UL) & 0xffffffff;
        return m_uiHostSeed;
    }

    UINT m_uiHostSeed;

#pragma endregion

};

extern __constant__ class QLRandom* __r;

__DefineRandomFuncion(Real, F)

__DefineRandomFuncion(Real, GaussF)

__DefineRandomFuncion(QLComplex, C)

__DefineRandomFuncion(QLComplex, GaussC)

__device__ __inline__ static Real _deviceRandomIF(UINT uiThreadIdx, UINT uiN)
{
    return __r->_deviceRandomIF(uiThreadIdx, uiN);
}

__device__ __inline__ static QLComplex _deviceRandomZN(UINT uiThreadIdx, UINT uiN)
{
    return __r->_deviceRandomZN(uiThreadIdx, uiN);
}

class QLAPI QLRandomInitializer
{
public:
    QLRandomInitializer(ERandom eRandom = ERandom::ER_Schrage, UINT uiSeed = 0);
    ~QLRandomInitializer();

    QLRandom* m_pRandom;
    QLRandom* m_pDeviceRandom;
};

//extern QLRandomInitializer GRandom;

__END_NAMESPACE

#endif //#ifndef _QLRANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================
