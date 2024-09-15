//=============================================================================
// FILENAME : SimpleEncode.cu
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [06/05/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* 
*/
__global__ void _QL_LAUNCH_BOUND
_kernelSE_FetchDegrees(Real* YLst, Real* ZLst, const QLComplex * __restrict__ vBuffer, UINT vectorCount, UINT lengthCount)
{
    const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);

    if (idx < vectorCount * lengthCount)
    {
        UINT uiV = idx / lengthCount;
        UINT uiD = idx % lengthCount;

        YLst[uiD * vectorCount + uiV] = vBuffer[uiV * lengthCount + uiD].x * PI;
        ZLst[uiD * vectorCount + uiV] = -vBuffer[uiV * lengthCount + uiD].y * PI;
    }
}

#pragma endregion

/**
* ry rz CNOT
*/
QLGate QLAPI SimpleEncodeOneVector(const QLComplex* hostv, BYTE qubits, UINT uiVLength)
{
    QLGate ret;
    ret.AddQubits(qubits);

    QLGate cnot(EBasicOperation::EBO_CX);

    BYTE toaddbyte = 0;

    for (UINT i = 0; i < uiVLength; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, hostv[i].x * PI);
        QLGate rz(EBasicOperation::EBO_RZ, hostv[i].y * PI);

        ret.AppendGate(ry, toaddbyte);
        ret.AppendGate(rz, toaddbyte);

        toaddbyte = toaddbyte + 1;
        if (toaddbyte == qubits && i != (uiVLength - 1))
        {
            for (BYTE j = 0; j < qubits; ++j)
            {
                if (0 == j)
                {
                    ret.AppendGate(cnot, j, qubits - 1);
                }
                else
                {
                    ret.AppendGate(cnot, j, j - 1);
                }
            }
            toaddbyte = 0;
        }
    }

    return ret;
}

QLGate QLAPI SimpleEncodeVectors(const QLComplex* hostv, BYTE vectorCountPower, BYTE qubits, UINT uiVLength)
{
    QLGate ret;
    //const BYTE uiVPower = static_cast<BYTE>(MostSignificantPowerTwo(uiVLength));
    ret.AddQubits(vectorCountPower + qubits);
    const UINT uiVectorCount = static_cast<UINT>(1UL << vectorCountPower);
    const UINT uiNumberCount = uiVectorCount * uiVLength;

    QLComplex* devicev = NULL;
    Real* devicedegreeY = NULL;
    Real* devicedegreeZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&devicev, sizeof(QLComplex) * uiNumberCount));
    checkCudaErrors(cudaMalloc((void**)&devicedegreeY, sizeof(Real) * uiNumberCount));
    checkCudaErrors(cudaMalloc((void**)&devicedegreeZ, sizeof(Real) * uiNumberCount));
    checkCudaErrors(cudaMemcpy(devicev, hostv, sizeof(QLComplex) * uiNumberCount, cudaMemcpyHostToDevice));
    UINT uiBlock = Ceil(uiNumberCount, _QL_LAUNCH_MAX_THREAD);
    _kernelSE_FetchDegrees << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (devicedegreeY, devicedegreeZ, devicev, uiVectorCount, uiVLength);

    Real* hostdegreeY = (Real*)malloc(sizeof(Real) * uiNumberCount);
    Real* hostdegreeZ = (Real*)malloc(sizeof(Real) * uiNumberCount);
    checkCudaErrors(cudaMemcpy(hostdegreeY, devicedegreeY, sizeof(Real) * uiNumberCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(hostdegreeZ, devicedegreeZ, sizeof(Real) * uiNumberCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(devicev));
    checkCudaErrors(cudaFree(devicedegreeY));
    checkCudaErrors(cudaFree(devicedegreeZ));

    //QLGate fryz = FRyz(hostY + idx, hostZ + idx, static_cast<UINT>(bits.Num()));
    QLGate cnot(EBasicOperation::EBO_CX);
    QLGate h(EBasicOperation::EBO_H);

    BYTE toaddbyte = 0;
    //UBOOL bLastLevelCnotAdded = FALSE;
    TArray<BYTE> bits;
    for (BYTE b = 0; b < vectorCountPower; ++b)
    {
        bits.AddItem(vectorCountPower + qubits - 1 - b);
        ret.AppendGate(h, vectorCountPower + qubits - 1 - b);
    }
    bits.AddItem(0);

    for (UINT i = 0; i < uiVLength; ++i)
    {
        appParanoiac(_T("encode feature %d\n"), i);
        QLGate fryz = FRyz(hostdegreeY + i * uiVectorCount, hostdegreeZ + i * uiVectorCount, static_cast<UINT>(vectorCountPower + 1));
        //QLGate fry = FRy(hostdegreeY + i * uiVectorCount, static_cast<UINT>(vectorCountPower + 1));

        //TArray<Real> degrees;
        //for (UINT k = 0; k < uiVectorCount; ++k)
        //{
        //    degrees.AddItem(hostdegreeY[i * uiVectorCount + k]);
        //    degrees.AddItem(hostdegreeZ[i * uiVectorCount + k]);
        //}
        //appGeneral(_T("degress: %s\n"), appToString(degrees));

        //For example, if we have 16 vectors (length-8 complex), we need
        //6543-0, 6543-1, 6543-2, ...
        
        bits[vectorCountPower] = toaddbyte;
        ret.AppendGate(fryz, bits);
        //ret.AppendGate(fry, bits);

        toaddbyte = toaddbyte + 1;
        if (toaddbyte == qubits && i != (uiVLength - 1))
        {
            for (BYTE j = 0; j < qubits; ++j)
            {
                if (0 == j)
                {
                    ret.AppendGate(cnot, j, qubits - 1);
                }
                else
                {
                    ret.AppendGate(cnot, j, j - 1);
                }
            }
            toaddbyte = 0;
            //bLastLevelCnotAdded = TRUE;
        }
        //else
        //{
        //    bLastLevelCnotAdded = FALSE;
        //}
    }

    return ret;
}

QLGate QLAPI SimpleEncodeOneVectorWithLinkStype(const QLComplex* hostv, ELinkStyle eStype, ELinkLayer eLayer, BYTE qubits, UINT uiVLength)
{
    QLGate ret;
    ret.AddQubits(qubits);

    BYTE toaddbyte = 0;
    UINT level = 0;
    for (UINT i = 0; i < uiVLength; ++i)
    {
        QLGate ry(EBasicOperation::EBO_RY, hostv[i].x * PI);
        QLGate rz(EBasicOperation::EBO_RZ, hostv[i].y * PI);

        ret.AppendGate(ry, toaddbyte);
        ret.AppendGate(rz, toaddbyte);

        toaddbyte = toaddbyte + 1;
        if (toaddbyte == qubits && i != (uiVLength - 1))
        {
            SimpleEncodeAddLinkLayerAll(eStype, eLayer, ret, level, ret.m_lstQubits);
            ++level;
            toaddbyte = 0;
        }
    }

    return ret;
}

QLGate QLAPI SimpleEncodeVectorsWithLinkStype(const QLComplex* hostv, ELinkStyle eStype, ELinkLayer eLayer, BYTE vectorCountPower, BYTE qubits, UINT uiVLength)
{
    QLGate ret;
    //const BYTE uiVPower = static_cast<BYTE>(MostSignificantPowerTwo(uiVLength));
    ret.AddQubits(vectorCountPower + qubits);
    const UINT uiVectorCount = static_cast<UINT>(1UL << vectorCountPower);
    const UINT uiNumberCount = uiVectorCount * uiVLength;
    TArray<BYTE> vectorbits;
    vectorbits.Append(ByteSequnce, qubits);

    QLComplex* devicev = NULL;
    Real* devicedegreeY = NULL;
    Real* devicedegreeZ = NULL;
    checkCudaErrors(cudaMalloc((void**)&devicev, sizeof(QLComplex) * uiNumberCount));
    checkCudaErrors(cudaMalloc((void**)&devicedegreeY, sizeof(Real) * uiNumberCount));
    checkCudaErrors(cudaMalloc((void**)&devicedegreeZ, sizeof(Real) * uiNumberCount));
    checkCudaErrors(cudaMemcpy(devicev, hostv, sizeof(QLComplex) * uiNumberCount, cudaMemcpyHostToDevice));
    UINT uiBlock = Ceil(uiNumberCount, _QL_LAUNCH_MAX_THREAD);
    _kernelSE_FetchDegrees << <uiBlock, _QL_LAUNCH_MAX_THREAD >> > (devicedegreeY, devicedegreeZ, devicev, uiVectorCount, uiVLength);

    Real* hostdegreeY = (Real*)malloc(sizeof(Real) * uiNumberCount);
    Real* hostdegreeZ = (Real*)malloc(sizeof(Real) * uiNumberCount);
    checkCudaErrors(cudaMemcpy(hostdegreeY, devicedegreeY, sizeof(Real) * uiNumberCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(hostdegreeZ, devicedegreeZ, sizeof(Real) * uiNumberCount, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(devicev));
    checkCudaErrors(cudaFree(devicedegreeY));
    checkCudaErrors(cudaFree(devicedegreeZ));

    //QLGate fryz = FRyz(hostY + idx, hostZ + idx, static_cast<UINT>(bits.Num()));
    QLGate h(EBasicOperation::EBO_H);

    BYTE toaddbyte = 0;
    TArray<BYTE> bits;
    for (BYTE b = 0; b < vectorCountPower; ++b)
    {
        bits.AddItem(vectorCountPower + qubits - 1 - b);
        ret.AppendGate(h, vectorCountPower + qubits - 1 - b);
    }
    bits.AddItem(0);

    UINT level = 0;
    for (UINT i = 0; i < uiVLength; ++i)
    {
        appParanoiac(_T("encode feature %d\n"), i);
        QLGate fryz = FRyz(hostdegreeY + i * uiVectorCount, hostdegreeZ + i * uiVectorCount, static_cast<UINT>(vectorCountPower + 1));
        //QLGate fry = FRy(hostdegreeY + i * uiVectorCount, static_cast<UINT>(vectorCountPower + 1));

        //TArray<Real> degrees;
        //for (UINT k = 0; k < uiVectorCount; ++k)
        //{
        //    degrees.AddItem(hostdegreeY[i * uiVectorCount + k]);
        //    degrees.AddItem(hostdegreeZ[i * uiVectorCount + k]);
        //}
        //appGeneral(_T("degress: %s\n"), appToString(degrees));

        //For example, if we have 16 vectors (length-8 complex), we need
        //6543-0, 6543-1, 6543-2, ...

        bits[vectorCountPower] = toaddbyte;
        ret.AppendGate(fryz, bits);
        //ret.AppendGate(fry, bits);

        toaddbyte = toaddbyte + 1;
        if (toaddbyte == qubits && i != (uiVLength - 1))
        {
            SimpleEncodeAddLinkLayerAll(eStype, eLayer, ret, level, vectorbits);
            ++level;
            toaddbyte = 0;
        }
    }

    return ret;
}

void QLAPI SimpleEncodeAddLinkLayerAll(ELinkStyle eStype, ELinkLayer eLinkLayer, QLGate& gate, UINT level, const TArray<BYTE>& linkQubits)
{
    if (ELinkLayer::CRX == eLinkLayer)
    {
        appWarning(_T("crx layer not supported, change to cz\n"));
    }
    UBOOL bCx = (ELinkLayer::CX == eLinkLayer);
    QLGate cz(EBasicOperation::EBO_CZ);
    QLGate cx(EBasicOperation::EBO_CX);

    switch (eStype)
    {
    case ELinkStyle::Full:
    {
        for (INT q1 = 0; q1 < linkQubits.Num(); ++q1)
        {
            for (INT q2 = q1 + 1; q2 < linkQubits.Num(); ++q2)
            {
                gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q2]);
            }
        }
    }
    break;
    case ELinkStyle::Linear:
    {
        for (INT q1 = 0; q1 < linkQubits.Num() - 1; ++q1)
        {
            gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q1 + 1]);
        }
    }
    break;
    case ELinkStyle::Circular:
    {
        for (INT q1 = 0; q1 < linkQubits.Num(); ++q1)
        {
            if (q1 == linkQubits.Num() - 1)
            {
                gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[0]);
            }
            else
            {
                gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q1 + 1]);
            }
        }
    }
    break;
    case ELinkStyle::PairWise:
    {
        if (1 & level)
        {
            for (INT q1 = 1; q1 < linkQubits.Num() - 1; q1 += 2)
            {
                gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q1 + 1]);
            }
        }
        else
        {
            //0-1, 2-3
            for (INT q1 = 0; q1 < linkQubits.Num() - 1; q1 += 2)
            {
                gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q1 + 1]);
            }
        }
    }
    break;
    case ELinkStyle::DoublePairWise:
    {
        for (INT q1 = 1; q1 < linkQubits.Num() - 1; q1 += 2)
        {
            gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q1 + 1]);
        }

        for (INT q1 = 0; q1 < linkQubits.Num() - 1; q1 += 2)
        {
            gate.AppendGate(bCx ? cx : cz, linkQubits[q1], linkQubits[q1 + 1]);
        }
    }
    break;
    case ELinkStyle::SCA:
    {
        for (INT q1 = 0; q1 < linkQubits.Num(); ++q1)
        {
            //level 0:
            //0 - end
            //1 - 0
            //2 - 1
            //level 1:
            //2 - 1
            //0 - end
            //1 - 0
            INT iq1 = -static_cast<INT>(level) + static_cast<INT>(q1);
            INT iq2 = iq1 - 1;
            while (iq1 < 0)
            {
                iq1 += linkQubits.Num();
            }
            while (iq2 < 0)
            {
                iq2 += linkQubits.Num();
            }
            while (iq1 >= linkQubits.Num())
            {
                iq1 -= linkQubits.Num();
            }
            while (iq2 >= linkQubits.Num())
            {
                iq2 -= linkQubits.Num();
            }

            gate.AppendGate(bCx ? cx : cz, linkQubits[iq1], linkQubits[iq2]);
        }
    }
    break;
    default:
        appCrucial(_T("Link style %d not implemented!\n"), eStype);
        break;
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================