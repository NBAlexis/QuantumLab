//=============================================================================
// FILENAME : AmplitudeEncode.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

TArray<QLComplex> QLAPI NormalizeV(const TArray<QLComplex>& v, UINT& lenPower)
{
    TArray<QLComplex> ret = v;
    UINT vLength = static_cast<UINT>(ret.Num());
    UINT vLengthPower = MostSignificantPowerTwo(vLength);
    if (vLengthPower < 1)
    {
        vLengthPower = 1;
    }
    UINT vLengthWanted = (1 << vLengthPower);
    for (UINT i = vLength; i < vLengthWanted; ++i)
    {
        ret.AddItem(_make_cuComplex(F(0.0), F(0.0)));
    }

    Real sum = F(0.0);
    for (INT i = 0; i < ret.Num(); ++i)
    {
        sum += __cuCabsSqf(ret[i]);
    }

    if (sum < _QL_FLT_MIN_)
    {
        appCrucial(_T("vector length too small!\n"));
        sum = _QL_FLT_MIN_;
    }
    sum = _sqrt(sum);
    for (INT i = 0; i < ret.Num(); ++i)
    {
        ret[i] = cuCdivf_cr_host(ret[i], sum);
    }

    lenPower = vLengthPower;
    return ret;
}

TArray<Real> QLAPI NormalizeVReal(const TArray<Real>& v, UINT& lenPower)
{
    TArray<Real> ret;
    UINT vLength = static_cast<UINT>(v.Num());
    UINT vLengthPower = MostSignificantPowerTwo(vLength);
    if (vLengthPower < 1)
    {
        vLengthPower = 1;
    }
    //UINT vLengthWanted = (1 << vLengthPower);

    Real sum = F(0.0);
    for (INT i = 0; i < v.Num(); ++i)
    {
        sum += v[i] * v[i];
    }

    if (sum < _QL_FLT_MIN_)
    {
        appCrucial(_T("vector length too small!\n"));
        sum = _QL_FLT_MIN_;
    }
    sum = _sqrt(sum);
    for (INT i = 0; i < v.Num(); ++i)
    {
        ret.AddItem(v[i] / sum);
    }

    lenPower = vLengthPower;
    return ret;
}

TArray<Real> CalculateVectorAbsoluteRotations(TArray<QLComplex> v, UINT lengthPower)
{
    TArray<Real> degreeList;
    for (UINT i = 0; i < lengthPower; ++i)
    {
        UINT degreeCount = 1U << (lengthPower - i - 1);
        UINT skip = 1U << i;
        UINT stride = 1U << (i + 1);
        for (UINT j = 0; j < degreeCount; ++j)
        {
            Real cs = _cuCabsf(v[j * stride]);
            Real sn = _cuCabsf(v[j * stride + skip]);
            Real degree = _atan2(sn, cs);
            degreeList.AddItem(degree);
            Real csd = _cos(degree);
            Real snd = _sin(degree);
            if (csd > _QL_FLT_MIN_)
            {
                for (UINT k = 0; k < skip; ++k)
                {
                    v[j * stride + k] = cuCdivf_cr_host(v[j * stride + k], csd);
                }
            }

            if (snd > _QL_FLT_MIN_)
            {
                for (UINT k = 0; k < skip; ++k)
                {
                    v[j * stride + skip + k] = cuCdivf_cr_host(v[j * stride + skip + k], snd);
                }
            }
        }
    }
    return degreeList;
}

TArray<Real> CalculateVectorAbsoluteRotationsReal(TArray<Real> v, UINT lengthPower)
{
    TArray<Real> degreeList;
    for (UINT i = 0; i < lengthPower; ++i)
    {
        UINT degreeCount = 1U << (lengthPower - i - 1);
        UINT skip = 1U << i;
        UINT stride = 1U << (i + 1);
        for (UINT j = 0; j < degreeCount; ++j)
        {
            Real cs = v[j * stride];
            Real sn = v[j * stride + skip];
            Real degree = _atan2(sn, cs);
            degreeList.AddItem(degree);
            Real csd = _cos(degree);
            Real snd = _sin(degree);
            if (csd > _QL_FLT_MIN_)
            {
                for (UINT k = 0; k < skip; ++k)
                {
                    v[j * stride + k] = v[j * stride + k] / csd;
                }
            }

            if (snd > _QL_FLT_MIN_)
            {
                for (UINT k = 0; k < skip; ++k)
                {
                    v[j * stride + skip + k] = v[j * stride + skip + k] / snd;
                }
            }
        }
    }
    return degreeList;
}

void MakeCircuitWithRotations(QLGate& gate, const TArray<Real>& degreeArray, UINT lenPower)
{
    UINT len = static_cast<UINT>(degreeArray.Num());
    UINT vLengthPower = MostSignificantPowerTwo(len);
    if (len != ((1U << vLengthPower) - 1))
    {
        appCrucial(_T("Some thing bad happens! MakeCircuitWithRotations"));
        return;
    }
    
    UINT degreeIndex = 0;
    for (UINT i = 0; i < vLengthPower; ++i)
    {
        if (0 == i)
        {
            Real degree = degreeArray[len - 1 - degreeIndex];
            degreeIndex = degreeIndex + 1;

            QLGate ry(EBasicOperation::EBO_RY, degree * F(2.0));
            TArray<BYTE> bit;
            bit.AddItem(static_cast<BYTE>(vLengthPower - 1));
            gate.AppendGate(ry, bit);
        }
        else
        {
            UINT degreeCount = 1U << i;
            TArray<BYTE> controller;
            TArray<Real> degrees;
            for (UINT j = 0; j < degreeCount; ++j)
            {
                degrees.InsertAt(0, degreeArray[len - 1 - degreeIndex] * F(2.0));
                degreeIndex = degreeIndex + 1;
            }

            for (UINT k = 0; k < i; ++k)
            {
                controller.AddItem(static_cast<BYTE>(vLengthPower - 1 - k));
            }
            controller.AddItem(static_cast<BYTE>(vLengthPower - 1 - i));
            QLGate fry = FRy(degrees, static_cast<UINT>(controller.Num()));
            gate.AppendGate(fry, controller);
        }
    }
}

void ApplyPhase(QLGate& gate, const TArray<QLComplex>& v, UINT lenPower)
{
    TArray<Real> phase;
    for (INT i = 0; i < v.Num(); ++i)
    {
        phase.AddItem(__cuCargf(v[i]));
    }
    CSDMatrix csd;
    csd.SetAsDiagonalUMatrix(lenPower, lenPower, phase);
    TArray<CSDMatrix> csdRotations = csd.DecomposeUIntoRZ();
    for (INT i = 0; i < csdRotations.Num(); ++i)
    {
        csdRotations[csdRotations.Num() - 1 - i].Gate(gate);
    }
}

//QLGate QLAPI AmplitudeEncode(const TArray<QLComplex>& v)
//{
//    UINT vLengthPower = 0;
//    TArray<QLComplex> varray = NormalizeV(v, vLengthPower);
//
//    QLGate ret;
//    ret.AddQubits(static_cast<BYTE>(vLengthPower));
//    ret.m_sName = _T("AmpEnc");
//    MakeCircuitWithRotations(ret, CalculateVectorAbsoluteRotations(varray, vLengthPower), vLengthPower);
//    ApplyPhase(ret, varray, vLengthPower);
//    return ret;
//}
//
//QLGate QLAPI AmplitudeEncodeReal(const TArray<Real>& v)
//{
//    UINT vLengthPower = 0;
//    TArray<Real> varray = NormalizeVReal(v, vLengthPower);
//
//    QLGate ret;
//    ret.AddQubits(static_cast<BYTE>(vLengthPower));
//    ret.m_sName = _T("AmpEnc");
//    MakeCircuitWithRotations(ret, CalculateVectorAbsoluteRotationsReal(varray, vLengthPower), vLengthPower);
//    return ret;
//}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================