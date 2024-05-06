//=============================================================================
// FILENAME : QuantumFit.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/10/2022 nbale]
//=============================================================================

#include "SwapTest.h"

TArray<QLComplex> ExchangeComplexBuffer(const TArray<QLComplex>& row)
{
    TArray<QLComplex> ret;
    for (INT i = 0; i < (row.Num() + 1) / 2; ++i)
    {
        if (row.Num() >= 2 * (i + 1))
        {
            ret.AddItem(_make_cuComplex(row[2 * i].x, row[2 * i + 1].x));
        }
        else
        {
            ret.AddItem(_make_cuComplex(row[2 * i].x, F(0.0)));
        }
    }
    return ret;
}

TArray<Real> QuantumOverlap(const TArray<QLComplex>& row1, const TArray<QLComplex>& row2, UBOOL bComplex, const QLSimulatorParametersDensityMatrix& param, INT iRealMeasure, UBOOL bEnableNoise)
{
    QLSimulatorParametersDensityMatrix param2 = param;
    if (bComplex)
    {
        QLGate zerotest = ZeroTest(row1, row2);
        param2.m_byQubitCount = static_cast<BYTE>(zerotest.m_lstQubits.Num());
        param2.m_MasterGate = zerotest;
        param2.m_lstMeasureQubits.Append(zerotest.m_lstQubits);
    }
    else
    {
        QLGate zerotest = ZeroTestReal(row1, row2);
        param2.m_byQubitCount = static_cast<BYTE>(zerotest.m_lstQubits.Num());
        param2.m_MasterGate = zerotest;
        param2.m_lstMeasureQubits.Append(zerotest.m_lstQubits);
    }

    TArray<UINT> gates = param2.m_MasterGate.SummarizeGateCounts();
    appGeneral(_T("Single: %d, TwoQubit: %d, Noise: %d, Measure: %d, Other: %d\n"), gates[0], gates[1], gates[2], gates[3], gates[4]);

    Real fAveragePurity = F(1.0);
    Real fAverageFidelity = F(1.0);
    Real fRes = F(0.0);

    if (bEnableNoise)
    {
        param2.m_iMeasureTimes = iRealMeasure;
        QLSimulatorOutputDensityMatrix out;
        QLSimulatorDensityMatrix sim;
        sim.Simulate(&param2, &out);

        fAveragePurity = out.AveragePurity();
        fAverageFidelity = out.AverageFidelity();
        fRes = out.m_lstMeasureOutcomes[0];
    }
    else
    {
        QLSimulatorMeasure sim;
        QLSimulatorOutputMeasure out;
        QLSimulatorParametersMeasure param3;
        param3.m_bPrint = FALSE;
        param3.m_iRepeat = -1;
        param3.m_iMeasureUntil = -1;
        param3.m_byQubitCount = param2.m_byQubitCount;
        param3.m_MasterGate = param2.m_MasterGate;
        param3.m_lstMeasureBits = param2.m_lstMeasureQubits;

        sim.Simulate(&param3, &out);

        fRes = out.m_lstMeasureOutcomes[0];

        if (iRealMeasure > 0)
        {
            UINT iMeasure0 = 0;
            for (INT i = 0; i < iRealMeasure; ++i)
            {
                if (RandomF() < out.m_lstMeasureOutcomes[0])
                {
                    ++iMeasure0;
                }
            }
            fRes = iMeasure0 / static_cast<Real>(iRealMeasure);
        }
    }

    if (fRes < _QL_FLT_MIN)
    {
        fRes = F(0.0);
    }
    else
    {
        fRes = _sqrt(fRes);
    }

    appGeneral(_T("average purity: %f, average fidelity: %f, res: %f \n"), fAveragePurity, fAverageFidelity, fRes);

    TArray<Real> ret;
    ret.AddItem(fAveragePurity);
    ret.AddItem(fAverageFidelity);
    ret.AddItem(fRes);

    return ret;
}

Real ClassicalDot(const TArray<QLComplex>& row1, const TArray<QLComplex>& row2)
{
    QLMatrix m1 = QLMatrix::CopyCreate(1, row1.Num(), row1.GetData());
    QLMatrix m2 = QLMatrix::CopyCreate(1, row2.Num(), row2.GetData());

    m1 = m1 / _sqrt(_cuCabsf(m1.VectorDot(m1)));
    m2 = m2 / _sqrt(_cuCabsf(m2.VectorDot(m2)));

    return _cuCabsf(m1.VectorDot(m2));
}

int main()
{
    QLRandomInitializer random(ERandom::ER_XORWOW, appGetTimeStamp());
    CParameters params;
    CYAMLParser::ParseFile(_T("../SwapTest.yaml"), params);
    params.Dump();

    CCString sValues;
    __FetchStringWithDefault(_T("FileName1"), _T(""));

    //Load Data File
    CCString sFile1 = sValues;
    QLMatrix m1 = ReadCSVR(sValues);
    m1.Print("v1");

    __FetchStringWithDefault(_T("FileName2"), _T(""));
    QLMatrix m2 = ReadCSVR(sValues);
    m2.Print("v2");

    UBOOL bSameFile = sFile1 == sValues;

    QLSimulatorParametersDensityMatrix simulateParam;
    simulateParam.m_bMeasureFidelity = TRUE;
    simulateParam.m_bMeasurePurity = TRUE;
    simulateParam.m_bPrint = FALSE;
    
    Real fValues = F(0.0);
    INT iValues = 0;

    __FetchRealWithDefault(_T("DampingAfterGate"), F(0.0));
    simulateParam.m_fDampingAfterGate = fValues;
    __FetchRealWithDefault(_T("DepolarisingAfterGate"), F(0.0));
    simulateParam.m_fDepolarisingAfterGate = fValues;
    __FetchRealWithDefault(_T("DephasingAfterGate"), F(0.0));
    simulateParam.m_fDephaseAfterGate = fValues;
    __FetchRealWithDefault(_T("MixPauliAfterGate"), F(0.0));
    simulateParam.m_fMixPauliAfterGate = fValues;
    __FetchRealWithDefault(_T("TwoQubitDephasing"), F(0.0));
    simulateParam.m_fTwoDephaseAfterGate = fValues;
    __FetchRealWithDefault(_T("TwoQubitDepolarising"), F(0.0));
    simulateParam.m_fTwoDepolarisingAfterGate = fValues;

    __FetchRealWithDefault(_T("DampingBeforeMeasure"), F(0.0));
    simulateParam.m_fDampingBeforeMeasure = fValues;
    __FetchRealWithDefault(_T("DepolarisingBeforeMeasure"), F(0.0));
    simulateParam.m_fDepolarisingBeforeMeasure = fValues;
    __FetchRealWithDefault(_T("DephasingBeforeMeasure"), F(0.0));
    simulateParam.m_fDephaseBeforeMeasure = fValues;
    __FetchRealWithDefault(_T("MixPauliBeforeMeasure"), F(0.0));
    simulateParam.m_fMixPauliBeforeMeasure = fValues;

    __FetchIntWithDefault(_T("MeasureRepeat"), 100);
    simulateParam.m_iMeasureTimes = -1;
    INT iRealMeasure = iValues;

    __FetchIntWithDefault(_T("Complex"), 1);
    UBOOL bComplex = (0 != iValues);

    __FetchIntWithDefault(_T("EnableClassical"), 1);
    UBOOL bHasClassical = (0 != iValues);

    __FetchIntWithDefault(_T("EnableNoise"), 1);
    UBOOL bEnableNoise = (0 != iValues);

    QLComplex* expectedoutput = NULL;
    if (bHasClassical)
    {
        expectedoutput = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m1.Y() * m2.Y()));
    }
    QLComplex* output = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m1.Y() * m2.Y()));
    QLComplex* averagePurity = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m1.Y() * m2.Y()));
    QLComplex* averageFidelity = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m1.Y() * m2.Y()));
    UINT iNow = 0;
    for (INT i = 0; i < static_cast<INT>(m1.Y()); ++i)
    {
        for (INT j = 0; j < static_cast<INT>(m2.Y()); ++j)
        {
            UINT idx1 = i * m2.Y() + j;
            UINT idx2 = j * m2.Y() + i;
            if (bSameFile)
            {
                if (i == j)
                {
                    if (bHasClassical)
                    {
                        expectedoutput[idx1] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    output[idx1] = _make_cuComplex(F(1.0), F(0.0));
                    averagePurity[idx1] = _make_cuComplex(F(1.0), F(0.0));
                    averageFidelity[idx1] = _make_cuComplex(F(1.0), F(0.0));
                    continue;
                }
                else if (j < i)
                {
                    if (bHasClassical)
                    {
                        expectedoutput[idx1] = expectedoutput[idx2];
                    }
                    output[idx1] = output[idx2];
                    averagePurity[idx1] = averagePurity[idx2];
                    averageFidelity[idx1] = averageFidelity[idx2];
                    continue;
                }
            }
            TArray<QLComplex> v1 = m1.GetLine(i);
            TArray<QLComplex> v2 = m2.GetLine(j);
            if (bComplex)
            {
                v1 = ExchangeComplexBuffer(v1);
                v2 = ExchangeComplexBuffer(v2);
            }

            if (bHasClassical)
            {
                Real fClassical = ClassicalDot(v1, v2);
                ++iNow;
                appGeneral(_T("=============\nCalculating v1[%d] * v2[%d] (progress %d/%d), expected result: %f\n"),
                    i,
                    j,
                    iNow,
                    bSameFile ? ((m1.Y() * m1.Y() - m1.Y()) / 2) : m1.Y() * m2.Y(),
                    fClassical);
                expectedoutput[idx1] = _make_cuComplex(fClassical, F(0.0));
            }
            else
            {
                ++iNow;
                appGeneral(_T("=============\nCalculating v1[%d] * v2[%d] (progress %d/%d)\n"),
                    i,
                    j,
                    iNow,
                    bSameFile ? ((m1.Y() * m1.Y() - m1.Y()) / 2) : m1.Y() * m2.Y());
            }

            TArray<Real> quantumRes = QuantumOverlap(v1, v2, bComplex, simulateParam, iRealMeasure, bEnableNoise);
            
            output[idx1] = _make_cuComplex(quantumRes[2], F(0.0));
            averagePurity[idx1] = _make_cuComplex(quantumRes[0], F(0.0));
            averageFidelity[idx1] = _make_cuComplex(quantumRes[1], F(0.0));
        }
    }

    if (bHasClassical)
    {
        QLMatrix expm(m1.Y(), m2.Y(), expectedoutput);
        if (m1.Y() * m2.Y() < 10000)
        {
            expm.Print("ClassicalResult");
        }
        __FetchStringWithDefault(_T("ClassicalResSaveFileName"), _T(""));
        if (!sValues.IsEmpty())
        {
            SaveCSVR(expm, sValues);
            appGeneral(_T("%f file saved...\n"), sValues.c_str());
        }
    }
    
    QLMatrix resm(m1.Y(), m2.Y(), output);
    QLMatrix purity(m1.Y(), m2.Y(), averagePurity);
    QLMatrix fidelity(m1.Y(), m2.Y(), averageFidelity);

    if (m1.Y() * m2.Y() < 10000)
    {
        resm.Print("QuantumResult");
        purity.Print("AveragePurity");
        fidelity.Print("AverageFidelity");
    }

    __FetchStringWithDefault(_T("QuantumResSaveFileName"), _T(""));
    if (!sValues.IsEmpty())
    {
        SaveCSVR(resm, sValues);
        appGeneral(_T("%f file saved...\n"), sValues.c_str());
    }
    __FetchStringWithDefault(_T("AveragePuritySaveFileName"), _T(""));
    if (!sValues.IsEmpty())
    {
        SaveCSVR(purity, sValues);
        appGeneral(_T("%f file saved...\n"), sValues.c_str());
    }
    __FetchStringWithDefault(_T("AverageFidelitySaveFileName"), _T(""));
    if (!sValues.IsEmpty())
    {
        SaveCSVR(fidelity, sValues);
        appGeneral(_T("%f file saved...\n"), sValues.c_str());
    }

    return 0;
}


//=============================================================================
// END OF FILE
//=============================================================================