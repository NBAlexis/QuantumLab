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
#include <random>

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

TArray<Real> QuantumOverlap(const TArray<QLComplex>& row1, const TArray<QLComplex>& row2, UBOOL bComplex, const QLSimulatorParametersDensityMatrix& param, INT iRealMeasure)
{
    QLSimulatorParametersDensityMatrix param2 = param;
    if (bComplex)
    {
        QLGate cswap = CreateSwapTest(row1, row2);
        param2.m_byQubitCount = static_cast<BYTE>(cswap.m_lstQubits.Num());
        param2.m_MasterGate = cswap;
    }
    else
    {
        QLGate cswap = CreateSwapTestReal(row1, row2);
        param2.m_byQubitCount = static_cast<BYTE>(cswap.m_lstQubits.Num());
        param2.m_MasterGate = cswap;
    }

    param2.m_lstMeasureQubits.AddItem(0);

    QLSimulatorOutputDensityMatrix out;
    QLSimulatorDensityMatrix sim;
    sim.Simulate(&param2, &out);

    Real fAveragePurity = out.AveragePurity();
    Real fAverageFidelity = out.AverageFidelity();
    Real fRes = 2 * out.m_lstMeasureOutcomes[0] - 1;
    if (iRealMeasure > 0)
    {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0, 1.0);

        UINT iMeasure0 = 0;
        for (INT i = 0; i < iRealMeasure; ++i)
        {
            //appGeneral(_T("%f %f\n"), dis(gen), out.m_lstMeasureOutcomes[0]);
            if (dis(gen) < out.m_lstMeasureOutcomes[0])
            {
                ++iMeasure0;
            }
        }
        fRes = iMeasure0 / static_cast<Real>(iRealMeasure);
        fRes = 2 * fRes - 1;
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

    m1 = m1 / _sqrt(cuCabsf(m1.VectorDot(m1)));
    m2 = m2 / _sqrt(cuCabsf(m2.VectorDot(m2)));

    return cuCabsf(m1.VectorDot(m2));
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

    TArray<QLComplex> expectedoutput;
    TArray<QLComplex> output;
    TArray<QLComplex> averagePurity;
    TArray<QLComplex> averageFidelity;
    UINT iNow = 0;
    for (INT i = 0; i < static_cast<INT>(m1.Y()); ++i)
    {
        for (INT j = 0; j < static_cast<INT>(m2.Y()); ++j)
        {
            if (bSameFile)
            {
                if (i == j)
                {
                    expectedoutput.AddItem(_make_cuComplex(F(1.0), F(0.0)));
                    output.AddItem(_make_cuComplex(F(1.0), F(0.0)));
                    averagePurity.AddItem(_make_cuComplex(F(1.0), F(0.0)));
                    averageFidelity.AddItem(_make_cuComplex(F(1.0), F(0.0)));
                    continue;
                }
                else if (j < i)
                {
                    UINT idx = j * m2.Y() + i;
                    expectedoutput.AddItem(expectedoutput[idx]);
                    output.AddItem(output[idx]);
                    averagePurity.AddItem(averagePurity[idx]);
                    averageFidelity.AddItem(averageFidelity[idx]);
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
            Real fClassical = ClassicalDot(v1, v2);
            ++iNow;
            appGeneral(_T("=============\nCalculating v1[%d] * v2[%d] (progress %d/%d), expected result: %f\n"),
                i, 
                j, 
                iNow,
                bSameFile ? ((m1.Y() * m1.Y() - m1.Y()) / 2) : m1.Y() * m2.Y(),
                fClassical);

            TArray<Real> quantumRes = QuantumOverlap(v1, v2, bComplex, simulateParam, iRealMeasure);

            expectedoutput.AddItem(_make_cuComplex(fClassical, F(0.0)));
            output.AddItem(_make_cuComplex(quantumRes[2], F(0.0)));
            averagePurity.AddItem(_make_cuComplex(quantumRes[0], F(0.0)));
            averageFidelity.AddItem(_make_cuComplex(quantumRes[1], F(0.0)));
        }
    }

    QLMatrix expm = QLMatrix::CopyCreate(m1.Y(), m2.Y(), expectedoutput.GetData());
    QLMatrix resm = QLMatrix::CopyCreate(m1.Y(), m2.Y(), output.GetData());
    QLMatrix purity = QLMatrix::CopyCreate(m1.Y(), m2.Y(), averagePurity.GetData());
    QLMatrix fidelity = QLMatrix::CopyCreate(m1.Y(), m2.Y(), averageFidelity.GetData());

    expm.Print("ClassicalResult");
    resm.Print("QuantumResult");
    purity.Print("AveragePurity");
    fidelity.Print("AverageFidelity");

    __FetchStringWithDefault(_T("ClassicalResSaveFileName"), _T(""));
    if (!sValues.IsEmpty())
    {
        SaveCSVR(expm, sValues);
        appGeneral(_T("%f file saved...\n"), sValues.c_str());
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