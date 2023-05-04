//=============================================================================
// FILENAME : CircuitTest.cpp
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [16/03/2023 nbale]
//=============================================================================

#include "QKMeans.h"

/**
*
*/


void TestDistance()
{
    TestDistance(_T("p510d150.csv"), 510, F(15.0), 100000);
}

void TestProbalityToBuild()
{
    Real fBuildProb, fMeasureProb;
    TestProbalityToBuildAmplitude(6, 8, F(1.5), fBuildProb, fMeasureProb, TRUE);
}

void TestProbalityMeanP6()
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../QKMeans.yaml"), params);
    TArray<UINT> testCount;
    params.FetchValueArrayUINT(_T("P6TestList"), testCount);

    UINT vectorLen = 6;
    Real fFactor = F(1.5);
    appPushLogDate(FALSE);
    for (UINT uiK = 4; uiK < 17; ++uiK)
    {
        TArray<Real> probs;
        appGeneral(_T("uiK = %d\n"), uiK);
        for (UINT i = 0; i < testCount[uiK - 4]; ++i)
        {
            Real fBuildProb, fMeasureProb;
            TestProbalityToBuildAmplitude(vectorLen, uiK, fFactor, fBuildProb, fMeasureProb, FALSE);
            probs.AddItem(fBuildProb);
            probs.AddItem(fMeasureProb);
            if (0 == i % 50)
            {
                appGeneral(_T("\n"));
            }
            appGeneral(_T("="));
        }
        appGeneral(_T("\n"));
        CCString sFileName;
        sFileName.Format(_T("p%dk%d.csv"), vectorLen, uiK);
        SaveCSVAR(probs.GetData(), 2, testCount[uiK - 4], sFileName);
    }

    appPopLogDate();
}

/*
void TestProbalityMeanK6()
{
    //CParameters params;
    //CYAMLParser::ParseFile(_T("../QKMeans.yaml"), params);
    //TArray<UINT> testCount;
    //params.FetchValueArrayUINT(_T("P6TestList"), testCount);
    UINT uiK = 6;
    UINT vectorLen = 6;
    Real fFactor = F(1.5);
    appPushLogDate(FALSE);
    for (UINT uiPPower = 2; uiPPower < 10; ++uiPPower)
    {
        TArray<Real> probs;
        appGeneral(_T("uiK = %d\n"), uiK);
        for (UINT i = 0; i < testCount[uiK - 4]; ++i)
        {
            Real fBuildProb, fMeasureProb;
            TestProbalityToBuildAmplitude(vectorLen, uiK, fFactor, fBuildProb, fMeasureProb, FALSE);
            probs.AddItem(fBuildProb);
            probs.AddItem(fMeasureProb);
            if (0 == i % 50)
            {
                appGeneral(_T("\n"));
            }
            appGeneral(_T("="));
        }
        appGeneral(_T("\n"));
        CCString sFileName;
        sFileName.Format(_T("p%dk%d.csv"), vectorLen, uiK);
        SaveCSVAR(probs.GetData(), 2, testCount[uiK - 4], sFileName);
    }

    appPopLogDate();
}
*/

//=============================================================================
// END OF FILE
//=============================================================================