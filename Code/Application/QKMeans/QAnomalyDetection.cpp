//=============================================================================
// FILENAME : QAnomalyDetection.cpp
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#include "QKMeans.h"

void TestFitPointSet(CParameters& params)
{
    CCString sValues;
    INT iValues;
    Real fValues;

    CCString sTrainingPoints;
    __FetchStringWithDefault(_T("PointFile"), _T(""));
    sTrainingPoints = sValues;

    CCString sAnsatzFile;
    __FetchStringWithDefault(_T("AnsatzFile"), _T(""));
    sAnsatzFile = sValues;

    CCString sHistory;
    __FetchStringWithDefault(_T("HistoryFile"), _T(""));
    sHistory = sValues;

    UINT uiLevel = 1;
    __FetchIntWithDefault(_T("AnsatzLevel"), 2);
    uiLevel = static_cast<UINT>(iValues);

    UINT uiMaxStep = 10000;
    __FetchIntWithDefault(_T("MaxStep"), 10000);
    uiMaxStep = static_cast<UINT>(iValues);

    UBOOL bOnlyReal = FALSE;
    __FetchIntWithDefault(_T("OnlyReal"), 0);
    bOnlyReal = (0 != iValues);

    Real fLearnRate = F(0.001);
    __FetchRealWithDefault(_T("LearnRate"), F(0.001));
    fLearnRate = fValues;

    Real fGoal = F(0.1);
    __FetchRealWithDefault(_T("Goal"), F(0.1));
    fGoal = fValues;

    FitAE(sTrainingPoints, sAnsatzFile, sHistory, uiLevel, fLearnRate, fGoal, uiMaxStep, bOnlyReal);
}

void TestFitPointSetAdap(CParameters& params)
{
    CCString sValues;
    INT iValues;
    Real fValues;

    CCString sTrainingPoints;
    __FetchStringWithDefault(_T("PointFile"), _T(""));
    sTrainingPoints = sValues;

    CCString sAnsatzFile;
    __FetchStringWithDefault(_T("AnsatzFile"), _T(""));
    sAnsatzFile = sValues;

    CCString sHistory;
    __FetchStringWithDefault(_T("HistoryFile"), _T(""));
    sHistory = sValues;

    UINT uiAdap = 200;
    __FetchIntWithDefault(_T("AdaptiveWait"), 200);
    uiAdap = static_cast<UINT>(iValues);

    UINT uiMaxStep = 10000;
    __FetchIntWithDefault(_T("MaxStep"), 10000);
    uiMaxStep = static_cast<UINT>(iValues);

    UINT uiMaxLayer = 100;
    __FetchIntWithDefault(_T("MaxLayer"), 100);
    uiMaxLayer = static_cast<UINT>(iValues);

    UBOOL bOnlyReal = FALSE;
    __FetchIntWithDefault(_T("OnlyReal"), 0);
    bOnlyReal = (0 != iValues);

    Real fLearnRate = F(0.001);
    __FetchRealWithDefault(_T("LearnRate"), F(0.001));
    fLearnRate = fValues;

    Real fAdapEps = F(0.001);
    __FetchRealWithDefault(_T("AdaptiveEps"), F(0.001));
    fAdapEps = fValues;

    Real fGoal = F(0.1);
    __FetchRealWithDefault(_T("Goal"), F(0.1));
    fGoal = fValues;

    FitAE(sTrainingPoints, sAnsatzFile, sHistory, fLearnRate, fGoal, uiMaxStep, uiMaxLayer, uiAdap, fAdapEps, bOnlyReal);
}

//=============================================================================
// END OF FILE
//=============================================================================