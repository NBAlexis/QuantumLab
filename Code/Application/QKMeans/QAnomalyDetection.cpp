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

    UBOOL bAbsolute = TRUE;
    __FetchIntWithDefault(_T("UseAbsorlute"), 1);
    bAbsolute = (0 != iValues);

    Real fLearnRate = F(0.001);
    __FetchRealWithDefault(_T("LearnRate"), F(0.001));
    fLearnRate = fValues;

    Real fGoal = F(0.1);
    __FetchRealWithDefault(_T("Goal"), F(0.1));
    fGoal = fValues;

    __FetchStringWithDefault(_T("AnsatzStyle"), _T("PairWise"));
    ELinkStyle ansatzStyle = __STRING_TO_ENUM(ELinkStyle, sValues);

    __FetchStringWithDefault(_T("AnsatzLayer"), _T("CZ"));
    ELinkLayer ansatzLayer = __STRING_TO_ENUM(ELinkLayer, sValues);

    __FetchStringWithDefault(_T("AnsatzSingleLayer"), _T("RYRZ"));
    ESingleLayer ansatzSingleLayer = __STRING_TO_ENUM(ESingleLayer, sValues);

    __FetchStringWithDefault(_T("AnsatzInitial"), _T("MBL"));
    EAnsatzInitial ansatzInital = __STRING_TO_ENUM(EAnsatzInitial, sValues);

    CSpliteAngleBufferHelper helper;

    FitAE(sTrainingPoints, sAnsatzFile, sHistory, uiLevel, ansatzStyle, ansatzSingleLayer, ansatzLayer, ansatzInital, fLearnRate, fGoal, uiMaxStep, bAbsolute);
}

void TestFitPointSetSE(CParameters& params)
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

    BYTE uiEncodeBits = 4;
    __FetchIntWithDefault(_T("NumberOfEncode"), 4);
    uiEncodeBits = static_cast<BYTE>(iValues);

    UINT uiLevel = 1;
    __FetchIntWithDefault(_T("AnsatzLevel"), 2);
    uiLevel = static_cast<UINT>(iValues);

    UINT uiMaxStep = 10000;
    __FetchIntWithDefault(_T("MaxStep"), 10000);
    uiMaxStep = static_cast<UINT>(iValues);

    UBOOL bAbsolute = TRUE;
    __FetchIntWithDefault(_T("UseAbsorlute"), 1);
    bAbsolute = (0 != iValues);

    Real fLearnRate = F(0.001);
    __FetchRealWithDefault(_T("LearnRate"), F(0.001));
    fLearnRate = fValues;

    Real fGoal = F(0.1);
    __FetchRealWithDefault(_T("Goal"), F(0.1));
    fGoal = fValues;

    __FetchStringWithDefault(_T("AnsatzStyle"), _T("PairWise"));
    ELinkStyle ansatzStyle = __STRING_TO_ENUM(ELinkStyle, sValues);

    __FetchStringWithDefault(_T("SEStyle"), _T("PairWise"));
    ELinkStyle simpleencodeStyle = __STRING_TO_ENUM(ELinkStyle, sValues);

    __FetchStringWithDefault(_T("AnsatzLayer"), _T("CZ"));
    ELinkLayer ansatzLayer = __STRING_TO_ENUM(ELinkLayer, sValues);

    __FetchStringWithDefault(_T("SELayer"), _T("CZ"));
    ELinkLayer simpleencodeLayer = __STRING_TO_ENUM(ELinkLayer, sValues);

    __FetchStringWithDefault(_T("AnsatzSingleLayer"), _T("RYRZ"));
    ESingleLayer ansatzSingleLayer = __STRING_TO_ENUM(ESingleLayer, sValues);

    __FetchStringWithDefault(_T("AnsatzInitial"), _T("MBL"));
    EAnsatzInitial ansatzInital = __STRING_TO_ENUM(EAnsatzInitial, sValues);

    CSpliteAngleBufferHelper helper;

    FitSE(sTrainingPoints, sAnsatzFile, sHistory, uiEncodeBits, uiLevel, ansatzStyle, ansatzSingleLayer, ansatzLayer, ansatzInital, simpleencodeStyle, simpleencodeLayer, fLearnRate, fGoal, uiMaxStep, bAbsolute);
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

    UBOOL bAbsolute = TRUE;
    __FetchIntWithDefault(_T("UseAbsorlute"), 1);
    bAbsolute = (0 != iValues);

    Real fLearnRate = F(0.001);
    __FetchRealWithDefault(_T("LearnRate"), F(0.001));
    fLearnRate = fValues;

    Real fAdapEps = F(0.001);
    __FetchRealWithDefault(_T("AdaptiveEps"), F(0.001));
    fAdapEps = fValues;

    Real fGoal = F(0.1);
    __FetchRealWithDefault(_T("Goal"), F(0.1));
    fGoal = fValues;

    __FetchStringWithDefault(_T("AnsatzStyle"), _T("PairWise"));
    ELinkStyle ansatzStyle = __STRING_TO_ENUM(ELinkStyle, sValues);

    __FetchStringWithDefault(_T("AnsatzLayer"), _T("CZ"));
    ELinkLayer ansatzLayer = __STRING_TO_ENUM(ELinkLayer, sValues);

    __FetchStringWithDefault(_T("AnsatzSingleLayer"), _T("RYRZ"));
    ESingleLayer ansatzSingleLayer = __STRING_TO_ENUM(ESingleLayer, sValues);

    __FetchStringWithDefault(_T("AnsatzInitial"), _T("MBL"));
    EAnsatzInitial ansatzInital = __STRING_TO_ENUM(EAnsatzInitial, sValues);

    CSpliteAngleBufferHelper helper;

    FitAE(sTrainingPoints, sAnsatzFile, sHistory, ansatzStyle, ansatzSingleLayer, ansatzLayer, ansatzInital, fLearnRate, fGoal, uiMaxStep, uiMaxLayer, uiAdap, fAdapEps, bAbsolute);
}

//=============================================================================
// END OF FILE
//=============================================================================