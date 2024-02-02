//=============================================================================
// FILENAME : QuantumKmeans.cpp
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [21/12/2022 nbale]
//=============================================================================

#include "QKMeans.h"

void TestKMeans2D(CParameters& params)
{
    CCString sValues;
    INT iValues;

    CCString sLoadFile;
    __FetchStringWithDefault(_T("LoadFile"), _T(""));
    sLoadFile = sValues;

    CCString sKFileHead;
    __FetchStringWithDefault(_T("KFileHead"), _T(""));
    sKFileHead = sValues;

    CCString sCenterFileHead;
    __FetchStringWithDefault(_T("CenterFileHead"), _T(""));
    sCenterFileHead = sValues;

    CCString sReepatFileHead;
    __FetchStringWithDefault(_T("RepeatFileHead"), _T(""));
    sReepatFileHead = sValues;

    UINT uiK = 4;
    __FetchIntWithDefault(_T("K"), 4);
    uiK = static_cast<UINT>(iValues);

    UINT uiMinHit = 2;
    __FetchIntWithDefault(_T("Hit"), 2);
    uiMinHit = static_cast<UINT>(iValues);

    UINT uiIteration = 100;
    __FetchIntWithDefault(_T("Iteration"), 100);
    uiIteration = static_cast<UINT>(iValues);

    QLQuantumKmeans::Kmeans2D(sLoadFile, sKFileHead, sCenterFileHead, sReepatFileHead, uiK, uiIteration, uiMinHit);
}

void TestKNN2D(CParameters& params)
{
    CCString sValues;
    INT iValues;

    CCString sTrainingPoints;
    __FetchStringWithDefault(_T("TrainingSet"), _T(""));
    sTrainingPoints = sValues;

    CCString sTestPoints;
    __FetchStringWithDefault(_T("TestSet"), _T(""));
    sTestPoints = sValues;

    CCString sTrainingClusters;
    __FetchStringWithDefault(_T("TrainingCluster"), _T(""));
    sTrainingClusters = sValues;

    CCString sTestingClusters;
    __FetchStringWithDefault(_T("TestSave"), _T(""));
    sTestingClusters = sValues;

    CCString sReepatFile;
    __FetchStringWithDefault(_T("RepeatSave"), _T(""));
    sReepatFile = sValues;

    //UINT uiPointHit = 2;
    //__FetchIntWithDefault(_T("PointHit"), 2);
    //uiPointHit = static_cast<UINT>(iValues);

    UINT uiKHit = 2;
    __FetchIntWithDefault(_T("KHit"), 2);
    uiKHit = static_cast<UINT>(iValues);

    UINT uiMaxCluster = 4;
    __FetchIntWithDefault(_T("MaxCluster"), 4);
    uiMaxCluster = static_cast<UINT>(iValues);

    QLQuantumKmeans::KNN2D(sTrainingPoints, sTestPoints, sTrainingClusters, sTestingClusters, sReepatFile, uiKHit, uiMaxCluster);
}

void TestKNN3D(CParameters& params)
{
    CCString sValues;
    INT iValues;

    CCString sTrainingPoints;
    __FetchStringWithDefault(_T("TrainingSet"), _T(""));
    sTrainingPoints = sValues;

    CCString sTestPoints;
    __FetchStringWithDefault(_T("TestSet"), _T(""));
    sTestPoints = sValues;

    CCString sTrainingClusters;
    __FetchStringWithDefault(_T("TrainingCluster"), _T(""));
    sTrainingClusters = sValues;

    CCString sTestingClusters;
    __FetchStringWithDefault(_T("TestSave"), _T(""));
    sTestingClusters = sValues;

    CCString sReepatFile;
    __FetchStringWithDefault(_T("RepeatSave"), _T(""));
    sReepatFile = sValues;

    //UINT uiPointHit = 2;
    //__FetchIntWithDefault(_T("PointHit"), 2);
    //uiPointHit = static_cast<UINT>(iValues);

    UINT uiKHit = 2;
    __FetchIntWithDefault(_T("KHit"), 2);
    uiKHit = static_cast<UINT>(iValues);

    UINT uiMaxCluster = 4;
    __FetchIntWithDefault(_T("MaxCluster"), 4);
    uiMaxCluster = static_cast<UINT>(iValues);

    QLQuantumKmeans::KNN3D(sTrainingPoints, sTestPoints, sTrainingClusters, sTestingClusters, sReepatFile, uiKHit, uiMaxCluster);
}


//=============================================================================
// END OF FILE
//=============================================================================