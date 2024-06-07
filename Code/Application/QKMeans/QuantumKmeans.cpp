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

    BYTE uiK = 4;
    __FetchIntWithDefault(_T("K"), 4);
    uiK = static_cast<BYTE>(iValues);

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

void TestKNN2DAnsatz(CParameters& params)
{
    CCString sValues;
    INT iValues;

    CCString sAnsatz;
    __FetchStringWithDefault(_T("Ansatz"), _T(""));
    sAnsatz = sValues;

    CCString sTestPoints;
    __FetchStringWithDefault(_T("TestSet"), _T(""));
    sTestPoints = sValues;

    CCString sTestingClusters;
    __FetchStringWithDefault(_T("TestSave"), _T(""));
    sTestingClusters = sValues;

    CCString sReepatFile;
    __FetchStringWithDefault(_T("RepeatSave"), _T(""));
    sReepatFile = sValues;

    UINT uiKHit = 2;
    __FetchIntWithDefault(_T("KHit"), 2);
    uiKHit = static_cast<UINT>(iValues);

    QLQuantumKmeans::KNN2DAnsatz(sAnsatz, sTestPoints, sTestingClusters, sReepatFile, uiKHit);
}

void TestKNNAnsatz(CParameters& params)
{
    CCString sValues;
    INT iValues;

    CCString sAnsatz;
    __FetchStringWithDefault(_T("Ansatz"), _T(""));
    sAnsatz = sValues;

    CCString sTestPoints;
    __FetchStringWithDefault(_T("TestSet"), _T(""));
    sTestPoints = sValues;

    CCString sTestingClusters;
    __FetchStringWithDefault(_T("TestSave"), _T(""));
    sTestingClusters = sValues;

    BYTE ansatzQubit = 0;
    __FetchIntWithDefault(_T("AnsatzQubits"), 0);
    ansatzQubit = static_cast<BYTE>(iValues);

    UINT uiAnsatzLevel = 1;
    __FetchIntWithDefault(_T("AnsatzLevel"), 1);
    uiAnsatzLevel = static_cast<UINT>(iValues);

    BYTE uiMeasureBits = 1;
    __FetchIntWithDefault(_T("NumberOfMeasure"), 1);
    uiMeasureBits = static_cast<BYTE>(iValues);

    __FetchIntWithDefault(_T("IsAdaptive"), 0);
    UBOOL bIsAdap = (iValues != 0);

    UINT uiRepeat = 100;
    __FetchIntWithDefault(_T("Repeat"), 2);
    uiRepeat = static_cast<UINT>(iValues);

    QLQuantumKmeans::KNNAnsatz(sAnsatz, sTestPoints, sTestingClusters, ansatzQubit, uiMeasureBits, bIsAdap, uiAnsatzLevel, uiRepeat);
}

void TestKNNAE(CParameters& params)
{
    CCString sValues;
    INT iValues;

    CCString sPoints;
    __FetchStringWithDefault(_T("PointSet"), _T(""));
    sPoints = sValues;

    CCString sTestPoints;
    __FetchStringWithDefault(_T("TestSet"), _T(""));
    sTestPoints = sValues;

    CCString sTestingClusters;
    __FetchStringWithDefault(_T("TestSave"), _T(""));
    sTestingClusters = sValues;

    BYTE uiMeasureBits = 1;
    __FetchIntWithDefault(_T("NumberOfMeasure"), 1);
    uiMeasureBits = static_cast<BYTE>(iValues);

    UINT uiRepeat = 100;
    __FetchIntWithDefault(_T("Repeat"), 2);
    uiRepeat = static_cast<UINT>(iValues);

    QLQuantumKmeans::KNNAE(sPoints, sTestPoints, sTestingClusters, uiMeasureBits, uiRepeat);
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


void QAnomaly2D(CParameters& params)
{
    CCString sValues;
    Real fValues;

    CCString sPointFile;
    __FetchStringWithDefault(_T("PointFile"), _T(""));
    sPointFile = sValues;

    CCString sReferenceFile;
    __FetchStringWithDefault(_T("ReferenceFile"), _T(""));
    sReferenceFile = sValues;

    CCString sBuildRate;
    __FetchStringWithDefault(_T("BuildRateSave"), _T(""));
    sBuildRate = sValues;

    Real fMinX = F(-5.0);
    __FetchRealWithDefault(_T("MinX"), F(-5.0));
    fMinX = fValues;

    Real fMinY = F(-5.0);
    __FetchRealWithDefault(_T("MinY"), F(-5.0));
    fMinY = fValues;

    Real fMaxX = F(5.0);
    __FetchRealWithDefault(_T("MaxX"), F(5.0));
    fMaxX = fValues;

    Real fMaxY = F(5.0);
    __FetchRealWithDefault(_T("MaxY"), F(5.0));
    fMaxY = fValues;

    QLQuantumKmeans::QAnomaly2D(sReferenceFile, sPointFile, sBuildRate, fMinX, fMaxX, fMinY, fMaxY);
}

void QAnomaly3D(CParameters& params)
{
    CCString sValues;

    CCString sPointFile;
    __FetchStringWithDefault(_T("PointFile"), _T(""));
    sPointFile = sValues;

    CCString sReferenceFile;
    __FetchStringWithDefault(_T("ReferenceFile"), _T(""));
    sReferenceFile = sValues;

    CCString sBuildRate;
    __FetchStringWithDefault(_T("BuildRateSave"), _T(""));
    sBuildRate = sValues;

    QLQuantumKmeans::QAnomaly3D(sReferenceFile, sPointFile, sBuildRate);
}

//=============================================================================
// END OF FILE
//=============================================================================