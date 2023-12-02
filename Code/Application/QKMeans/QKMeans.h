//=============================================================================
// FILENAME : QKMeans.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/12/2022 nbale]
//=============================================================================

#ifndef _QL_QKMEANS_H_
#define _QL_QKMEANS_H_

#include "QuantumLab.h"

__DEFINE_ENUM(EQKmeansJob,
    EQKJ_Kmeans,
    EQKJ_TestDifferentDimension,
    EQKJ_TestDifferentVectorNumber,
    EQKJ_TestQKmeans,
    EQKJ_TestAnomalyDetection
    )


extern void ClassicalKMeans(CParameters& yamlFile);
extern void QuantumKMeans(CParameters& yamlFile);

extern void TestDistance();
extern void TestProbalityToBuild();
extern void TestProbalityMeanP6();

extern void TestProbabilityDifferentDimension(CParameters& params);

#endif //#ifndef _QL_QKMEANS_H_

//=============================================================================
// END OF FILE
//=============================================================================