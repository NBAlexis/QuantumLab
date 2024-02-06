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
    EQKJ_TestBuildStateOnce,
    EQKJ_TestBuildStateFromFile,
    EQKJ_TestBuildStateRandom,
    EQKJ_TestQKmeans2D,
    EQKJ_TestQKNN2D,
    EQKJ_TestQKNN3D,
    EQKJ_TestFitPointSet,
    EQKJ_QAnomaly2D,
    EQKJ_QAnomaly3D,
    EQKJ_TestQKNNAnsatz2D,
    EQKJ_TestQKNNAnsatz,
    EQKJ_TestQKNNAE,
    EQKJ_TestAnomalyDetection
    )


extern void ClassicalKMeans(CParameters& yamlFile);

/**
* 
*/
extern void TestProbabilityDifferentDimension(CParameters& params);

/**
* randomly choose the vectors from a file
* for example choose 32 vectors as vi, and 1 vector as u
* to build \sum _i <vi|u>|i>
* and check the probability to build the state
*/
extern void TestProbabilityToBuildStateFromFile(CParameters& params);

extern void TestProbabilityToBuildStateRandom(CParameters& params);

extern void QAnomaly2D(CParameters& params);

extern void QAnomaly3D(CParameters& params);

extern void TestKNN2DAnsatz(CParameters& params);

extern void TestKNNAnsatz(CParameters& params);

extern void TestKNNAE(CParameters& params);


/**
* test whether the circuit can build the correct state
*/
extern void TestCircuitBuildStateOnce(CParameters& params);

/**
* test KMeans 2D
*/
extern void TestKMeans2D(CParameters& params);

/**
* test KNN
*/
extern void TestKNN2D(CParameters& params);
extern void TestKNN3D(CParameters& params);

/**
* Fit point set
*/
extern void TestFitPointSet(CParameters& params);

#endif //#ifndef _QL_QKMEANS_H_

//=============================================================================
// END OF FILE
//=============================================================================