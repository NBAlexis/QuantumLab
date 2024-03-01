//=============================================================================
// FILENAME : QKMeans.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/12/2022 nbale]
//=============================================================================

#include "QKMeans.h"

int main()
{
    QLRandomInitializer random(ERandom::ER_XORWOW, appGetTimeStamp());
    CParameters params;
    CYAMLParser::ParseFile(_T("../QKMeans.yaml"), params);
    params.Dump();


    CCString sValues;
    __FetchStringWithDefault(_T("KMeansJob"), _T("EQKJ_TestDifferentDimension"));
    EQKmeansJob eJob = __STRING_TO_ENUM(EQKmeansJob, sValues);
    switch (eJob)
    {
    case EQKJ_Kmeans:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("GPUKmeans"), jobparam))
            {
                ClassicalKMeans(jobparam);
            }
            else
            {
                appCrucial(_T("GPUKmeans not found!\n"));
            }
        }
        break;
    case EQKJ_TestBuildStateOnce:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("CircuitTestDim"), jobparam))
            {
                TestCircuitBuildStateOnce(jobparam);
            }
            else
            {
                appCrucial(_T("CircuitTestDim not found!\n"));
            }
        }
        break;
    case EQKJ_TestBuildStateRandom:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("CircuitTestDim"), jobparam))
            {
                TestProbabilityToBuildStateRandom(jobparam);
            }
            else
            {
                appCrucial(_T("CircuitTestDim not found!\n"));
            }
        }
        break;
    case EQKJ_TestBuildStateFromFile:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("CircuitTestDim"), jobparam))
            {
                TestProbabilityToBuildStateFromFile(jobparam);
            }
            else
            {
                appCrucial(_T("CircuitTestDim not found!\n"));
            }
        }
        break;
    case EQKJ_TestQKmeans2D:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QKMeans2D"), jobparam))
            {
                TestKMeans2D(jobparam);
            }
            else
            {
                appCrucial(_T("QKMeans2D not found!\n"));
            }
        }
        break;
    case EQKJ_TestQKNN2D:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QKNN2D"), jobparam))
            {
                TestKNN2D(jobparam);
            }
            else
            {
                appCrucial(_T("QKNN2D not found!\n"));
            }
        }
        break;
    case EQKJ_TestQKNN3D:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QKNN3D"), jobparam))
            {
                TestKNN3D(jobparam);
            }
            else
            {
                appCrucial(_T("QKNN3D not found!\n"));
            }
        }
        break;
    case EQKJ_TestFitPointSet:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("FitWaveFunction"), jobparam))
            {
                TestFitPointSet(jobparam);
            }
            else
            {
                appCrucial(_T("FitWaveFunction not found!\n"));
            }
        }
        break;
    case EQKJ_TestFitPointSetAdap:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("FitWaveFunction"), jobparam))
            {
                TestFitPointSetAdap(jobparam);
            }
            else
            {
                appCrucial(_T("FitWaveFunction not found!\n"));
            }
        }
        break;
    case EQKJ_QAnomaly2D:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QAnomaly2D"), jobparam))
            {
                QAnomaly2D(jobparam);
            }
            else
            {
                appCrucial(_T("QAnomaly2D not found!\n"));
            }
        }
        break;
    case EQKJ_QAnomaly3D:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QAnomaly3D"), jobparam))
            {
                QAnomaly3D(jobparam);
            }
            else
            {
                appCrucial(_T("QAnomaly3D not found!\n"));
            }
        }
        break;
    case EQKJ_TestQKNNAnsatz2D:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QKNN2DAnsatz"), jobparam))
            {
                TestKNN2DAnsatz(jobparam);
            }
            else
            {
                appCrucial(_T("QKNN2DAnsatz not found!\n"));
            }
        }
        break;
    case EQKJ_TestQKNNAnsatz:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QKNNAnsatz"), jobparam))
            {
                TestKNNAnsatz(jobparam);
            }
            else
            {
                appCrucial(_T("QKNN2DAnsatz not found!\n"));
            }
        }
        break;
    case EQKJ_TestQKNNAE:
        {
            CParameters jobparam;
            if (params.FetchParameterValue(_T("QKNNAnsatz"), jobparam))
            {
                TestKNNAE(jobparam);
            }
            else
            {
                appCrucial(_T("QKNN2DAnsatz not found!\n"));
            }
        }
        break;
    default:
        break;
    }

    return 0;
}


//=============================================================================
// END OF FILE
//=============================================================================