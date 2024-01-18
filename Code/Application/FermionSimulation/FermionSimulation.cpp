//=============================================================================
// FILENAME : FermionSimulation.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [07/01/2024 nbale]
//=============================================================================

#include "FermionSimulation.h"
#include <random>

int main()
{
    QLRandomInitializer random(ERandom::ER_XORWOW, appGetTimeStamp());
    CParameters params;
    CYAMLParser::ParseFile(_T("../FermionSimulation.yaml"), params);
    params.Dump();

    CCString sValues;
    __FetchStringWithDefault(_T("JobType"), _T("EFJ_Simulation"));
    EFermionJob eJob = __STRING_TO_ENUM(EFermionJob, sValues);

    switch (eJob)
    {
    case EFJ_Simulation:
    {
        CParameters jobparam;
        if (params.FetchParameterValue(_T("JobSimulation"), jobparam))
        {
            Fermion1DSimulation(jobparam);
        }
        else
        {
            appCrucial(_T("JobSimulation not found!\n"));
        }
    }
    break;
    case EFJ_Measure:
    {
        CParameters jobparam;
        if (params.FetchParameterValue(_T("JobMeasure"), jobparam))
        {
            Fermion1DMeasure(jobparam);
        }
        else
        {
            appCrucial(_T("JobMeasure not found!\n"));
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