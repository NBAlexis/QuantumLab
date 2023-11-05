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

void TestProbabilityDifferentDimension(CParameters& params)
{
    //build circuit

    //update measurement propability
    CCString sValues;
    INT iValues;

    CCString sFileName;
    __FetchStringWithDefault(_T("CSVFileName"), _T(""));
    sFileName = sValues;

    CCString sAmplitudeFileName;
    __FetchStringWithDefault(_T("AmplitudeSave"), _T(""));
    sAmplitudeFileName = sValues;

    CCString sMeasureFileName;
    __FetchStringWithDefault(_T("MeasureSave"), _T(""));
    sMeasureFileName = sValues;

    UINT uiCount = 0;
    UINT uiRepeat = 0;
    __FetchIntWithDefault(_T("Count"), 0);
    uiCount = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("Repeat"), 0);
    uiRepeat = static_cast<UINT>(iValues);

    QLQuantumKmeans::TestCircuitBuildState(sFileName, sAmplitudeFileName, sMeasureFileName, uiCount, uiRepeat);

}

//=============================================================================
// END OF FILE
//=============================================================================