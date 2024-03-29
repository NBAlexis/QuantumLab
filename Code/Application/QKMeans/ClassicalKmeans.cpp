//=============================================================================
// FILENAME : ClassicalKmeans.cpp
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [21/12/2022 nbale]
//=============================================================================

#include "QKMeans.h"

void ClassicalKMeans(CParameters& params)
{
    INT iValues;
    __FetchIntWithDefault(_T("DebugChange"), 0);
    UBOOL bDebugChange = (iValues != 0);

    __FetchIntWithDefault(_T("CK"), 0);
    BYTE k = static_cast<BYTE>(iValues);

    __FetchIntWithDefault(_T("Start"), 0);
    UINT uiStart = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("End"), 100);
    UINT uiEnd = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("Stop"), 300);
    UINT uiStop = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("HasCoefficientIdx"), 0);
    UBOOL bHasC = (iValues != 0);

    __FetchIntWithDefault(_T("CoefficientIdxStart"), 0);
    UINT uiCStart = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("CoefficientIdxEnd"), 20);
    UINT uiCEnd = static_cast<UINT>(iValues);

    TArray<INT> energy;
    params.FetchValueArrayINT(_T("Energy"), energy);

    CCString sFileHead;
    params.FetchStringValue(_T("CFileName"), sFileHead);

    QLRandomInitializer random;

    for (INT i = 0; i < energy.Num(); ++i)
    {
        if (bHasC)
        {
            for (UINT c = uiCStart; c <= uiCEnd; ++c)
            {
                CCString sFileToLoad;
                sFileToLoad.Format(_T("%s-%d-%d.csv"), sFileHead.c_str(), energy[i], c);
                QLGPUKmeans kmeans(sFileToLoad, 12, k, bDebugChange);

                for (UINT j = uiStart; j <= uiEnd; ++j)
                {
                    CCString sFileToSave;
                    sFileToSave.Format(_T("%s-%d-%d-%d-%d.csv"), sFileHead.c_str(), energy[i], k, c, j);
                    kmeans.Build(uiStop);
                    kmeans.Save(sFileToSave);
                    appGeneral(_T("%s saved.\n"), sFileToSave.c_str());
                }
            }
        }
        else
        {
            CCString sFileToLoad;
            sFileToLoad.Format(_T("%s-%d.csv"), sFileHead.c_str(), energy[i]);
            QLGPUKmeans kmeans(sFileToLoad, 12, k, bDebugChange);

            for (UINT j = uiStart; j <= uiEnd; ++j)
            {
                CCString sFileToSave;
                sFileToSave.Format(_T("%s-%d-%d-%d.csv"), sFileHead.c_str(), energy[i], k, j);
                kmeans.Build(uiStop);
                kmeans.Save(sFileToSave);
                appGeneral(_T("%s saved.\n"), sFileToSave.c_str());
            }
        }
    }
}


//=============================================================================
// END OF FILE
//=============================================================================