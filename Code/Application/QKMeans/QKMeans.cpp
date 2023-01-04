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
#if 0

    CParameters params;
    CYAMLParser::ParseFile(_T("../QKMeans.yaml"), params);

    params.Dump();

    INT iVaules;
    __FetchIntWithDefault(_T("DebugChange"), 0);
    UBOOL bDebugChange = (iVaules != 0);

    __FetchIntWithDefault(_T("K"), 0);
    BYTE k = static_cast<BYTE>(iVaules);

    __FetchIntWithDefault(_T("Start"), 0);
    UINT uiStart = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("End"), 100);
    UINT uiEnd = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("Stop"), 300);
    UINT uiStop = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("HasCoefficientIdx"), 0);
    UBOOL bHasC = (iVaules != 0);

    __FetchIntWithDefault(_T("CoefficientIdxStart"), 0);
    UINT uiCStart = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("CoefficientIdxEnd"), 20);
    UINT uiCEnd = static_cast<UINT>(iVaules);

    TArray<INT> energy;
    params.FetchValueArrayINT(_T("Energy"), energy);

    CCString sFileHead;
    params.FetchStringValue(_T("FileName"), sFileHead);

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
                    appGeneral(_T("%s saved.\n"), sFileToSave);
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
                appGeneral(_T("%s saved.\n"), sFileToSave);
            }
        }
    }
#endif

    QuantumKMeans(_T("../QKmeans.yaml"));

    return 0;
}


//=============================================================================
// END OF FILE
//=============================================================================