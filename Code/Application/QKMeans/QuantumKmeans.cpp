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

void QuantumKMeans(const CCString& yamlFile)
{
#if 0

    CParameters params;
    CYAMLParser::ParseFile(yamlFile, params);

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

    Real testvectorlist[9][7] = {
        //centeroids
        {0.9086684979650952, 0.09587600692175047, 0.0562573312587942, -0.31005085705469415, -0.15293267084482415, 0.04402456229032295, 0.20126193610640858},
        {0.8925470708509142, -0.021797211311818256, 0.3258217801818233, 0.08912894145530632, 0.0034999919895174537, -0.2553038415268419, -0.1535854999874524},
        {0.9015507038060578, 0.07125583500775895, -0.007990222824483216, -0.35121702415010936, -0.08645558718659387, 0.05401878012506663, 0.21981604938795143},
        {0.8905231008619802, 0.08151762673553516, 0.0675498273157747, -0.3412510941203642, -0.0541644534357144, -0.03822915258476314, 0.2737022813761226},
        {0.9145060717730416, -0.05878468469975029, -0.05904610045470057, -0.33314892328569107, 0.12533873820948094, 0.07510097189512603, 0.15619988218347758},
        {0.8961016607770689, -0.08215813401151334, 0.0032743030710981813, 0.3514906953171573, 0.0522045780250354, -0.020136920735897144, -0.2521202313771924},
        {0.8772951369798326, -0.08306440997213464, 0.14354702673486164, -0.2884249264558623, 0.13105589224280814, -0.11947064515513744, 0.29700164544440777},
        {0.9036870650924249, 0.08298068982105113, 0.12234409879145151, 0.3071848703857185, -0.07201038243325047, 0.010626733448849443, -0.24866613713270744},

        //v
        {0.9149709636292653, 0.020665636751607288, -0.10320308230145858, 0.3151349105483805, -0.052553362395221785, -0.061261278369092356, -0.21430207419357325}
    };

    QLQuantumKmeans qkmeans(8, NULL);

    QLGate gate = qkmeans.CompareCircuit((Real*)testvectorlist);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 5;
    param.m_MasterGate = gate;
    param.m_bPrint = TRUE;
    QLSimulatorOutputVector out;
    QLSimulatorVector sim;
    sim.Simulate(&param, &out);
}


//=============================================================================
// END OF FILE
//=============================================================================