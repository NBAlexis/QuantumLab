//=============================================================================
// FILENAME : Fermion1D.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/01/2024 nbale]
//=============================================================================

#include "FermionSimulation.h"

void SaveWaveFunction(const CCString& sFile, Real* real, Real* imag, UINT len)
{
    std::ofstream file(sFile);
    CCString sOut;

    for (UINT i = 0; i < len; ++i)
    {
        sOut.Format(_T("%f, "), real[i]);
        file << sOut.c_str();
        sOut.Format(_T("%f"), imag[i]);
        file << sOut.c_str();
        file << std::endl;
    }

    file.flush();
    file.close();
}

void LoadWaveFunction(const CCString& sFile, Real** real, Real** imag, UINT& len)
{
    std::ifstream file(sFile);
    if (!file.good())
    {
        appWarning(_T("file not exist: %s\n"), sFile.c_str());
        *real = NULL;
        *imag = NULL;
        len = 0;
        return;
    }

    const SIZE_T buf_size = 1024;
    static TCHAR buf[buf_size];
    memset(buf, 0, sizeof(TCHAR) * buf_size);
    TArray<Real> reallst;
    TArray<Real> imaglst;
    while (file.getline(buf, buf_size))
    {
        CCString sLine(buf);
        TArray<CCString> sep = appGetStringList(sLine, _T(','), EGSLF_IgnorTabSpace | EGSLF_IgnorTabSpaceInSide);

        if (2 != sep.Num())
        {
            appCrucial(_T("CSV file format not good!\n"));
            file.close();
            *real = NULL;
            *imag = NULL;
            len = 0;
            return;
        }

#if _QL_DOUBLEFLOAT
        const Real realpart = appStoD(sep[0]);
        const Real imagpart = appStoD(sep[1]);
#else
        const Real realpart = appStoF(sep[0]);
        const Real imagpart = appStoF(sep[1]);
#endif
        reallst.AddItem(realpart);
        imaglst.AddItem(imagpart);
        memset(buf, 0, sizeof(TCHAR) * buf_size);
    }

    file.close();


    *real = (Real*)malloc(sizeof(Real) * reallst.Num());
    *imag = (Real*)malloc(sizeof(Real) * reallst.Num());
    len = reallst.Num();
    memcpy(*real, reallst.GetData(), sizeof(Real) * reallst.Num());
    memcpy(*imag, imaglst.GetData(), sizeof(Real) * reallst.Num());
}

void Fermion1DSimulation(CParameters& params)
{
    Real fValues = F(0.0);
    INT iValues = 0;

    __FetchRealWithDefault(_T("Spacing"), F(0.0));
    Real fSpacing = fValues;

    __FetchRealWithDefault(_T("Mass"), F(0.0));
    Real fMass = fValues;

    __FetchRealWithDefault(_T("G0"), F(0.0));
    Real fG0 = fValues;

    __FetchRealWithDefault(_T("G1"), F(0.0));
    Real fG1 = fValues;

    __FetchRealWithDefault(_T("F1sq"), F(0.0));
    Real fF1Sq = fValues;

    __FetchRealWithDefault(_T("Time"), F(0.0));
    Real fTime = fValues;

    __FetchIntWithDefault(_T("Site"), 8);
    UINT iSite = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("DecomposeStep"), 10);
    UINT iDecomposeStep = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("Step"), 10);
    INT iMaxStep = iValues;

    CCString sValues;
    __FetchStringWithDefault(_T("SaveFile"), _T(""));
    const CCString sFileName = sValues;


    SJordanWeigner1DTerms param;
    param.m_fLatticeSpacing = fSpacing;
    param.m_fMass = fMass;
    param.m_fG0 = fG0;
    param.m_fG1 = fG1;
    param.m_f1Sq = fF1Sq;

    CHamitonianFermion1D hamiltonian(iSite, param);
    
    QLGate simulationGate = hamiltonian.BuildSimulationCircuit(fTime, iDecomposeStep);

    //simulationGate.DebugPrint(4);

    UINT veclen = static_cast<UINT>(1U << simulationGate.m_lstQubits.Num());

    QLSimulatorParametersVector simulateParam;
    simulateParam.m_MasterGate = simulationGate;
    simulateParam.m_byQubitCount = static_cast<BYTE>(simulationGate.m_lstQubits.Num());
    simulateParam.m_bPrint = FALSE;

    QLSimulatorOutputVector output;
    output.m_bOutputToBuffer = TRUE;
    output.m_pRealBuffer = (Real*)malloc(sizeof(Real) * veclen);
    output.m_pImageBuffer = (Real*)malloc(sizeof(Real) * veclen);

    QLSimulatorVector simulator;
    CCString sFullFileName;
    for (INT i = 0; i < iMaxStep; ++i)
    {
        simulator.Simulate(&simulateParam, &output);

        simulateParam.m_bHasInitial = TRUE;
        simulateParam.m_pRealWavefunction = output.m_pRealBuffer;
        simulateParam.m_pImageWavefunction = output.m_pImageBuffer;

        sFullFileName.Format(sFileName.c_str(), i);

        //Note that, the initial of ancilla is |0>, so we use the first half of real and imag
        SaveWaveFunction(sFullFileName, output.m_pRealBuffer, output.m_pImageBuffer, veclen >> 1);
    }

    appSafeFree(output.m_pRealBuffer);
    appSafeFree(output.m_pImageBuffer);
}

void Fermion1DMeasure(CParameters& params)
{
    Real fValues = F(0.0);
    INT iValues = 0;

    __FetchRealWithDefault(_T("Spacing"), F(0.0));
    Real fSpacing = fValues;

    __FetchRealWithDefault(_T("Mass"), F(0.0));
    Real fMass = fValues;

    __FetchRealWithDefault(_T("G0"), F(0.0));
    Real fG0 = fValues;

    __FetchRealWithDefault(_T("G1"), F(0.0));
    Real fG1 = fValues;

    __FetchRealWithDefault(_T("F1sq"), F(0.0));
    Real fF1Sq = fValues;

    __FetchIntWithDefault(_T("Site"), 8);
    UINT iSite = static_cast<UINT>(iValues);

    __FetchIntWithDefault(_T("Repeat"), -1);
    INT iRepeat = iValues;

    CCString sValues;
    __FetchStringWithDefault(_T("SaveFile"), _T(""));
    const CCString sFileName = sValues;

    UINT len;
    Real* pReal = NULL;
    Real* pImag = NULL;
    LoadWaveFunction(sFileName, &pReal, &pImag, len);

    appGeneral(_T("%s file loaded\n"), sFileName.c_str());

    SJordanWeigner1DTerms param;
    param.m_fLatticeSpacing = fSpacing;
    param.m_fMass = fMass;
    param.m_fG0 = fG0;
    param.m_fG1 = fG1;
    param.m_f1Sq = fF1Sq;

    CHamitonianFermion1D hamiltonian(iSite, param);

    Real fRes = hamiltonian.Measure(pReal, pImag, iRepeat);

    appGeneral(_T("res = %f\n"), fRes);
}




//=============================================================================
// END OF FILE
//=============================================================================