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
    CParameters params;
    CYAMLParser::ParseFile(yamlFile, params);

    params.Dump();

    INT iVaules;
    __FetchIntWithDefault(_T("K"), 0);
    BYTE k = static_cast<BYTE>(iVaules);

    __FetchIntWithDefault(_T("Count"), 0);
    UINT uiCount = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("Measure"), 10);
    UINT uiMeasure = static_cast<UINT>(iVaules);

    //__FetchIntWithDefault(_T("Stop"), 300);
    //UINT uiStop = static_cast<UINT>(iVaules);

    //__FetchIntWithDefault(_T("HasCoefficientIdx"), 0);
    //UBOOL bHasC = (iVaules != 0);

    //__FetchIntWithDefault(_T("CoefficientIdxStart"), 0);
    //UINT uiCStart = static_cast<UINT>(iVaules);

    //__FetchIntWithDefault(_T("CoefficientIdxEnd"), 20);
    //UINT uiCEnd = static_cast<UINT>(iVaules);

    //TArray<INT> energy;
    //params.FetchValueArrayINT(_T("Energy"), energy);

    CCString sFileHead;
    params.FetchStringValue(_T("FileName"), sFileHead);


    Real testvectorlist[17][7] = {
        //centeroids
        {0.7984549572443157, 0.15824669704747435, 0.22671818545994732, 0.2370551931947057, -0.3711218641271573, 0.22654828653545345, 0.20193018526546913},
        {0.772146218399182, 0.4556218177339974, 0.3560005244611919, 0.19263353668288288, 0.05538862206650115, -0.11901347046259368, 0.12297486791631627},
        {0.7339442776271119, 0.0947055795599421, 0.1343762423609536, 0.08636859514849261, -0.4470209946934188, 0.297687799005556, 0.3720139062268942},
        {0.7478277640658433, 0.14389911808408826, 0.09011404398246027, 0.08905156026810167, -0.37511040645542587, 0.4259626944119037, 0.28608377048697325},
        {0.7593812738040595, 0.1119190850920762, 0.19647411098632378, 0.125438530702432, -0.423025994927114, 0.1897440676145073, 0.376196085221859},
        {0.7161528642498367, 0.10662050108995731, 0.07247424289346104, 0.07898324419704841, 0.41078528293370553, -0.4239624956359024, -0.34026097138294226},
        {0.7768483635124902, 0.20534682441745902, 0.10148396852415725, 0.2782590127729307, 0.30535319675130884, -0.4041683085900973, -0.10009810846379419},
        {0.7914085955068865, 0.333848826579712, 0.3455273591219377, 0.1592593444387407, -0.19036014750592278, 0.11016407649455838, 0.26285279482621676},
        {0.783438088907364, 0.34756714397184985, 0.19784078529492766, 0.3688134492921155, 0.14477935883008058, -0.26156386328320486, 0.02967806019070182},
        {0.7745322386205985, 0.288832977560113, 0.3657389702588672, 0.12767216787414543, -0.19208399517790806, 0.04564576659604052, 0.3572538990462893},
        {0.7815686156757491, 0.2893541717903871, 0.31148005155550107, 0.10972104797348028, 0.22446650642207475, -0.11400480873990233, -0.36466947013922496},
        {0.7626971941498496, 0.42449929008380777, 0.4215444778980757, 0.19774621751375338, -0.007728664900718448, 0.05654128532630878, -0.13428842153898918},
        {0.7932737634873663, 0.2070432349173152, 0.1589045700145796, 0.3243836521241661, 0.3077902202330975, -0.2826169994717652, -0.15088814015379282},
        {0.7941964816042992, 0.17034758634391103, 0.17480464089809844, 0.2805995886782965, 0.33864537975607245, -0.34012787985506726, -0.023941363857886037},
        {0.7884063300241814, 0.21825952493437248, 0.37196623341506146, 0.29366036131298445, 0.2297941649910477, -0.02081392711263835, -0.23009644525715744},
        {0.7589005902029682, 0.17044186522994276, 0.0770059485075803, 0.15956227079865187, 0.3312854043099351, -0.5031004979318865, 0.02773626273384198},
        {0.8054699093415323, 0.20950054846425206, 0.19970119940095096, 0.22360755838515164, 0.3236915075673832, -0.2317605419145365, -0.24281205750458354}
    };

    QLRandomInitializer random;

    QLQuantumKmeans qkmeans(k);
    //qkmeans.TestCircuit(reinterpret_cast<Real*>(testvectorlist));

    CCString sLoadFileName;
    sLoadFileName.Format(_T("%s.csv"), sFileHead);
    qkmeans.Prepare(sLoadFileName, uiCount);
    qkmeans.KMeans(sFileHead, 1, uiMeasure, TRUE);
    //QLGate gate = qkmeans.CompareCircuit((Real*)testvectorlist);

    //QLSimulatorParametersVector param;
    //param.m_byQubitCount = 5;
    //param.m_MasterGate = gate;
    //param.m_bPrint = TRUE;
    //QLSimulatorOutputVector out;
    //QLSimulatorVector sim;
    //sim.Simulate(&param, &out);
}


//=============================================================================
// END OF FILE
//=============================================================================