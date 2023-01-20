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
    CParameters params;
    CYAMLParser::ParseFile(_T("../QKMeans.yaml"), params);
    INT iVaules;
    __FetchIntWithDefault(_T("QKMeans"), 0);
    UBOOL bQKMeans = (iVaules != 0);

    if (bQKMeans)
    {
        QuantumKMeans(_T("../QKMeans.yaml"));
    }
    else
    {
        ClassicalKMeans(_T("../QKMeans.yaml"));
    }

    return 0;
}


//=============================================================================
// END OF FILE
//=============================================================================