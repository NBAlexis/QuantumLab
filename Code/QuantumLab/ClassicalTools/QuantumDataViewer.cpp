//=============================================================================
// FILENAME : QuantumDataViewer.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [08/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

/**
* upLimit: the full list of upper
* lowerLimit: the full list of lower
* for example, at some qubit, 
*    if it want index = 1, it should be lower=1, upper=2, 
*    if it want index = 0, it should be lower=0, upper=1, 
*    if it want all, it should be lower=0, upper=2, 
* addedIndex: start position in stateData
* qubitIndex: the qubit this recursive is reading
*/
TArray<QLComplex> RecursivelyGetElement(const QLComplex* stateData, TArray<BYTE> upLimit, TArray<BYTE> lowerLimit, UINT addedIndex, BYTE qubitIndex)
{
    TArray<QLComplex> ret;
    //factor is stride
    UINT factor = 1U << (upLimit.Num() - 1 - qubitIndex);

    if (1 == factor)
    {
        //this is qubit 0
        for (BYTE i = lowerLimit[qubitIndex]; i < upLimit[qubitIndex]; ++i)
        {
            ret.AddItem(stateData[addedIndex + i]);
        }
        return ret;
    }

    for (BYTE i = lowerLimit[qubitIndex]; i < upLimit[qubitIndex]; ++i)
    {
        ret.Append(RecursivelyGetElement(stateData, upLimit, lowerLimit, addedIndex + i * factor, qubitIndex + 1));
    }
    return ret;
}

QLMatrix QLAPI ShowStateVectorDetail(const QLComplex* statedata, TArray<BYTE> index, UBOOL bNormalize)
{
    TArray<BYTE> indexUpperLimit;
    TArray<BYTE> indexLowerLimit;
    INT expectedCount = 1;
    for (INT i = 0; i < index.Num(); ++i)
    {
        if (2 == index[i])
        {
            indexUpperLimit.InsertAt(0, 2);
            indexLowerLimit.InsertAt(0, static_cast<BYTE>(0));
            expectedCount = (expectedCount << 1);
        }
        else
        {
            indexUpperLimit.InsertAt(0, index[i] + 1);
            indexLowerLimit.InsertAt(0, index[i]);
        }
    }

    TArray<QLComplex> res = RecursivelyGetElement(statedata, indexUpperLimit, indexLowerLimit, 0, 0);
    if (expectedCount != res.Num())
    {
        appCrucial(_T("bad thing happens in ShowStateVectorDetail!\n"));
        return QLMatrix();
    }
    if (bNormalize)
    {
        UINT len;
        res = NormalizeV(res, len);
    }
    QLMatrix ret = QLMatrix::CopyCreate(static_cast<UINT>(expectedCount), 1, res.GetData());
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================