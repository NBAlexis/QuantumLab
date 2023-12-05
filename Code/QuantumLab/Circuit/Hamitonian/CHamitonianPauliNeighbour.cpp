//=============================================================================
// FILENAME : CHamitonianPauliNeighbour.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [05/12/2023 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLGate CHamitonianPauliNeighbour::BuildCircuit(const CLattice* pLattice, Real fTrotterTime) const
{
    TArray<CLatticeSiteData> pairs = pLattice->GetNeighbourPairs();
    UINT uiControllerCount = pLattice->GetControllerCount();

    TArray<BYTE> toAdd;
    toAdd.Append(ByteSequnce, uiControllerCount + 1);
    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(uiControllerCount + 1));

    for (INT i = 0; i < pairs.Num(); ++i)
    {
        TArray<BYTE> pauliProduct;
        for (UINT uiQubit = 0; uiQubit < uiControllerCount; ++uiQubit)
        {
            if (static_cast<UINT>(pairs[i].m_lstSites[0]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(m_eType1));
            }
            else if (static_cast<UINT>(pairs[i].m_lstSites[1]) == uiQubit)
            {
                pauliProduct.AddItem(static_cast<BYTE>(m_eType2));
            }
            else
            {
                pauliProduct.AddItem(static_cast<BYTE>(EPauliType::EPT_I));
            }
        }

        ret.AppendGate(OneTermOneStep(pauliProduct, fTrotterTime), toAdd);
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================