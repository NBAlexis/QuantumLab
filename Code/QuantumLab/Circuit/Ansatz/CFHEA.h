//=============================================================================
// FILENAME : CFHEA.h
// 
// DESCRIPTION:
// Fixed-structure layered Hardware-Efficient Ansatz
// in
// 1909.05820
// 
// it is
// a column of Ry
// then:
// 
// layers with 4 columns:
// 1 - cz for 01, 23, 45, ...
// 2 - ry on each qubit
// 3 - cz for 12, 34, 56, ...
// 4 - ry on qubits except for the first and last
// 
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#ifndef _CFHEA_H_
#define _CFHEA_H_

__BEGIN_NAMESPACE

class QLAPI CFHEA : public CAnsatz
{
public:

    CFHEA(BYTE qubits, UINT uiLayerCount, UBOOL bOnlyReal = FALSE);
    QLGate BuildState(const TArray<Real>& params) const override
    {
        if (m_bOnlyReal)
        {
            return BuildStateReal(params);
        }
        return BuildStateCmp(params);
    }

protected:

    QLGate BuildStateReal(const TArray<Real>& params) const;
    QLGate BuildStateCmp(const TArray<Real>& params) const;

    UINT m_uiLayerCount;
    UBOOL m_bOnlyReal;
};

__END_NAMESPACE


#endif //#ifndef _CFHEA_H_

//=============================================================================
// END OF FILE
//=============================================================================