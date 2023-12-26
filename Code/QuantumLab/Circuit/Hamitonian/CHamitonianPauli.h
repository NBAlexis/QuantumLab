//=============================================================================
// FILENAME : CHamitonianPauli.h
// 
// DESCRIPTION:
// 
// sigma_i(n) term
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CHAMITONIANPAULI_H_
#define _CHAMITONIANPAULI_H_

__BEGIN_NAMESPACE

enum class EPauliType : BYTE
{
    EPT_I = 0,
    EPT_X = 1,
    EPT_Y = 2,
    EPT_Z = 3,
};

class QLAPI CHamitonianPauli : public CHamitonianTerm
{
public:

    CHamitonianPauli(Real fCoeff = F(1.0), EPauliType eType = EPauliType::EPT_X)
        : CHamitonianTerm(fCoeff)
        , m_eType(eType)
    {
        
    }

    TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const override;

    EPauliType m_eType;

};


__END_NAMESPACE


#endif //#ifndef _CHAMITONIANPAULI_H_

//=============================================================================
// END OF FILE
//=============================================================================