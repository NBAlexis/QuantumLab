//=============================================================================
// FILENAME : CHamitonianPauliNeighbour.h
// 
// DESCRIPTION:
// 
// \sum _{<m,n>} sigma_i (m)sigma _j(n)
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CHAMITONIANPAULINEIGHBOUR_H_
#define _CHAMITONIANPAULINEIGHBOUR_H_

__BEGIN_NAMESPACE

class QLAPI CHamitonianPauliNeighbour : public CHamitonianTerm
{
public:

    CHamitonianPauliNeighbour(Real fCoeff = F(1.0), EPauliType eType1 = EPauliType::EPT_Z, EPauliType eType2 = EPauliType::EPT_Z)
        : CHamitonianTerm(fCoeff)
        , m_eType1(eType1)
        , m_eType2(eType2)
    {

    }

    TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const override;

    EPauliType m_eType1;
    EPauliType m_eType2;
};

__END_NAMESPACE


#endif //#ifndef _CHAMITONIANPAULINEIGHBOUR_H_

//=============================================================================
// END OF FILE
//=============================================================================