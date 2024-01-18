//=============================================================================
// FILENAME : CHamitonianTerm.h
// 
// DESCRIPTION:
// This is one term in Hamitonian which only touch 1-site once
//
// REVISION: [dd/mm/yy]
//  [04/12/2023 nbale]
//=============================================================================

#ifndef _CHAMITONIANTERM_H_
#define _CHAMITONIANTERM_H_

__BEGIN_NAMESPACE

class QLAPI CHamitonianTerm
{
public:

    CHamitonianTerm(Real fCoeff) : m_fCoefficient(fCoeff) {}
    virtual ~CHamitonianTerm() {}

    virtual TArray<PauliProduct> GetAllTerms(const CLattice* pLattice) const = 0;

    /**
    * build a circuit to implement Controlled-Exp[i c t H]
    * where c is the coefficient
    * 
    */
    virtual QLGate BuildCircuit(const CLattice* pLattice, Real fTrotterTime) const;

    /**
    * When H = sigma1 otimes sigma2 ..... otimes sigman, where n is the number of sites
    * Use Pauli product to do one step simulation
    * 
    * The ancilla is the last qubit
    * The number of qubits is n+1
    */
    inline QLGate OneTermOneStep(const TArray<BYTE>& pauliType, Real fTrotterTime) const
    {
        return OneTermOneStep(pauliType, fTrotterTime, m_fCoefficient);
    }

    inline static QLGate OneTermOneStep(const TArray<BYTE>& pauliType, Real fTrotterTime, Real fCoefficient)
    {
        PauliProduct prod(pauliType, fCoefficient);
        return prod.OneStepGate(fTrotterTime);
    }

    /**
    * when iRepeat <= 0, measure exact
    */
    inline static Real OneTermMeasure(const TArray<BYTE>& pauliType, const Real* hostWaveFunctionReal, const Real* hostWaveFunctionImagin, INT iRepeat);

    virtual Real Measure(const CLattice* pLattice, const Real* hostWaveFunctionReal, const Real* hostWaveFunctionImagin, INT iRepeat) const;

    Real m_fCoefficient;
};


__END_NAMESPACE


#endif //#ifndef _CHAMITONIANTERM_H_

//=============================================================================
// END OF FILE
//=============================================================================