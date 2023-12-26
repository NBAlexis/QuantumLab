//=============================================================================
// FILENAME : PauliSimulate.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [04/10/2022 nbale]
//=============================================================================

#ifndef _PAULISIMULATE_H_
#define _PAULISIMULATE_H_

__BEGIN_NAMESPACE

static inline QLMatrix GetPauliProductMatrix(TArray<BYTE> lstPauli)
{
    if (lstPauli.Num() < 1)
    {
        return _I2;
    }

    const QLMatrix paulis[4] = { _I2, _PauliX, _PauliY, _PauliZ };

    QLMatrix ret = paulis[lstPauli[0]];

    for (INT i = 1; i < lstPauli.Num(); ++i)
    {
        ret = paulis[lstPauli[i]].KroneckerProduct(ret);
    }

    return ret;
}

class QLAPI PauliProduct
{
public:

    PauliProduct();

    PauliProduct(const TArray<BYTE>& pauliType, Real fCoefficient);

    QLMatrix GetMatrix() const 
    { 
        return GetPauliProductMatrix(m_lstPauliType) * m_fCoefficient;
    }

    /**
    * ancilla is the last qubit
    */
    QLGate OneStepGate(Real fTrotterTime) const;

    /**
    * project to this base
    */
    QLGate Project() const;

    UINT m_iOrder;
    TArray<BYTE> m_lstPauliType;
    Real m_fCoefficient;
};

extern QLMatrix QLAPI PrintPauliDecompsedMatrix(const TArray<PauliProduct>& lstPauli);

/**
* h must be Hermitian nxn matrix
* 
*/
extern TArray<PauliProduct> QLAPI DecomposePauliProducts(const QLMatrix& h, Real fMinimalKept = F(0.000001));

/**
* h must be Hermitian nxn matrix
* This gate is exp(i t h)
* where the last qubit is ancilla
*/
extern QLGate QLAPI PauliSimulateGate(const QLMatrix& h, Real t, UINT trotterStep, Real fMinimalKept = F(0.000001));

/**
* similar as above but use second order trotter
*/
extern QLGate QLAPI PauliSimulateGateLeapfrog(const QLMatrix& h, Real t, UINT trotterStep, Real fMinimalKept = F(0.000001));

extern QLGate QLAPI PauliSimulateGateLeapfrog(BYTE totalOrder, const TArray<PauliProduct>& h, Real t, UINT trotterStep);

__END_NAMESPACE


#endif //#ifndef _PAULISIMULATE_H_

//=============================================================================
// END OF FILE
//=============================================================================