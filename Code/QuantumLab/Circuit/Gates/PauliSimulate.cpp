//=============================================================================
// FILENAME : PauliSimulate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [04/10/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

PauliProduct::PauliProduct()
    : m_iOrder(0)
    , m_fCoefficient(F(0.0))
{

}

PauliProduct::PauliProduct(const TArray<BYTE>& pauliType, Real fCoefficient)
    : m_iOrder(static_cast<UINT>(pauliType.Num()))
    , m_lstPauliType(pauliType)
    , m_fCoefficient(fCoefficient)
{

}

QLGate PauliProduct::OneStepGate(Real fTrotterTime) const
{
    QLGate ret;
    ret.AddQubits(m_iOrder + 1);
    ret.m_sName = _T("Trotter");

    BYTE byAncilla = static_cast<BYTE>(m_iOrder);
    TArray<BYTE> abit;
    abit.AddItem(byAncilla);

    QLGate prepare;
    prepare.AddQubits(m_iOrder + 1);
    prepare.m_sName = _T("Prepare");

    for (INT i = 0; i < m_lstPauliType.Num(); ++i)
    {
        TArray<BYTE> bit;
        bit.AddItem(static_cast<BYTE>(i));
        TArray<BYTE> cbit;
        cbit.AddItem(static_cast<BYTE>(i));
        cbit.AddItem(byAncilla);
        QLGate cnot(EBasicOperation::EBO_CX);
        QLGate h(EBasicOperation::EBO_H);

        switch (m_lstPauliType[i])
        {
        case 1:
            {
                prepare.AppendGate(h, bit);
                prepare.AppendGate(cnot, cbit);
            }
            break;
        case 2:
            {
                QLGate p(EBasicOperation::EBO_P, -PI / 2);
                prepare.AppendGate(p, bit);
                prepare.AppendGate(h, bit);
                prepare.AppendGate(cnot, cbit);
            }
            break;
        case 3:
            {
                prepare.AppendGate(cnot, cbit);
            }
            break;
        }
    }

    ret.AppendGate(prepare, ret.m_lstQubits);
    QLGate rz(EBasicOperation::EBO_RZ, -fTrotterTime * m_fCoefficient * F(2.0));
    ret.AppendGate(rz, abit);
    prepare.Dagger();
    ret.AppendGate(prepare, ret.m_lstQubits);

    return ret;
}

QLMatrix QLAPI PrintPauliDecompsedMatrix(const TArray<PauliProduct>& lstPauli)
{
    QLMatrix ret = lstPauli[0].GetMatrix();
    for (INT i = 1; i < lstPauli.Num(); ++i)
    {
        ret = ret + lstPauli[i].GetMatrix();
    }
    return ret;
}

TArray<PauliProduct> QLAPI DecomposePauliProducts(const QLMatrix& h, Real fMinimalKept)
{
    UINT n = Log2(h.X());
    TArray<PauliProduct> ret;
    Real norm = F(1.0) / (1U << n);

    for (UINT idx = 0; idx < (1U << (2U * n)); ++idx)
    {
        TArray<BYTE> lstId;
        for (UINT j = 0; j < n; ++j)
        {
            lstId.AddItem(3U & (idx >> (2U * j)));
        }
        QLMatrix mtr = GetPauliProductMatrix(lstId);

        Real coefficient = mtr.VectorDot(h).x * norm;

        if (0 == idx || abs(coefficient) > fMinimalKept)
        {
            ret.AddItem(PauliProduct(lstId, coefficient));
        }
    }

    return ret;
}

QLGate QLAPI PauliSimulateGate(const QLMatrix& h, Real t, UINT trotterStep, Real fMinimalKept)
{
    QLGate ret;
    BYTE totalOrder = static_cast<BYTE>(Log2(h.X()));
    ret.AddQubits(totalOrder + 1);
    ret.m_sName = _T("Exp(iHt)");

    Real fOneTrotterStep = t / trotterStep;
    TArray<PauliProduct> decomposed = DecomposePauliProducts(h, fMinimalKept);

    //first one is phase, just applied to the ancilla
    QLGate rz(EBasicOperation::EBO_RZ, -F(2.0) * decomposed[0].m_fCoefficient * t);
    TArray<BYTE> ancillabit;
    ancillabit.AddItem(totalOrder);
    ret.AppendGate(rz, ancillabit);

    for (UINT i = 0; i < trotterStep; ++i)
    {
        for (INT j = 1; j < decomposed.Num(); ++j)
        {
            ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep), ret.m_lstQubits);
        }
    }
    return ret;
}

QLGate QLAPI PauliSimulateGateLeapfrog(const QLMatrix& h, Real t, UINT trotterStep, Real fMinimalKept)
{
    QLGate ret;
    BYTE totalOrder = static_cast<BYTE>(Log2(h.X()));
    ret.AddQubits(totalOrder + 1);
    ret.m_sName = _T("Exp(iHt)");
    Real fOneTrotterStep = t / trotterStep;
    TArray<PauliProduct> decomposed = DecomposePauliProducts(h, fMinimalKept);

    for (INT i = 1; i < decomposed.Num(); ++i)
    {
        for (INT j = i + 1; j < decomposed.Num(); ++j)
        {
            if (abs(decomposed[i].m_fCoefficient) < abs(decomposed[j].m_fCoefficient))
            {
                PauliProduct temp = decomposed[i];
                decomposed[i] = decomposed[j];
                decomposed[j] = temp;
            }
        }
    }

    //first one is phase, just applied to the ancilla
    QLGate rz(EBasicOperation::EBO_RZ, -F(2.0) * decomposed[0].m_fCoefficient * t);
    TArray<BYTE> ancillabit;
    ancillabit.AddItem(totalOrder);
    ret.AppendGate(rz, ancillabit);

    for (UINT i = 0; i < trotterStep; ++i)
    {
        for (INT j = 1; j < decomposed.Num(); ++j)
        {
            if (1 == j)
            {
                if (0 == i)
                {
                    ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep * F(0.5)), ret.m_lstQubits);
                }
                else
                {
                    ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep), ret.m_lstQubits);
                }
            }
            else if (decomposed.Num() == (j + 1))
            {
                ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep), ret.m_lstQubits);
            }
            else
            {
                ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep * F(0.5)), ret.m_lstQubits);
            }
        }

        for (INT j = decomposed.Num() - 1; j > 0; --j)
        {
            if (1 == j)
            {
                if (trotterStep == i + 1)
                {
                    ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep * F(0.5)), ret.m_lstQubits);
                }
            }
            else
            {
                ret.AppendGate(decomposed[j].OneStepGate(fOneTrotterStep * F(0.5)), ret.m_lstQubits);
            }
        }
    }
    return ret;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================