//=============================================================================
// FILENAME : QLSimulator.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [11/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

/**
* We use the leapfrog decompose instead of Trotter,
* 
* H = H1 + H2 + H3 + ... + Hn
* 
* exp(iHt) = exp(iH1t/2m)exp(iH2t/2m)...exp(iHn-1t/2m) exp(iHn/m)  for (m-1) steps exp(iH1t/m)exp(iH2t/m)...exp(iHnt/m) exp(iHn-1t/2m)...exp(iH2t/2m)exp(iH1t/2m)
*/
QLGate CHamitonianList::BuildSimulationCircuit(Real fTime, UINT uiTrotterStep) const
{
    switch (m_eDecompose)
    {
    case ETimeDecomposeType::ETDT_ABCBA:
        return BuildSimulationCircuitABCBA(fTime, uiTrotterStep);
    case ETimeDecomposeType::ETDT_Trotter:
    default:
        return BuildSimulationCircuitTrotter(fTime, uiTrotterStep);
    }
}

CHamitonianList::~CHamitonianList()
{
    for (INT i = 0; i < m_lstTerms.Num(); ++i)
    {
        appSafeDelete(m_lstTerms[i]);
    }
    appSafeDelete(m_pLattice);
}

/**
* H = H1 + H2 + H3 + ... + Hn
* exp(iHt) = [exp(iH1t/m)exp(iH2t/m)...exp(iHnt/m)]^m
*/
QLGate CHamitonianList::BuildSimulationCircuitTrotter(Real fTime, UINT uiTrotterStep) const
{
    QLGate ret;
    UINT uiControllerCount = m_pLattice->GetControllerCount();
    ret.AddQubits(static_cast<BYTE>(uiControllerCount + 1));
    ret.m_sName = _T("Exp(iHt)");
    Real fOneTrotterStep = fTime / uiTrotterStep;

    TArray<BYTE> addEachTerm;
    addEachTerm.Append(ByteSequnce, uiControllerCount + 1);

    for (UINT i = 0; i < uiTrotterStep; ++i)
    {
        for (INT j = 0; j < m_lstTerms.Num(); ++j)
        {
            ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
        }
    }
    return ret;
}

/**
* ABCBA ABCBA ABCBA
*/
QLGate CHamitonianList::BuildSimulationCircuitABCBA(Real fTime, UINT uiTrotterStep) const
{
    QLGate ret;
    UINT uiControllerCount = m_pLattice->GetControllerCount();
    ret.AddQubits(static_cast<BYTE>(uiControllerCount + 1));
    ret.m_sName = _T("Exp(iHt)");
    Real fOneTrotterStep = fTime / uiTrotterStep;

    TArray<BYTE> addEachTerm;
    addEachTerm.Append(ByteSequnce, uiControllerCount + 1);

    for (UINT i = 0; i < uiTrotterStep; ++i)
    {
        if (0 == i)
        {
            if (1 == uiTrotterStep)
            {
                // a b C b a
                for (INT j = 0; j < m_lstTerms.Num() - 1; ++j)
                {
                    ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                }

                ret.AppendGate(m_lstTerms[m_lstTerms.Num() - 1]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

                for (INT j = m_lstTerms.Num() - 2; j >= 0; --j)
                {
                    ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                }
            }
            else
            {
                // a b C b A
                for (INT j = 0; j < m_lstTerms.Num() - 1; ++j)
                {
                    ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                }

                ret.AppendGate(m_lstTerms[m_lstTerms.Num() - 1]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

                for (INT j = m_lstTerms.Num() - 2; j > 0; --j)
                {
                    ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                }

                ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
            }
        }
        else if (i == uiTrotterStep - 1)
        {
            //b C b a
            for (INT j = 1; j < m_lstTerms.Num() - 1; ++j)
            {
                ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            }

            ret.AppendGate(m_lstTerms[m_lstTerms.Num() - 1]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

            for (INT j = m_lstTerms.Num() - 2; j >= 0; --j)
            {
                ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            }
        }
        else
        {
            //b C b A
            for (INT j = 1; j < m_lstTerms.Num() - 1; ++j)
            {
                ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            }

            ret.AppendGate(m_lstTerms[m_lstTerms.Num() - 1]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

            for (INT j = m_lstTerms.Num() - 2; j > 0; --j)
            {
                ret.AppendGate(m_lstTerms[j]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            }

            ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
        }
    }
    return ret;
}

CHamitonianTIM1D::CHamitonianTIM1D(UINT uiSiteCount, Real fLambda)
    : CHamitonianList()
{
    m_pLattice = new CLattice1D(uiSiteCount);
    CHamitonianPauliNeighbour* pPauliNeighbour = new CHamitonianPauliNeighbour();
    m_lstTerms.AddItem(pPauliNeighbour);
    CHamitonianPauli* pTransverseField = new CHamitonianPauli(fLambda);
    m_lstTerms.AddItem(pTransverseField);
}

void CHamitonianTIM1D::SetLambda(Real fLambda)
{
    m_lstTerms[1]->m_fCoefficient = fLambda;
}

CHamitonianNJL1D::CHamitonianNJL1D(UINT uiSiteCount, Real fMass, Real fMu, Real fg, Real fLatticeSpacing)
    : CHamitonianList()
    , m_fG(fg)
{
    Real fOneOverLatticeSpacing = F(1.0) / fLatticeSpacing;
    m_pLattice = new CLattice1D(uiSiteCount);

    CHamitonianStaggeredJordanWigner1DKinetic* pKineticX = new CHamitonianStaggeredJordanWigner1DKinetic(fOneOverLatticeSpacing, EPauliType::EPT_X);
    CHamitonianStaggeredJordanWigner1DKinetic* pKineticY = new CHamitonianStaggeredJordanWigner1DKinetic(fOneOverLatticeSpacing, EPauliType::EPT_Y);

    CHamitonianStaggeredJordanWigner1DPsibarPsi* pMassTerm = new CHamitonianStaggeredJordanWigner1DPsibarPsi(fMass);

    CHamitonianStaggeredJordanWigner1DPsibarG0Psi* pChemical = new CHamitonianStaggeredJordanWigner1DPsibarG0Psi(fMu);

    CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure* pContact = new CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure(fg * fOneOverLatticeSpacing);

    m_lstTerms.AddItem(pContact);
    m_lstTerms.AddItem(pKineticX);
    m_lstTerms.AddItem(pKineticY);
    m_lstTerms.AddItem(pMassTerm);
    m_lstTerms.AddItem(pChemical);

}

void CHamitonianNJL1D::SetLatticeSpacing(Real fa)
{
    Real fOneOverLatticeSpacing = F(1.0) / fa;
    m_lstTerms[1]->m_fCoefficient = fOneOverLatticeSpacing;
    m_lstTerms[2]->m_fCoefficient = fOneOverLatticeSpacing;
    m_lstTerms[0]->m_fCoefficient = m_fG * fOneOverLatticeSpacing;
}

void CHamitonianNJL1D::SetMass(Real fm)
{
    m_lstTerms[3]->m_fCoefficient = fm;
}

void CHamitonianNJL1D::SetChemicalPotential(Real fMu)
{
    m_lstTerms[4]->m_fCoefficient = fMu;
}

void CHamitonianNJL1D::SetContactCoupling(Real fg)
{
    m_fG = fg;
    m_lstTerms[0]->m_fCoefficient = fg * m_lstTerms[1]->m_fCoefficient;
}

/**
* simplify by notice that mass term and chemical term are commutable
*/
QLGate CHamitonianNJL1D::BuildSimulationCircuitABCBA(Real fTime, UINT uiTrotterStep) const
{
    QLGate ret;
    UINT uiControllerCount = m_pLattice->GetControllerCount();
    ret.AddQubits(static_cast<BYTE>(uiControllerCount + 1));
    ret.m_sName = _T("Exp(iHt)");
    Real fOneTrotterStep = fTime / uiTrotterStep;

    TArray<BYTE> addEachTerm;
    addEachTerm.Append(ByteSequnce, uiControllerCount + 1);

    for (UINT i = 0; i < uiTrotterStep; ++i)
    {
        if (0 == i)
        {
            if (1 == uiTrotterStep)
            {
                // a b C b a
                ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

                ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

                ret.AppendGate(m_lstTerms[3]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
                ret.AppendGate(m_lstTerms[4]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

                ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

                ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            }
            else
            {
                // a b C b A
                ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

                ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

                ret.AppendGate(m_lstTerms[3]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
                ret.AppendGate(m_lstTerms[4]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

                ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
                ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

                ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
            }
        }
        else if (i == uiTrotterStep - 1)
        {
            //b C b a
            ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

            ret.AppendGate(m_lstTerms[3]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
            ret.AppendGate(m_lstTerms[4]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

            ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

            ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
        }
        else
        {
            //b C b A
            ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

            ret.AppendGate(m_lstTerms[3]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
            ret.AppendGate(m_lstTerms[4]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);

            ret.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);
            ret.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fOneTrotterStep * F(0.5)), addEachTerm);

            ret.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fOneTrotterStep), addEachTerm);
        }
    }
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================