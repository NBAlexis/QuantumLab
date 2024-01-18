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

#pragma region CHamitonianList

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
* 
*/
Real CHamitonianList::Measure(const Real* hostWaveFunctionReal, const Real* hostWaveFunctionImagin, INT iRepeat) const
{
    Real ret = F(0.0);
    for (INT i = 0; i < m_lstTerms.Num(); ++i)
    {
        ret += m_lstTerms[i]->Measure(m_pLattice, hostWaveFunctionReal, hostWaveFunctionImagin, iRepeat);
    }
    return ret;
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

#pragma endregion

#pragma region CHamitonianTIM1D

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

#pragma endregion

#pragma region CHamitonianFermion1D

CHamitonianFermion1D::CHamitonianFermion1D(UINT uiSiteCount, const SJordanWeigner1DTerms& param)
    : CHamitonianList()
    , m_sParam(param)
{
    Real fOneOverLatticeSpacing = F(1.0) / m_sParam.m_fLatticeSpacing;
    m_pLattice = new CLattice1D(uiSiteCount);

    CHamitonianStaggeredJordanWigner1DKinetic* pKineticX = new CHamitonianStaggeredJordanWigner1DKinetic(fOneOverLatticeSpacing, EPauliType::EPT_X);
    CHamitonianStaggeredJordanWigner1DKinetic* pKineticY = new CHamitonianStaggeredJordanWigner1DKinetic(fOneOverLatticeSpacing, EPauliType::EPT_Y);

    m_lstTerms.AddItem(pKineticX);
    m_lstTerms.AddItem(pKineticY);

    //TArray<CHamitonianTerm*> terms;
    CHamitonianStaggeredJordanWigner1DPsibarG1Psi* pG1a = new CHamitonianStaggeredJordanWigner1DPsibarG1Psi(param.m_fG1, EPauliType::EPT_Y);
    CHamitonianStaggeredJordanWigner1DPsibarG1Psi* pG1b = new CHamitonianStaggeredJordanWigner1DPsibarG1Psi(param.m_fG1, EPauliType::EPT_X);
    m_lstTerms.AddItem(pG1a);
    m_lstTerms.AddItem(pG1b);
    if (abs(param.m_fG1) > F(1.0e-6))
    {
        m_sParam.m_bHasG1 = TRUE;
    }

    //we put mass g0 and f1 at last because [g0, f1] = [mass, g0] = [mass, f1] = 0, they are commutable
    CHamitonianStaggeredJordanWigner1DPsibarPsi* pMassTerm = new CHamitonianStaggeredJordanWigner1DPsibarPsi(param.m_fMass);
    m_lstTerms.AddItem(pMassTerm);
    if (abs(param.m_fMass) > F(1.0e-6))
    {
        m_sParam.m_bHasMass = TRUE;
    }

    CHamitonianStaggeredJordanWigner1DPsibarG0Psi* pG0 = new CHamitonianStaggeredJordanWigner1DPsibarG0Psi(param.m_fG0);
    m_lstTerms.AddItem(pG0);
    if (abs(param.m_fG0) > F(1.0e-6))
    {
        m_sParam.m_bHasG0 = TRUE;
    }

    CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure* pContact = new CHamitonianStaggeredJordanWigner1DPsibarPsiSqaure(param.m_f1Sq * fOneOverLatticeSpacing);
    m_lstTerms.AddItem(pContact);
    if (abs(param.m_f1Sq) > F(1.0e-6))
    {
        m_sParam.m_bHasF1sq = TRUE;
    }
}

void CHamitonianFermion1D::QLGateA(QLGate& gate, Real fTime, UINT uiControllerCount) const
{
    TArray<BYTE> addEachTerm;
    addEachTerm.Append(ByteSequnce, uiControllerCount + 1);

    gate.AppendGate(m_lstTerms[0]->BuildCircuit(m_pLattice, fTime), addEachTerm);
}

void CHamitonianFermion1D::QLGateB(QLGate& gate, Real fTime, UINT uiControllerCount) const
{
    TArray<BYTE> addEachTerm;
    addEachTerm.Append(ByteSequnce, uiControllerCount + 1);

    gate.AppendGate(m_lstTerms[1]->BuildCircuit(m_pLattice, fTime), addEachTerm);

    if (m_sParam.m_bHasG0 || m_sParam.m_bHasF1sq || m_sParam.m_bHasMass)
    {
        if (m_sParam.m_bHasG1)
        {
            gate.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fTime), addEachTerm);
            gate.AppendGate(m_lstTerms[3]->BuildCircuit(m_pLattice, fTime), addEachTerm);
        }
    }
    else
    {
        gate.AppendGate(m_lstTerms[2]->BuildCircuit(m_pLattice, fTime), addEachTerm);
    }

}

void CHamitonianFermion1D::QLGateC(QLGate& gate, Real fTime, UINT uiControllerCount) const
{
    TArray<BYTE> addEachTerm;
    addEachTerm.Append(ByteSequnce, uiControllerCount + 1);

    if (m_sParam.m_bHasG0 || m_sParam.m_bHasF1sq || m_sParam.m_bHasMass)
    {
        if (m_sParam.m_bHasG0)
        {
            gate.AppendGate(m_lstTerms[4]->BuildCircuit(m_pLattice, fTime), addEachTerm);
        }
        if (m_sParam.m_bHasF1sq)
        {
            gate.AppendGate(m_lstTerms[5]->BuildCircuit(m_pLattice, fTime), addEachTerm);
        }
        if (m_sParam.m_bHasMass)
        {
            gate.AppendGate(m_lstTerms[6]->BuildCircuit(m_pLattice, fTime), addEachTerm);
        }
    }
    else
    {
        gate.AppendGate(m_lstTerms[3]->BuildCircuit(m_pLattice, fTime), addEachTerm);
    }
}

/**
* simplify by notice that mass term and chemical term are commutable
*/
QLGate CHamitonianFermion1D::BuildSimulationCircuitABCBA(Real fTime, UINT uiTrotterStep) const
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
                QLGateA(ret, fOneTrotterStep * F(0.5), uiControllerCount);
                QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
                QLGateC(ret, fOneTrotterStep, uiControllerCount);
                QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
                QLGateA(ret, fOneTrotterStep * F(0.5), uiControllerCount);
            }
            else
            {
                // a b C b A
                QLGateA(ret, fOneTrotterStep * F(0.5), uiControllerCount);
                QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
                QLGateC(ret, fOneTrotterStep, uiControllerCount);
                QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
                QLGateA(ret, fOneTrotterStep, uiControllerCount);
            }
        }
        else if (i == uiTrotterStep - 1)
        {
            //b C b a
            QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
            QLGateC(ret, fOneTrotterStep, uiControllerCount);
            QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
            QLGateA(ret, fOneTrotterStep * F(0.5), uiControllerCount);
        }
        else
        {
            //b C b A
            QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
            QLGateC(ret, fOneTrotterStep, uiControllerCount);
            QLGateB(ret, fOneTrotterStep * F(0.5), uiControllerCount);
            QLGateA(ret, fOneTrotterStep, uiControllerCount);
        }
    }
    return ret;
}

#pragma endregion

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================