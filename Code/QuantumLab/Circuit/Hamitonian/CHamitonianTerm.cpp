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

Real CHamitonianTerm::OneTermMeasure(const TArray<BYTE>& pauliType, const Real* hostWaveFunctionReal, const Real* hostWaveFunctionImagin, INT iRepeat)
{
    PauliProduct prod(pauliType, F(1.0));
    QLGate projection = prod.Project();
    TArray<SBasicOperation> ops = projection.GetOperation(projection.m_lstQubits);
    SIZE_T opssize = ops.Num();

    QuESTEnv evn = createQuESTEnv();
    UINT uiQubit = static_cast<UINT>(projection.m_lstQubits.Num());
    Qureg vec = createQureg(uiQubit, evn);

    LONGLONG veclen = 1LL << static_cast<UINT>(uiQubit);
    memcpy(vec.stateVec.real, hostWaveFunctionReal, sizeof(Real) * veclen);
    memcpy(vec.stateVec.imag, hostWaveFunctionImagin, sizeof(Real) * veclen);

    copyStateToGPU(vec);
    syncQuESTEnv(evn);

    for (SIZE_T i = 0; i < opssize; ++i)
    {
        QLGate::PerformBasicOperation(vec, ops[static_cast<INT>(i)]);
    }

    UINT uiOdd = 0;
    UINT uiEven = 0;

    if (iRepeat <= 0)
    {
        TArray<INT> tomeasure;
        for (INT qubit = 0; qubit < static_cast<INT>(projection.m_lstQubits.Num()); ++qubit)
        {
            tomeasure.AddItem(qubit);
        }
        Real* resprob = (Real*)(malloc(sizeof(Real) * veclen));
        calcProbOfAllOutcomes(resprob, vec, tomeasure.GetData(), tomeasure.Num());
        Real fRes = F(0.0);
        for (UINT measureRes = 0; measureRes < veclen; ++measureRes)
        {
            if (measureRes & 1)
            {
                fRes -= resprob[measureRes];
            }
            else
            {
                fRes += resprob[measureRes];
            }
        }
        destroyQureg(vec, evn);
        destroyQuESTEnv(evn);
        return fRes;
    }

    syncQuESTEnv(evn);
    copyStateFromGPU(vec);

    for (INT i = 0; i < iRepeat; ++i)
    {
        if (0 != i)
        {
            copyStateToGPU(vec);
            syncQuESTEnv(evn);
        }

        UINT measureRes = 0;
        for (INT j = 0; j < static_cast<INT>(uiQubit); ++j)
        {
            INT out = measure(vec, j);
            if (1 == out)
            {
                ++measureRes;
            }
        }

        if (measureRes & 1)
        {
            ++uiOdd;
        }
        else
        {
            ++uiEven;
        }
    }

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);

    return (uiEven - uiOdd) / static_cast<Real>(uiOdd + uiEven);
}

QLGate CHamitonianTerm::BuildCircuit(const CLattice* pLattice, Real fTrotterTime) const
{
    const TArray<PauliProduct> allTerms = GetAllTerms(pLattice);
    UINT uiControllerCount = pLattice->GetControllerCount();

    TArray<BYTE> toAdd;
    toAdd.Append(ByteSequnce, uiControllerCount + 1);
    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(uiControllerCount + 1));

    for (INT i = 0; i < allTerms.Num(); ++i)
    {
        ret.AppendGate(allTerms[i].OneStepGate(fTrotterTime), toAdd);
    }

    return ret;
}

Real CHamitonianTerm::Measure(const CLattice* pLattice, const Real* hostWaveFunctionReal, const Real* hostWaveFunctionImagin, INT iRepeat) const
{
    const TArray<PauliProduct> allTerms = GetAllTerms(pLattice);
    //UINT uiControllerCount = pLattice->GetControllerCount();

    Real ret = F(0.0);
    for (INT i = 0; i < allTerms.Num(); ++i)
    {
        ret += allTerms[i].m_fCoefficient * OneTermMeasure(allTerms[i].m_lstPauliType, hostWaveFunctionReal, hostWaveFunctionImagin, iRepeat);
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================