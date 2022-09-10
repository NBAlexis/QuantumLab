//=============================================================================
// FILENAME : QLSimulatorMatrix.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [10/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

void QLSimulatorMatrix::Simulate(const QLSimulatorParameters * params) const
{
    const QLSimulatorParametersMatrix* param = dynamic_cast<const QLSimulatorParametersMatrix*>(params);

    QuESTEnv evn = createQuESTEnv();

    Qureg mat = createDensityQureg(1, evn);

    hadamard(mat, 0);
    hadamard(mat, 0);
    //collapseToOutcome(mat, 0, 1);

    copyStateFromGPU(mat);
    syncQuESTEnv(evn);

    for (LONGLONG i = 0; i < mat.numAmpsPerChunk; ++i)
    {
        *m_pStdOut << mat.stateVec.real[i];
        *m_pStdOut << mat.stateVec.imag[i];
    }

    destroyQureg(mat, evn);
    destroyQuESTEnv(evn);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================