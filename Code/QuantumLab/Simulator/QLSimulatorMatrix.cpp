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
    std::vector<BYTE> qubits;
    for (BYTE byQ = 0; byQ < param->m_byQubitCount; ++byQ)
    {
        qubits.push_back(byQ);
    }
    std::vector<SBasicOperation> ops = param->m_MasterGate.GetOperation(qubits);
    SIZE_T opssize = ops.size();

    QuESTEnv evn = createQuESTEnv();
    Qureg vec = createQureg(param->m_byQubitCount, evn);
    

    //This is a lazy slow implement, I need to use cuda to improve it
    *m_pStdOut << "{\n";
    LONGLONG veclen = 1LL << param->m_byQubitCount;
    for (LONGLONG line = 0; line < veclen; ++line)
    {
        syncQuESTEnv(evn);
        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            if (line2 == line)
            {
                vec.stateVec.real[line2] = 1.0;
                vec.stateVec.imag[line2] = 0.0;
            }
            else
            {
                vec.stateVec.real[line2] = 0.0;
                vec.stateVec.imag[line2] = 0.0;
            }
        }
        copyStateToGPU(vec);

        for (SIZE_T i = 0; i < opssize; ++i)
        {
            QLGate::PerformBasicOperation(vec, ops[i]);
        }
        syncQuESTEnv(evn);
        copyStateFromGPU(vec);

        for (LONGLONG line2 = 0; line2 < veclen; ++line2)
        {
            if (0 == line2)
            {
                *m_pStdOut << "{";
            }
            else 
            {
                *m_pStdOut << ", ";
            }
            *m_pStdOut << PrintComplex(vec.stateVec.real[line2], vec.stateVec.imag[line2]);
        }
        if (line == (veclen - 1))
        {
            *m_pStdOut << "}\n}\n";
        }
        else
        {
            *m_pStdOut << "},\n";
        }
    }

    destroyQureg(vec, evn);
    destroyQuESTEnv(evn);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================