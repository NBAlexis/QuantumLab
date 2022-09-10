//=============================================================================
// FILENAME : SimpleTest.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [10/09/2022 nbale]
//=============================================================================

#include "QuantumLab.h"

int main()
{
    QLGate ha = QLGate(EBasicOperation::EBO_H);
    QLSimulatorMatrix sim(NULL);
    sim.Simulate(NULL);
}


//=============================================================================
// END OF FILE
//=============================================================================