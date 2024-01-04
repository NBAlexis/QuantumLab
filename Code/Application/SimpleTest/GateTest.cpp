//=============================================================================
// FILENAME : GateTest.h
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [04/01/2024 nbale]
//=============================================================================

#include "SimpleTest.h"

void ShowMatrixOfGate()
{
    QLGate ch = CreateControlledSwap(1);
    ch.DebugPrint(2);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;

    QLSimulatorOutputMatrix output;

    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    output.m_OutputMatrix.Print("");
}

void CSDOfGate()
{
    QLComplex m[64] = {
        _onec,  _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _zeroc,
        _zeroc, _onec,  _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _zeroc,
        _zeroc, _zeroc, _onec,  _zeroc, _zeroc, _zeroc, _zeroc, _zeroc,
        _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _onec,  _zeroc, _zeroc,
        _zeroc, _zeroc, _zeroc, _zeroc, _onec,  _zeroc, _zeroc, _zeroc,
        _zeroc, _zeroc, _zeroc, _onec,  _zeroc, _zeroc, _zeroc, _zeroc,
        _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _onec,  _zeroc,
        _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _zeroc, _onec
    };

    QLMatrix mtr = QLMatrix::CopyCreate(8, 8, m);
    mtr.Print();

    QLGate ch = CSDDecompose(mtr, 3);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;

    QLSimulatorOutputMatrix output;

    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    output.m_OutputMatrix.Print();

    ch.DebugPrint(-1);
}


//=============================================================================
// END OF FILE
//=============================================================================