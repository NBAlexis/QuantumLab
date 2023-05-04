//=============================================================================
// FILENAME : AmplitudeEncodeTest.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [12/03/2023 nbale]
//=============================================================================

#include "SimpleTest.h"

void TestAmplitudeEncodeOneVector()
{
    QLMatrix m(1, 8);
    m.RandomOne();
    m = m / _sqrt(m.VectorDot(m).x);
    m.Print(_T("v"));
    TArray<QLComplex> v = m.ToVector();

    QLGate ae = AmplitudeEncodeOneVector(v.GetData(), 3, TRUE);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ae;

    QLSimulatorVector sim;
    sim.Simulate(&param);
}

void TestAmplitudeEncodeVectors()
{
    QLMatrix m1(8, 8);
    m1.RandomUnitary();
    QLMatrix m2(8, 8);
    m2.RandomUnitary();
    QLMatrix m = m1.BlockAdd(m2, 1);

    m.Print(_T("m"));

    QLSimulatorOutputVector out;
    QLGate ae = AmplitudeEncodeVectors(m.ToVector().GetData(), 4, 3, TRUE);
    QLSimulatorParametersVector param;
    param.m_byQubitCount = 7;
    param.m_MasterGate = ae;
    param.m_bPrint = FALSE;

    QLSimulatorVector sim;
    sim.Simulate(&param, &out);
    out.m_OutputMatrix = out.m_OutputMatrix * 4;
    out.m_OutputMatrix.ReShape(16, 8);
    out.m_OutputMatrix = out.m_OutputMatrix - m;
    out.m_OutputMatrix.Transpose();
    out.m_OutputMatrix.Print(_T("delta"));

    ae = AmplitudeEncodeVectors(m.ToVector().GetData(), 4, 3, FALSE);
    param.m_MasterGate = ae;
    sim.Simulate(&param, &out);
    out.m_OutputMatrix = out.m_OutputMatrix * 4;
    out.m_OutputMatrix.ReShape(16, 8);
    out.m_OutputMatrix.ElementDiv(m);
    out.m_OutputMatrix.Transpose();
    out.m_OutputMatrix.Print(_T("res"));

}


//=============================================================================
// END OF FILE
//=============================================================================