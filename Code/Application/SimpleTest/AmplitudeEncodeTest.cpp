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

void TestCnNOTAncilla()
{
    QLComplex _u[4] = {
        _make_cuComplex(F(0.55729), F(-0.591656)),
        _make_cuComplex(F(-0.58221), F(0.0200615)),
        _make_cuComplex(F(-0.365375), F(-0.453731)),
        _make_cuComplex(F(0.0834975), F(-0.808491))
    };
    QLMatrix u = QLMatrix::CopyCreate(2, 2, _u);

    //u.Print(_T("u"));
    
    QLGate cu = CreateControlledZYZGate(u);
    u = u * u;
    u.Dagger();
    QLGate cui2 = CreateControlledZYZGate(u);

    QLGate ch = CreateCnUWithAncilla(5, cu, cui2, EAncillaType::EAT_Zeroed);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 7;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    //ch.DebugPrint(1);

    output.m_OutputMatrix.Print("res");
}

void TestGroverGate()
{
    QLGate grover = GroverSXGate(3, 5);

    QLGate ch;
    ch.AddQubits(5);
    TArray<BYTE> toadd;
    toadd.Append(ByteSequnce + 2, 3);
    toadd.AddItem(1);
    toadd.AddItem(0);

    QLGate x(EBasicOperation::EBO_X);
    QLGate h(EBasicOperation::EBO_H);

    ch.AppendGate(x, 1);
    ch.AppendGate(h, 1);
    ch.AppendGate(grover, toadd);
    //ch.AppendGate(h, toadd2);
    //ch.AppendGate(x, toadd2);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 5;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    //ch.DebugPrint(1);

    output.m_OutputMatrix.Print("res");
}


//=============================================================================
// END OF FILE
//=============================================================================