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

void TestAmplitudeEncodeOneRealVector()
{
    QLMatrix m(1, 8);
    m.RandomOneReal();
    m = m / _sqrt(m.VectorDot(m).x);
    m.Print(_T("v"));
    TArray<QLComplex> v = m.ToVector();

    QLGate ae = AmplitudeEncodeOneVectorReal(v.GetData(), 3);

    QLMatrix res = QLSimulatorVector::ShowState(ae);
    res.Print(_T("resv"));
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

void TestSimpleEncode()
{
    //-2.219428427064474574e-01,-2.080759378357706357e-02,
    // 9.842508935813076842e-02,1.096014891232260924e-01,
    // -5.563325360209519371e-02,2.340083368805381814e-02,
    // -8.668671161025723326e-02,-2.495999396951170124e-01,
    // 2.267926831618269079e-01,1.711864065389313849e-02,
    // -8.363397026485545893e-02,-7.620323975045822928e-02

    QLComplex onev[6] = { 
        _make_cuComplex(F(-2.219428427064474574e-01), F(-2.080759378357706357e-02)),
        _make_cuComplex(F(9.842508935813076842e-02), F(1.096014891232260924e-01)),
        _make_cuComplex(F(-5.563325360209519371e-02), F(2.340083368805381814e-02)),
        _make_cuComplex(F(-8.668671161025723326e-02), F(-2.495999396951170124e-01)),
        _make_cuComplex(F(2.267926831618269079e-01), F(1.711864065389313849e-02)),
        _make_cuComplex(F(-8.363397026485545893e-02), F(-7.620323975045822928e-02))
    };
    QLGate simpleencodeOne = SimpleEncodeOneVector(onev, 6, 6);
    QLMatrix res = QLSimulatorVector::ShowState(simpleencodeOne);
    res.Print(_T("res"));
}

//=============================================================================
// END OF FILE
//=============================================================================