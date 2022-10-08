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

void TestCSD2By1()
{
    QLMatrix m(4, 10);

    m.RandomUnitary();
    m.Print("m");
    QLMatrix q1, q2, c, s, v;
    m.CSD2BY1(q1, q2, c, s, v, 6);
    QLMatrix q = q1.BlockAdd(q2, 0);
    q.Print("q");
    QLMatrix cs = c.BlockAdd(s, 2);
    cs.Print("cs");
    v.Print("v");

}


void TestCSD()
{
    QLMatrix m(8, 8);
    UINT xSep = 4;
    UINT ySep = 4;

    m.RandomUnitary();
    m.Print("m");
    QLMatrix q1, q2, c, s, v1, v2;
    m.CSD(q1, q2, c, s, v1, v2, xSep, ySep);

    QLMatrix q, cs, v;

    QLMatrix::CombinedCSD(q, cs, v, q1, q2, c, s, v1, v2, xSep, ySep);

    q.Print("q");
    cs.Print("cs");
    v.Print("v");
}

void TestRYGate()
{
    QLGate ch = QLGate(EBasicOperation::EBO_Phase, 0.3);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 1;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestZYZGate()
{
    QLMatrix m(2, 2);
    m.RandomUnitary();

    m.Print("m");

    QLGate ch = CreateZYZGate(m);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 1;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestControoledZYZ()
{
    QLMatrix m(2, 2);
    m.RandomUnitary();
    m.Print("m");
    QLGate ch = CreateControlledZYZGate(m);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 2;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestControoledHadamard()
{
    QLGate ha = QLGate(EBasicOperation::EBO_H);
    TArray<BYTE> testtarget0;
    testtarget0.AddItem(0);
    TArray<BYTE> testtarget1;
    testtarget1.AddItem(1);
    QLGate ch = ha.Controlled(1, testtarget0);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 2;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestSwap()
{
    QLGate ch = CreateSwapGate();
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 2;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestToffoli()
{
    QLGate ch = CreateToffoliGate();
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestMatrixSqrt()
{
    QLMatrix m(2, 2);
    m.RandomUnitary();
    m.Print("m");

    QLMatrix v = Sqrt2by2(m);
    v = v * v;
    v.Print("v");
}

UINT TestCnU()
{
    QLMatrix m(2, 2);
    m.RandomUnitary();

    /**
    * ccm expect:
    * 
    * 1 0 0 0 0 0 0
    * 0 1 0 0 0 0 0
    * 0 0 1 0 0 0 0
    * 0 0 0 a 0 0 b
    * 0 0 0 0 1 0 0
    * 0 0 0 0 0 1 0
    * 0 0 0 c 0 0 d
    */
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, m.Get(0, 0));
    correct.Set(15, 7, m.Get(1, 0));
    correct.Set(7, 15, m.Get(0, 1));
    correct.Set(15, 15, m.Get(1, 1));


    QLGate ch = CreateCnU(3, m);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnRX()
{
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));
    correct.Set(15, 7, _make_cuComplex(F(0.0), -_sin(F(0.1) / 2)));
    correct.Set(7, 15, _make_cuComplex(F(0.0), -_sin(F(0.1) / 2)));
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));

    QLGate ch = CreateCnRX(3, F(0.1));
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnRY()
{
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));
    correct.Set(15, 7, _make_cuComplex(-_sin(F(0.1) / 2), F(0.0)));
    correct.Set(7, 15, _make_cuComplex(_sin(F(0.1) / 2), F(0.0)));
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));

    QLGate ch = CreateCnRY(3, F(0.1));
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnRZ()
{
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, _make_cuComplex(_cos(F(0.1) / 2), -_sin(F(0.1) / 2)));
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1) / 2), _sin(F(0.1) / 2)));

    QLGate ch = CreateCnRZ(3, F(0.1));
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnP()
{
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1)), _sin(F(0.1))));

    QLGate ch = CreateCnP(3, F(0.1));
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnPh()
{
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, _make_cuComplex(_cos(F(0.1)), _sin(F(0.1))));
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1)), _sin(F(0.1))));

    QLGate ch = CreateCnPh(3, F(0.1));
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnU2()
{
    QLMatrix m(2, 2);
    m.RandomUnitary();

    /**
    * ccm expect:
    *
    * 1 0 0 0 0 0 0
    * 0 1 0 0 0 0 0
    * 0 0 1 0 0 0 0
    * 0 0 0 a 0 0 b
    * 0 0 0 0 1 0 0
    * 0 0 0 0 0 1 0
    * 0 0 0 c 0 0 d
    */
    //QLMatrix correct = QLMatrix::CreateEye(16, 16);
    //correct.Set(7, 7, m.Get(0, 0));
    //correct.Set(15, 7, m.Get(1, 0));
    //correct.Set(7, 15, m.Get(0, 1));
    //correct.Set(15, 15, m.Get(1, 1));
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, m.Get(0, 0));
    correct.Set(15, 7, m.Get(1, 0));
    correct.Set(7, 15, m.Get(0, 1));
    correct.Set(15, 15, m.Get(1, 1));

    QLGate ch = CreateZYZGate(m);
    QLGate ch1 = ch.CreateControlled().CreateControlled().CreateControlled();
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch1;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    //correct.Print("c");
    //output.m_OutputMatrix.Print("o");

    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    appGeneral(_T("delta = %f\n"), fDelta);
    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnRY2()
{
    QLMatrix correct = QLMatrix::CreateEye(16, 16);
    correct.Set(7, 7, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));
    correct.Set(15, 7, _make_cuComplex(-_sin(F(0.1) / 2), F(0.0)));
    correct.Set(7, 15, _make_cuComplex(_sin(F(0.1) / 2), F(0.0)));
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));

    QLGate ch0 = QLGate(EBasicOperation::EBO_RY, F(0.1)).CreateControlled().CreateControlled();
    QLGate ch = ch0.CreateControlled();

    //ch0.DebugPrint(1);
    //ch.DebugPrint(1);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnRY3()
{
    QLMatrix correct = QLMatrix::CreateEye(8, 8);
    correct.Set(3, 3, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));
    correct.Set(7, 3, _make_cuComplex(_sin(F(0.1) / 2), F(0.0)));
    correct.Set(3, 7, _make_cuComplex(-_sin(F(0.1) / 2), F(0.0)));
    correct.Set(7, 7, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));

    QLGate ch0 = QLGate(EBasicOperation::EBO_RY, F(0.1)).CreateControlled();
    //QLGate ch1 = ch0.CreateControlled();
    ch0.Dagger();
    QLGate ch = ch0.CreateControlled();

    //ch0.DebugPrint(1);
    //ch.DebugPrint(1);
    //ch1.DebugPrint(1);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

UINT TestCnRY4()
{
    QLMatrix correct = QLMatrix::CreateEye(32, 32);
    correct.Set(15, 15, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));
    correct.Set(31, 15, _make_cuComplex(-_sin(F(0.1) / 2), F(0.0)));
    correct.Set(15, 31, _make_cuComplex(_sin(F(0.1) / 2), F(0.0)));
    correct.Set(31, 31, _make_cuComplex(_cos(F(0.1) / 2), F(0.0)));

    QLGate ch = CreateCnRY(4, F(0.1));

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 5;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    correct = correct - output.m_OutputMatrix;
    Real fDelta = abs(correct.VectorDot(correct).x);

    return fDelta > F(0.0000001) ? 1 : 0;
}

void TestCP()
{
    QLGate ch = QLGate(EBasicOperation::EBO_CP, 0.3);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 2;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
}

void TestSqrtNot()
{
    QLMatrix sqrt = Sqrt2by2(_PauliX);
    sqrt.Print();

    QLMatrix back = sqrt * sqrt;
    back.Print("back");
}

void TestCCPhase()
{
    QLGate ch = CreateCnPh(3, F(0.1));
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
    //ch.DebugPrint(1);
}

void TestFRy()
{
    TArray <Real> angles;
    angles.AddItem(F(0.1));
    angles.AddItem(F(0.2));
    angles.AddItem(F(0.3));
    angles.AddItem(F(0.4));
    angles.AddItem(F(0.5));
    angles.AddItem(F(0.6));
    angles.AddItem(F(0.7));
    angles.AddItem(F(0.8));

    TArray <Real> angles2;
    angles2.AddItem(F(0.11));
    angles2.AddItem(F(0.12));
    angles2.AddItem(F(0.13));
    angles2.AddItem(F(0.14));
    angles2.AddItem(F(0.15));
    angles2.AddItem(F(0.16));
    angles2.AddItem(F(0.17));
    angles2.AddItem(F(0.18));

    QLGate fryz = FRyz(angles, angles2, 4);
    //QLGate fryz = FRy(angles, 4);

    /**
    TArray <BYTE> reverse;
    reverse.AddItem(3);
    reverse.AddItem(2);
    reverse.AddItem(1);
    reverse.AddItem(0);
    QLGate ch;
    ch.AddQubits(4);
    ch.AppendGate(fry, reverse);
    */

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    //param.m_MasterGate = ch;
    param.m_MasterGate = fryz;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
    //fry.DebugPrint(1);
}

void TestFRy2()
{
    TArray <Real> angles;
    angles.AddItem(F(0.1));
    angles.AddItem(F(0.2));
    angles.AddItem(F(0.3));
    angles.AddItem(F(0.4));
    angles.AddItem(F(0.5));
    angles.AddItem(F(0.6));
    angles.AddItem(F(0.7));
    angles.AddItem(F(0.8));

    QLGate fry = FRy(angles, 4);
    TArray<BYTE> applied;
    applied.AddItem(2);
    applied.AddItem(1);
    applied.AddItem(0);
    applied.AddItem(3);
    fry.ApplyOnQubits(applied);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = fry;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
    //fry.DebugPrint(1);

    QLGate x(EBasicOperation::EBO_X);
    QLGate fry2;
    fry2.AddQubits(4);
    for (UINT i = 0; i < 8; ++i)
    {
        for (BYTE j = 0; j < 3; ++j)
        {
            if (0 == (i & (1U << j)))
            {
                TArray<BYTE> xbit;
                xbit.AddItem(j);
                fry2.AppendGate(x, xbit);
            }
        }

        QLGate cnry = CreateCnRY(3, angles[i]);
        fry2.AppendGate(cnry, fry2.m_lstQubits);

        for (BYTE j = 0; j < 3; ++j)
        {
            if (0 == (i & (1U << j)))
            {
                TArray<BYTE> xbit;
                xbit.AddItem(j);
                fry2.AppendGate(x, xbit);
            }
        }
    }

    param.m_MasterGate = fry2;
    sim.Simulate(&param);
}

void TestCSDGate()
{
    QLMatrix m(32, 32);
    m.RandomUnitary();
    
    QLGate ch = CSDDecompose(m, 5);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 5;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;

    QLSimulatorOutputMatrix output;

    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
    
    m = m - output.m_OutputMatrix;
    QLComplex norm = m.VectorDot(m);
    appGeneral(_T("||d||=%f"), norm.x);
}

void TestAmplitudeEncode()
{
    QLMatrix m(1, 8);
    m.RandomOne();
    m = m / _sqrt(m.VectorDot(m).x);
    m.Print(_T("v"));
    TArray<QLComplex> v = m.ToVector();

    QLGate ae = AmplitudeEncode(v);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ae;

    QLSimulatorVector sim;
    sim.Simulate(&param);
}

void TestQFFT()
{
    QLMatrix m(1, 8);
    m.RandomOne();
    m = m / _sqrt(m.VectorDot(m).x);
    m.Print(_T("m"));
    QLMatrix fft = m.VectorFFT(FALSE);
    fft = fft / _sqrt(fft.VectorDot(fft).x);
    fft.Print(_T("fft"));

    QLGate ae = AmplitudeEncode(m.ToVector());

    QLGate fftg = QuantumFFTGate(3);

    QLGate all;
    all.AddQubits(3);
    TArray<BYTE> allqubits;
    allqubits.AddItem(0);
    allqubits.AddItem(1);
    allqubits.AddItem(2);
    all.AppendGate(ae, allqubits);
    all.AppendGate(fftg, allqubits);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = all;

    QLSimulatorVector sim;
    sim.Simulate(&param);
}

void TestQFFT2()
{
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = QuantumFFTGate(3);
    param.m_bPrint = TRUE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);
}

void TestEVD()
{
    QLMatrix m(8, 8);
    m.RandomOne();
    m.Print("m0");
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");


    m.MatrixIExp(F(0.5));

    m.Print("mexp");
}

void TestKron()
{
    QLComplex m1data[4] = { _mcr(1.0), _mcr(2.0), _mcr(3.0), _mcr(4.0) };
    QLComplex m2data[9] = { _mcr(1.0), _mcr(2.0), _mcr(3.0), _mcr(4.0), _mcr(5.0), _mcr(6.0), _mcr(7.0), _mcr(8.0), _mcr(9.0) };

    QLMatrix m1 = QLMatrix::CopyCreate(2, 2, m1data);
    QLMatrix m2 = QLMatrix::CopyCreate(3, 3, m2data);

    m1.Print("m1");
    m2.Print("m2");

    QLMatrix kron = m1.KroneckerProduct(m2);
    kron.Print("kron");
}

void TestPauliSimulate()
{
    QLMatrix m(8, 8);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");
    PrintPauliDecompsedMatrix(DecomposePauliProducts(m)).Print("dec");

    QLGate ch1 = PauliSimulateGate(m, F(0.5), 40);
    QLGate ch2 = PauliSimulateGateLeapfrog(m, F(0.5), 10, F(0.03));

    //TArray<BYTE> remap = ch.m_lstQubits;
    //remap.Pop();
    //remap.InsertAt(0, static_cast<BYTE>(remap.Num()));

    //ch.ApplyOnQubits(remap);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch1;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    QLMatrix p05 = output.m_OutputMatrix.GetBlock(0, 8, 0, 8);
    QLMatrix m05 = output.m_OutputMatrix.GetBlock(8, 8, 8, 8);

    appGeneral(_T("============= TROTTER ============\n"));
    m05.Print("m05");

    n = m;
    m.MatrixIExp(-F(0.5));
    m.Print("m05exp");

    p05.Print("p05");
    n.MatrixIExp(F(0.5));
    n.Print("p05exp");

    p05 = p05 - n;
    Real deltap05 = p05.VectorDot(p05).x;
    m05 = m05 - m;
    Real deltam05 = m05.VectorDot(m05).x;

    appGeneral(_T("============= TROTTER P05 delta = %f M05 delta = %f ============\n\n"), deltap05, deltam05);

    param.m_MasterGate = ch2;
    sim.Simulate(&param, &output);

    p05 = output.m_OutputMatrix.GetBlock(0, 8, 0, 8);
    m05 = output.m_OutputMatrix.GetBlock(8, 8, 8, 8);

    appGeneral(_T("============= 2ND TROTTER ============\n"));
    m05.Print("m05");

    m.Print("m05exp");

    p05.Print("p05");
    n.Print("p05exp");

    p05 = p05 - n;
    deltap05 = p05.VectorDot(p05).x;
    m05 = m05 - m;
    deltam05 = m05.VectorDot(m05).x;

    appGeneral(_T("============= 2ND TROTTER P05 delta = %f M05 delta = %f ============\n\n"), deltap05, deltam05);
}

void TestPauliSimulateController()
{
    QLMatrix m(4, 4);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");
    QLGate ch2 = PauliSimulateGateLeapfrog(m, F(0.5), 10, F(0.01));
    //QLGate ch2 = PauliSimulateGateLeapfrog(m, F(0.5), 2, F(0.8));
    QLGate cch2 = ch2.CreateControlled();

    //ch2.DebugPrint(4);

    //cch2.DebugPrint(3);

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = cch2;
    param.m_bPrint = TRUE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    n = m;
    m.MatrixIExp(-F(0.5));
    m.Print("m05exp");

    n.MatrixIExp(F(0.5));
    n.Print("p05exp");
}

void TestPhaseEstimate()
{
    BYTE phaseBit = 5;
    Real t = F(1.2);
    QLMatrix m(8, 8);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");

    QLMatrix v, w;
    m.EVD(v, w);
    Real eigenv[8];
    for (INT i = 0; i < 8; ++i)
    {
        eigenv[i] = w.Get(i, 0).x;
    }
    appGeneral(_T("expected eigen values: %f %f %f %f %f %f %f %f \n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3], eigenv[4], eigenv[5], eigenv[6], eigenv[7]);
    for (INT i = 0; i < 8; ++i)
    {
        eigenv[i] = eigenv[i] * t;
        while (eigenv[i] > PI2)
        {
            eigenv[i] = eigenv[i] - PI2;
        }
        while (eigenv[i] < 0)
        {
            eigenv[i] = eigenv[i] + PI2;
        }

        eigenv[i] = eigenv[i] * (1U << phaseBit) / PI2;
    }
    for (INT i = 0; i < 8; ++i)
    {
        for (INT j = i + 1; j < 8; ++j)
        {
            if (eigenv[i] > eigenv[j])
            {
                Real temp = eigenv[i];
                eigenv[i] = eigenv[j];
                eigenv[j] = temp;
            }
        }
    }
    appGeneral(_T("expected peaks: %f %f %f %f %f %f %f %f \n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3], eigenv[4], eigenv[5], eigenv[6], eigenv[7]);

    QLGate ch = QuantumPhaseEstimateWithHSimple(m, t, 20, phaseBit, F(0.02));

    //ch.DebugPrint(1);

    QLSimulatorParametersMeasure param;
    param.m_byQubitCount = 4 + phaseBit;
    param.m_MasterGate = ch;
    param.m_iRepeat = 10000;
    for (BYTE p = 0; p < phaseBit; ++p)
    {
        param.m_lstMeasureBits.AddItem(4 + p);
    }

    QLSimulatorMeasure sim;
    sim.Simulate(&param);
}

void TestConditionalHamitonianEvolution()
{
    Real t = F(1.0);
    QLMatrix m(8, 8);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");

    QLGate ch = ConditionalHamiltonianEvolution(m, 3, t, 5, F(0.02));

    // 0 is (1, 0), the smaller the index, the smaller the range, and the binary string is |54>
    // |0><0| is 4=0, 5=0, act on (1, 0, 0, 0)
    // |1><1| is 4=1, 5=0, act on (0, 1, 0, 0)
    // |2><2| is 4=0, 5=1, act on (0, 0, 1, 0)
    // |3><3| is 4=1, 5=1, act on (0, 0, 0, 1)
    // So it should be
    // exp(+-0)   0       0      0
    //     0    exp(2/4)  0      0
    //     0      0    exp(1/4)  0
    //     0      0       0   exp(3/4)

    QLMatrix correct = QLMatrix::CreateEye(128, 128);
    for (INT i = 1; i < 8; ++i)
    {
        n = m;
        n.MatrixIExp(F(0.125) * i);
        correct.SetBlock(16 * i, 8, 16 * i, 8, n.HostBuffer());
        n.Dagger();
        correct.SetBlock(16 * i + 8, 8, 16 * i + 8, 8, n.HostBuffer());
    }

    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 7;
    param.m_MasterGate = ch;
    param.m_bPrint = FALSE;
    QLSimulatorOutputMatrix output;
    QLSimulatorMatrix sim;
    sim.Simulate(&param, &output);

    correct = correct - output.m_OutputMatrix;
    Real delta = correct.VectorDot(correct).x;
    appGeneral(_T("delta = %f\n"), delta);
}

void TestImprovedQPEInitalState()
{
    QLGate ae = BuildImprovedQPEInitialState(3);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 3;
    param.m_MasterGate = ae;

    /**
    * the tau is |210>
    * so |0> = |000> = (1, 0, 0, 0, 0, 0, 0, 0)
    * so |1> = |001> = (0, 1, 0, 0, 0, 0, 0, 0)
    * so |2> = |010> = (0, 0, 1, 0, 0, 0, 0, 0)
    * ...
    */
    QLSimulatorVector sim;
    sim.Simulate(&param);

    appGeneral(_T("expecting\n"));
    for (INT i = 0; i < 8; ++i)
    {
        appGeneral(_T("%f\n"), _sqrt(2 / F(8.0)) * _sin(PI * (i + F(0.5)) / 8));
    }
}

void TestPhaseEstimateImproved()
{
    BYTE phaseBit = 5;
    Real t = F(1.5);
    QLMatrix m(8, 8);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");

    QLMatrix v, w;
    m.EVD(v, w);
    Real eigenv[8];
    for (INT i = 0; i < 8; ++i)
    {
        eigenv[i] = w.Get(i, 0).x;
    }
    appGeneral(_T("expected eigen values: %f %f %f %f %f %f %f %f \n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3], eigenv[4], eigenv[5], eigenv[6], eigenv[7]);
    for (INT i = 0; i < 8; ++i)
    {
        eigenv[i] = eigenv[i] * t;
        while (eigenv[i] > PI2)
        {
            eigenv[i] = eigenv[i] - PI2;
        }
        while (eigenv[i] < 0)
        {
            eigenv[i] = eigenv[i] + PI2;
        }

        eigenv[i] = eigenv[i] * (1U << phaseBit) / PI2;
    }
    for (INT i = 0; i < 8; ++i)
    {
        for (INT j = i + 1; j < 8; ++j)
        {
            if (eigenv[i] > eigenv[j])
            {
                Real temp = eigenv[i];
                eigenv[i] = eigenv[j];
                eigenv[j] = temp;
            }
        }
    }
    appGeneral(_T("expected peaks: %f %f %f %f %f %f %f %f \n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3], eigenv[4], eigenv[5], eigenv[6], eigenv[7]);

    QLGate ch = QuantumPhaseEstimateWithHImproved(m, t, 20, phaseBit, F(0.02));

    //ch.DebugPrint(2);

    QLSimulatorParametersMeasure param;
    param.m_byQubitCount = 4 + phaseBit;
    param.m_MasterGate = ch;
    param.m_iRepeat = 10000;
    for (BYTE p = 0; p < phaseBit; ++p)
    {
        param.m_lstMeasureBits.AddItem(4 + p);
    }

    QLSimulatorMeasure sim;
    sim.Simulate(&param);
}


void TestHHL()
{
    BYTE phaseBit = 4;
    QLMatrix m(4, 4);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");

    QLMatrix v, w;
    m.EVD(v, w);
    Real eigenv[4];
    for (INT i = 0; i < 4; ++i)
    {
        eigenv[i] = w.Get(i, 0).x;
    }
    appGeneral(_T("expected eigen values: %f %f %f %f\n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3]);
    Real maxAbsEigenValue = F(3.0);
    Real t = PI / maxAbsEigenValue;
    for (INT i = 0; i < 4; ++i)
    {
        eigenv[i] = eigenv[i] * t;
        while (eigenv[i] > PI2)
        {
            eigenv[i] = eigenv[i] - PI2;
        }
        while (eigenv[i] < 0)
        {
            eigenv[i] = eigenv[i] + PI2;
        }

        eigenv[i] = eigenv[i] * (1U << phaseBit) / PI2;
    }
    for (INT i = 0; i < 4; ++i)
    {
        for (INT j = i + 1; j < 4; ++j)
        {
            if (eigenv[i] > eigenv[j])
            {
                Real temp = eigenv[i];
                eigenv[i] = eigenv[j];
                eigenv[j] = temp;
            }
        }
    }
    appGeneral(_T("expected peaks: %f %f %f %f\n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3]);

    QLMatrix y(4, 1);
    y.RandomOne();
    y = y / _sqrt(y.VectorDot(y).x);
    y.Print("y");

    QLMatrix x = m.GELS(y);
    x = x / _sqrt(x.VectorDot(x).x);
    x.Print("x");

    QLGate hhl = HHLGate(m, y.ToVector(), 10, maxAbsEigenValue, phaseBit);
    //hhl.DebugPrint(1);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 4 + phaseBit;
    param.m_MasterGate = hhl;
    param.m_bPrint = FALSE;

    QLSimulatorOutputVector out;
    QLSimulatorVector sim;
    sim.Simulate(&param, &out);

    TArray<BYTE> aftermeasured;
    aftermeasured.AddItem(2);
    aftermeasured.AddItem(2);
    aftermeasured.AddItem(0);
    aftermeasured.AddItem(0);
    aftermeasured.AddItem(0);
    aftermeasured.AddItem(0);
    aftermeasured.AddItem(0);
    aftermeasured.AddItem(1);

    QLMatrix res = ShowStateVectorDetail(out.m_OutputMatrix.HostBuffer(), aftermeasured, FALSE);

    res.Print("res");
}

void TestHHLLarge()
{
    BYTE phaseBit = 6;
    QLMatrix m(8, 8);
    m.RandomOne();
    QLMatrix n = m;
    n.Dagger();
    m = m + n;
    m.Print("m");

    QLMatrix v, w;
    m.EVD(v, w);
    Real eigenv[8];
    for (INT i = 0; i < 8; ++i)
    {
        eigenv[i] = w.Get(i, 0).x;
    }
    appGeneral(_T("expected eigen values: %f %f %f %f %f %f %f %f \n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3], eigenv[4], eigenv[5], eigenv[6], eigenv[7]);
    Real maxAbsEigenValue = F(6.0);
    Real t = PI / maxAbsEigenValue;
    for (INT i = 0; i < 8; ++i)
    {
        eigenv[i] = eigenv[i] * t;
        while (eigenv[i] > PI2)
        {
            eigenv[i] = eigenv[i] - PI2;
        }
        while (eigenv[i] < 0)
        {
            eigenv[i] = eigenv[i] + PI2;
        }

        eigenv[i] = eigenv[i] * (1U << phaseBit) / PI2;
    }
    for (INT i = 0; i < 8; ++i)
    {
        for (INT j = i + 1; j < 8; ++j)
        {
            if (eigenv[i] > eigenv[j])
            {
                Real temp = eigenv[i];
                eigenv[i] = eigenv[j];
                eigenv[j] = temp;
            }
        }
    }
    appGeneral(_T("expected peaks: %f %f %f %f %f %f %f %f \n"), eigenv[0], eigenv[1], eigenv[2], eigenv[3], eigenv[4], eigenv[5], eigenv[6], eigenv[7]);

    QLMatrix y(8, 1);
    y.RandomOne();
    y = y / _sqrt(y.VectorDot(y).x);
    y.Print("y");

    QLMatrix x = m.GELS(y);
    x = x / _sqrt(x.VectorDot(x).x);
    x.Print("x");

    QLGate hhl = HHLGate(m, y.ToVector(), 20, maxAbsEigenValue, phaseBit);
    //hhl.DebugPrint(1);

    QLSimulatorParametersVector param;
    param.m_byQubitCount = 5 + phaseBit;
    param.m_MasterGate = hhl;
    param.m_bPrint = FALSE;

    QLSimulatorOutputVector out;
    QLSimulatorVector sim;
    sim.Simulate(&param, &out);

    TArray<BYTE> aftermeasured;
    aftermeasured.AddItem(2);
    aftermeasured.AddItem(2);
    aftermeasured.AddItem(2);
    aftermeasured.AddItem(0);

    for (BYTE i = 0; i < phaseBit; ++i)
    {
        aftermeasured.AddItem(0);
    }
    
    aftermeasured.AddItem(1);

    QLMatrix res = ShowStateVectorDetail(out.m_OutputMatrix.HostBuffer(), aftermeasured, FALSE);

    res.Print("res");
}

int main()
{
    QLRandomInitializer random;

    //appGeneral(_T("%d"), TestCnRY4());
    TestHHLLarge();
    //TestFRy2();

    //std::vector<QLComplex> l1;
    //l1.push_back(_make_cuComplex(-0.70876551, -0.66743494));
    //l1.push_back(_make_cuComplex(0.03215581, 0.22615937));
    //std::vector<QLComplex> l2;
    //l2.push_back(_make_cuComplex(0.19998138, -0.11040608));
    //l2.push_back(_make_cuComplex(0.95956715, 0.16446531));
    //std::vector<std::vector<QLComplex>> u;
    //u.push_back(l1);
    //u.push_back(l2);
    //QLGate ha = CreateZYZGate(u);
    //ha.Dagger();
    //QLSimulatorParametersMatrix param;
    //param.m_byQubitCount = 1;
    //param.m_MasterGate = ha;
    //QLSimulatorMatrix sim(NULL);
    //sim.Simulate(&param);
}


//=============================================================================
// END OF FILE
//=============================================================================