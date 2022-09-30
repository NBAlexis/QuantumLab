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

void TestCnU()
{
    QLMatrix m(2, 2);
    m.RandomUnitary();

    m.Print("m");

    QLGate ch = CreateCnU(3, m);
    QLSimulatorParametersMatrix param;
    param.m_byQubitCount = 4;
    param.m_MasterGate = ch;
    QLSimulatorMatrix sim;
    sim.Simulate(&param);
    ch.DebugPrint(1);
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
    m.Print(_T("Delta"));
}

int main()
{
    QLRandomInitializer random;
    TestCSDGate();




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