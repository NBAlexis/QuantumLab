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

int main()
{
    QLRandomInitializer random;
    TestCSD();


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