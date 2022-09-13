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
    QLRandomInitializer random;

    QLMatrix mtr(6, 8);
    mtr.RandomOne();
    printf(mtr.Print().c_str());

    QLMatrix* q;
    QLMatrix* r;
    mtr.QR(&q, &r);
    printf(q->Print().c_str());
    printf(r->Print().c_str());
    appSafeDelete(q);
    appSafeDelete(r);

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