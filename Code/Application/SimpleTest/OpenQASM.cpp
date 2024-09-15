//=============================================================================
// FILENAME : OpenQASM.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [28/08/2024 nbale]
//=============================================================================

#include "SimpleTest.h"

void TestZYZ()
{

    /**
    *   cos(t/2)    -i sin(t/2)
    * -i sin(t/2)    cos(t/2)
    */

    /** Note that the sign is same as Qiskit
    *   cos(t/2)   -sin(t/2)
    *   sin(t/2)    cos(t/2)
    */

    /**
    * U = [u11 u12]
    *     [u21 u22]
    *
    * U = exp(i delta) [ exp(i (alpha+beta)/2) cos(theta/2)  exp(i (alpha-beta)/2) sin(theta/2)] = exp(i delta) Rz(-alpha) Ry(-theta) Rz(-beta)
    *                  [-exp(-i(alpha-beta)/2) sin(theta/2)  exp(-i(alpha+beta)/2) cos(theta/2)]
    *
    * for OpenQASM
    * U(a,b,c) where a=-theta, b=-beta, c=-alpha
    */

    //degrees[2] = pi/2
    //degrees[1] = -pi/2
    //degrees[0] = t
    
    //QLGate rz1(EBasicOperation::EBO_RZ, Pi);
    //QLGate cnot1(EBasicOperation::EBO_CX);
    //QLGate ry1(EBasicOperation::EBO_RY, t);
    //QLGate cnot2(EBasicOperation::EBO_CX);
    //QLGate ry2(EBasicOperation::EBO_RY, -t);
    //QLGate rz3(EBasicOperation::EBO_RZ, -Pi);

    //retGate.AppendGate(rz1, target);
    //retGate.AppendGate(cnot1, retGate.m_lstQubits);
    //retGate.AppendGate(ry1, target);
    //retGate.AppendGate(cnot2, retGate.m_lstQubits);
    //retGate.AppendGate(ry2, target);
    //retGate.AppendGate(rz3, target);

    QLGate cx(EBasicOperation::EBO_CX);
    QLGate cy(EBasicOperation::EBO_CY);
    QLGate cz(EBasicOperation::EBO_CZ);
    QLGate cp(EBasicOperation::EBO_CP, F(1.2));

    QLGate crx(EBasicOperation::EBO_CRX, F(3.4));
    QLGate cry(EBasicOperation::EBO_CRY, F(5.6));
    QLGate crz(EBasicOperation::EBO_CRZ, F(7.8));

    QLGate h(EBasicOperation::EBO_H);
    QLGate x(EBasicOperation::EBO_X);
    QLGate y(EBasicOperation::EBO_Y);
    QLGate z(EBasicOperation::EBO_Z);
    QLGate rx(EBasicOperation::EBO_RX, F(1.2));
    QLGate ry(EBasicOperation::EBO_RY, F(3.4));
    QLGate rz(EBasicOperation::EBO_RZ, F(5.6));
    QLGate p(EBasicOperation::EBO_P, F(7.8));
    QLGate ph(EBasicOperation::EBO_Phase, F(1.2));
    QLGate swap = CreateSwapGate();
    QLGate toffoli = CreateToffoliGate();

    QLGate test;
    test.AddQubits(3);
    test.AppendGate(h, 0);
    test.AppendGate(h, 1);
    test.AppendGate(h, 2);
    test.AppendGate(cp, 1, 2);
    test.AppendGate(crx, 0, 1);
    test.AppendGate(rx, 0);
    test.AppendGate(ry, 1);
    test.AppendGate(cx, 1, 0);
    test.AppendGate(cry, 2, 0);
    test.AppendGate(crx, 2, 1);
    test.AppendGate(crz, 0, 1);
    test.AppendGate(cz, 0, 1);
    test.AppendGate(cy, 1, 2);
    test.AppendGate(swap, 1, 2);
    test.AppendGate(toffoli, 2, 1, 0);
    test.AppendGate(x, 2);
    test.AppendGate(y, 1);
    test.AppendGate(z, 0);
    test.AppendGate(cp, 0, 1);
    test.AppendGate(crz, 0, 1);
    test.AppendGate(ph, 2);

    TArray<SBasicOperation> opbuiltin = QLGate::ToBuiltIn(test.GetOperation());

    QLMatrix b1 = QLSimulatorMatrix::ShowMatrix(test);
    //QLMatrix b2 = QLSimulatorMatrix::ShowMatrix(opbuiltin, 3);

    //QLGate test;
    //test.AddQubits(1);
    //test.AppendGate(ry, 0);
    //TArray<SBasicOperation> opbuiltin = QLGate::ToBuiltIn(test.GetOperation());
    //QLMatrix b1 = QLSimulatorMatrix::ShowMatrix(test);
    //QLMatrix b2 = QLSimulatorMatrix::ShowMatrix(opbuiltin, 1);

    b1.Print(_T("b"), TRUE);
    //QLMatrix b2 = b1;
    //b2.Print(_T("b2"));

    //b1.ElementArg();
    //b1.Print(_T("argb"), TRUE);
    //b2.ElementAbs();
    //b2.Print(_T("absb"), TRUE);

    TArray<BYTE> measure;
    //measure.AddItem(0);
    appGeneral(test.ToOpenOASM(measure));
    appGeneral(test.ToQLISP(measure));


}

//=============================================================================
// END OF FILE
//=============================================================================