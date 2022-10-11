//=============================================================================
// FILENAME : UsefulFunctions.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/10/2022 nbale]
//=============================================================================

#include "QuantumFit.h"

void BuildIF(const QLMatrix& xy, const TArray<fitFunction_t>& fitFunctions, UINT size, QLMatrix& mif, QLMatrix& mifd)
{
    assert(fitFunctions.Num() > 0);
    assert(size >= xy.X() + static_cast<UINT>(fitFunctions.Num()));

    QLComplex* fmatrixData = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * xy.X() * fitFunctions.Num()));
    for (UINT i = 0; i < xy.X(); ++i)
    {
        for (UINT j = 0; j < static_cast<UINT>(fitFunctions.Num()); ++j)
        {
            fmatrixData[fitFunctions.Num() * i + j] = (*(fitFunctions[j]))(xy.Get(i, 0));
        }
    }

    QLMatrix fmatrix(xy.X(), static_cast<UINT>(fitFunctions.Num()), fmatrixData);

    //build mifd
    mifd = QLMatrix(size, size);
    mifd.SetBlock(fmatrix.Y(), fmatrix.X(), 0, fmatrix.Y(), fmatrix.HostBuffer());
    fmatrix.Dagger();
    mifd.SetBlock(0, fmatrix.X(), fmatrix.X(), fmatrix.Y(), fmatrix.HostBuffer());

    if (size > fmatrix.X() + fmatrix.Y())
    {
        QLMatrix fillEye = QLMatrix::CreateEye(size - fmatrix.X() - fmatrix.Y(), size - fmatrix.X() - fmatrix.Y());
        mifd.SetBlock(fmatrix.X() + fmatrix.Y(), fillEye.X(), fmatrix.X() + fmatrix.Y(), fillEye.Y(), fillEye.HostBuffer());
    }

    //build mif
    mif = QLMatrix(size, size);
    mif.SetBlock(fmatrix.Y(), fmatrix.X(), 0, fmatrix.Y(), fmatrix.HostBuffer());
    fmatrix.Dagger();
    mif.SetBlock(0, fmatrix.X(), fmatrix.X(), fmatrix.Y(), fmatrix.HostBuffer());

    if (size > fmatrix.X() + fmatrix.Y())
    {
        QLMatrix fillEye = QLMatrix::CreateEye(size - fmatrix.X() - fmatrix.Y(), size - fmatrix.X() - fmatrix.Y());
        mif.SetBlock(fmatrix.X() + fmatrix.Y(), fillEye.X(), fmatrix.X() + fmatrix.Y(), fillEye.Y(), fillEye.HostBuffer());
    }

    //===debug

    //mif.Print("mif");

    //mifd.Print("mifd");

    QLMatrix v, w;
    mif.EVD(v, w);
    Real max = -1.0;
    Real min = -1.0;
    for (UINT i = 0; i < size; ++i)
    {
        if (min < 0 || min > abs(w.Get(i, 0).x))
        {
            min = abs(w.Get(i, 0).x);
        }

        if (max < abs(w.Get(i, 0).x))
        {
            max = abs(w.Get(i, 0).x);
        }
    }
    appGeneral(_T("min max of mif: %f %f\n"), min, max);

    min = -1.0;
    max = -1.0;
    mifd.EVD(v, w);
    for (UINT i = 0; i < size; ++i)
    {
        if (min < 0 || min > abs(w.Get(i, 0).x))
        {
            min = abs(w.Get(i, 0).x);
        }

        if (max < abs(w.Get(i, 0).x))
        {
            max = abs(w.Get(i, 0).x);
        }
    }
    appGeneral(_T("min max of mifd: %f %f\n"), min, max);
}

QLMatrix LoadAMatrix(const CCString& sFileName)
{
    QLMatrix m = ReadCSVR(sFileName);
    return m;
}

void BuildIFUsingA(QLMatrix a, QLMatrix& mif, QLMatrix& mifd)
{
    UINT size = MostSignificantPowerTwo(a.X());
    size = 1U << size;
    mifd = QLMatrix(size, size);
    mifd.SetBlock(a.Y(), a.X(), 0, a.Y(), a.HostBuffer());
    a.Dagger();
    mifd.SetBlock(0, a.X(), a.X(), a.Y(), a.HostBuffer());

    mif = QLMatrix(size, size);
    mif.SetBlock(a.Y(), a.X(), 0, a.Y(), a.HostBuffer());
    a.Dagger();
    mif.SetBlock(0, a.X(), a.X(), a.Y(), a.HostBuffer());

    //mifd.Print("mifd");

    //mif.Print("mif");

    QLMatrix v, w;
    mif.EVD(v, w);
    Real max = -1.0;
    for (UINT i = 0; i < size; ++i)
    {
        if (max < abs(w.Get(i, 0).x))
        {
            max = abs(w.Get(i, 0).x);
        }
    }
    appGeneral(_T("max of mif: %f \n"), max);

    max = -1.0;
    mifd.EVD(v, w);
    for (UINT i = 0; i < size; ++i)
    {
        if (max < abs(w.Get(i, 0).x))
        {
            max = abs(w.Get(i, 0).x);
        }
    }
    appGeneral(_T("max of mifd: %f \n"), max);
}

QLGate BuildQuantumFitCircuitWithAY(const QLMatrix& ay, UINT trotter, BYTE phaseQubitNum, Real maxAbsoluteEigen, Real fTruncate)
{
    QLMatrix a = ay.GetBlock(0, ay.X(), 0, ay.Y() - 1);
    QLMatrix y = ay.GetBlock(0, ay.X(), ay.Y() - 1, 1);

    a.Print("a");
    y.Print("y");

    QLMatrix mif, mifd;
    BuildIFUsingA(a, mif, mifd);

    //mifd.Print("mifd");

    UINT matrixqubit = MostSignificantPowerTwo(a.X());
    UINT matrixSize = 1U << matrixqubit;
    assert(matrixSize == a.X() + a.Y());
    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(matrixqubit + phaseQubitNum + 2));

    //step 1, put y
    TArray<QLComplex> ylist;
    for (UINT i = 0; i < a.Y(); ++i)
    {
        ylist.AddItem(_mcr(0.0001));
    }
    for (UINT i = 0; i < y.X(); ++i)
    {
        ylist.AddItem(y.Get(i, 0));
    }

    //for (INT i = 0; i < ylist.Num(); ++i)
    //{
    //    appGeneral(_T("%f\n"), ylist[i] );
    //}

    QLGate y_gate = AmplitudeEncode(ylist);
    ret.AppendGate(y_gate, y_gate.m_lstQubits);

    QLGate ifd_gate = MatrixPowerGate(mifd, 1, trotter, maxAbsoluteEigen, phaseQubitNum, fTruncate);
    ret.AppendGate(ifd_gate, ifd_gate.m_lstQubits);

    //QLGate if_inverse_gate = MatrixPowerGate(mif, -2, trotter, maxAbsoluteEigen, phaseQubitNum, fTruncate);
    //ret.AppendGate(if_inverse_gate, if_inverse_gate.m_lstQubits);

    return ret;
}

QLGate BuildQuantumFitCircuit(const QLMatrix& xy, const TArray<fitFunction_t>& fitFunctions, UINT trotter, BYTE phaseQubitNum, Real maxAbsoluteEigen, Real fTruncate)
{
    UINT matrixqubit = MostSignificantPowerTwo(xy.X() + static_cast<UINT>(fitFunctions.Num()));
    UINT matrixSize = 1U << matrixqubit;
    QLGate ret;
    ret.AddQubits(static_cast<BYTE>(matrixqubit + phaseQubitNum + 2));
    //build if
    QLMatrix mif, mifd;
    BuildIF(xy, fitFunctions, matrixSize, mif, mifd);

    //step 1, put y
    QLMatrix y = xy.GetBlock(0, xy.X(), 1, 1);
    TArray<Real> ylist;
    for (INT i = 0; i < fitFunctions.Num(); ++i)
    {
        ylist.AddItem(0.0);
    }
    for (INT i = 0; i < xy.X(); ++i)
    {
        ylist.AddItem(xy.Get(i, 1).x);
    }
    for (INT i = 0; i < static_cast<INT>(matrixSize) - xy.X() - fitFunctions.Num(); ++i)
    {
        ylist.AddItem(0.0);
    }
    QLGate y_gate = AmplitudeEncodeReal(ylist);
    ret.AppendGate(y_gate, y_gate.m_lstQubits);
    
    //step 2, apply ifdagger
    //step 3, apply if^-2
    QLGate ifd_gate = MatrixPowerGate(mifd, 1, trotter, maxAbsoluteEigen, phaseQubitNum, fTruncate);
    ret.AppendGate(ifd_gate, ifd_gate.m_lstQubits);

    //QLGate if_inverse_gate = MatrixPowerGate(mif, -2, trotter, maxAbsoluteEigen, phaseQubitNum, fTruncate);
    //ret.AppendGate(if_inverse_gate, if_inverse_gate.m_lstQubits);

    return ret;
}

void SimulateQuantumFit(const QLGate& gate, UINT size, BYTE phaseQubit)
{
    UINT matrixQubit = MostSignificantPowerTwo(size);
    assert(gate.m_lstQubits.Num() == static_cast<INT>(matrixQubit + phaseQubit + 2));

    QLSimulatorParametersVector param;
    param.m_byQubitCount = static_cast<UINT>(matrixQubit + phaseQubit + 2);
    param.m_MasterGate = gate;
    param.m_bPrint = FALSE;

    QLSimulatorOutputVector out;
    QLSimulatorVector sim;
    sim.Simulate(&param, &out);

    TArray<BYTE> aftermeasured;
    for (INT i = 0; i < matrixQubit; ++i)
    {
        aftermeasured.AddItem(2);
    }
    
    aftermeasured.AddItem(0);

    for (BYTE i = 0; i < phaseQubit; ++i)
    {
        aftermeasured.AddItem(0);
    }
    aftermeasured.AddItem(0);

    QLMatrix res = ShowStateVectorDetail(out.m_OutputMatrix.HostBuffer(), aftermeasured, FALSE);

    res.Print("res");
}


//=============================================================================
// END OF FILE
//=============================================================================