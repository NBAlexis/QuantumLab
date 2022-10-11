//=============================================================================
// FILENAME : QuantumFit.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/10/2022 nbale]
//=============================================================================

#include "QuantumFit.h"
#include <random>

enum { _ktestdataCount = 14 };

void SaveTempFitData()
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::uniform_real_distribution<double> distributiony(-0.1, 0.1);

    TArray<Real> xs;
    TArray<Real> ys;
    //a 512 fit as an example
    for (INT i = 0; i < _ktestdataCount; ++i)
    {
        Real x = i * 0.01;
        Real y = 0.3 + 0.5 * x;
        xs.AddItem(x);
        ys.AddItem(y);
    }

    for (INT i = 0; i < _ktestdataCount; ++i)
    {
        for (INT j = i + 1; j < _ktestdataCount; ++j)
        {
            if (xs[i] > xs[j])
            {
                Real temp = xs[i];
                xs[i] = xs[j];
                xs[j] = temp;

                temp = ys[i];
                ys[i] = ys[j];
                ys[j] = temp;
            }
        }
    }

    TArray<QLComplex> data;
    for (INT i = 0; i < _ktestdataCount; ++i)
    {
        data.AddItem(_make_cuComplex(xs[i], F(0.0)));
        data.AddItem(_make_cuComplex(ys[i], F(0.0)));
    }

    QLMatrix m = QLMatrix::CopyCreate(_ktestdataCount, 2, data.GetData());

    SaveCSVR(m, _T("testdata1.csv"));
}

int main()
{
    QLRandomInitializer random;

    BYTE phaseQubit = 6;
    //SaveTempFitData(); 
    //QLMatrix xy = ReadCSVR(_T("testdata1.csv"));
    //xy.Print("xy");
    //TArray<fitFunction_t> fitFunctions;
    //fitFunctions.AddItem(one);
    //fitFunctions.AddItem(x);
    //fitFunctions.AddItem(x2);
    //fitFunctions.AddItem(expm);

    //QLGate gate = BuildQuantumFitCircuit(xy, fitFunctions, 20, phaseQubit, 4.0, 0.1);

    //gate.DebugPrint(2);

    //SimulateQuantumFit(gate, xy, phaseQubit);
    //appSetFloatFormat(_T("%.18f"));

    //_I2.Print("i2");
    //_PauliX.Print("x");
    //_PauliY.Print("y");
    //_PauliZ.Print("z");

    //QLMatrix i2x = _I2.KroneckerProduct(_PauliX);
    //i2x.Print("i2x");

    QLMatrix ay = LoadAMatrix("testdata2.csv");
    QLGate gate = BuildQuantumFitCircuitWithAY(ay, 10, phaseQubit, 11.0, 0.1);

    gate.DebugPrint(2);

    SimulateQuantumFit(gate, ay.X(), phaseQubit);

    return 0;
}


//=============================================================================
// END OF FILE
//=============================================================================