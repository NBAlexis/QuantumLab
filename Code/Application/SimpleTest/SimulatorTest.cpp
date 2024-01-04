//=============================================================================
// FILENAME : SimulatorTest.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [03/01/2024 nbale]
//=============================================================================

#include "SimpleTest.h"

void TestNoiseSimulatorUsingPhaseEstimation()
{
    BYTE phaseBit = 5;
    Real t = F(1.5);
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

    QLGate ch = QuantumPhaseEstimateWithHImproved(m, t, 25, phaseBit, F(0.01));

    //ch.DebugPrint(2);

    QLSimulatorParametersDensityMatrix param;
    param.m_byQubitCount = 3 + phaseBit;
    param.m_MasterGate = ch;
    param.m_iMeasureTimes = 10000;
    param.m_bMeasurePurity = TRUE;
    param.m_bMeasureFidelity = TRUE;
    param.m_bPrint = TRUE;

    param.m_fMixPauliAfterGate = F(0.0000001);
    param.m_fDampingBeforeMeasure = F(0.01);

    for (BYTE p = 0; p < phaseBit; ++p)
    {
        param.m_lstMeasureQubits.AddItem(3 + p);
    }

    QLSimulatorDensityMatrix sim;
    sim.Simulate(&param);
}


//=============================================================================
// END OF FILE
//=============================================================================