//=============================================================================
// FILENAME : FermionSimulation.h
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [07/01/2024 nbale]
//=============================================================================

#ifndef _APP_FERMION_SIMULATION_H_
#define _APP_FERMION_SIMULATION_H_

#include "QuantumLab.h"

__DEFINE_ENUM(EFermionJob,
    EFJ_Simulation,
    EFJ_Measure
)

extern void SaveWaveFunction(const CCString& sFile, Real* real, Real* imag, UINT len);
extern void LoadWaveFunction(const CCString& sFile, Real** real, Real** imag, UINT& len);

extern void Fermion1DSimulation(CParameters& params);
extern void Fermion1DMeasure(CParameters& params);


#endif //#ifndef _APP_FERMION_SIMULATION_H_

//=============================================================================
// END OF FILE
//=============================================================================