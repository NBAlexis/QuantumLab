//=============================================================================
// FILENAME : Grover.h
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [11/28/2023 nbale]
//=============================================================================

#ifndef _GROVER_H_
#define _GROVER_H_

__BEGIN_NAMESPACE

extern QLGate QLAPI GroverSXGate(BYTE numberOfQubits, UINT elementToFlip, UBOOL bUseAncilla = TRUE);

__END_NAMESPACE


#endif //#ifndef _GROVER_H_

//=============================================================================
// END OF FILE
//=============================================================================