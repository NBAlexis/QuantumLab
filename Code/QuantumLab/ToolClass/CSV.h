//=============================================================================
// FILENAME : CSV.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [08/10/2022 nbale]
//=============================================================================

#ifndef _CSV_H_
#define _CSV_H_

__BEGIN_NAMESPACE

extern QLMatrix QLAPI ReadCSV(const CCString& fileName);

extern void QLAPI SaveCSV(const QLMatrix& m, const CCString& fileName);

extern QLMatrix QLAPI ReadCSVR(const CCString& fileName);

extern void QLAPI SaveCSVR(const QLMatrix& m, const CCString& fileName);

__END_NAMESPACE


#endif //#ifndef _CSV_H_

//=============================================================================
// END OF FILE
//=============================================================================