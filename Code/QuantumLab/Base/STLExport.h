//=============================================================================
// FILENAME : STLExport.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [13/09/2022 nbale]
//=============================================================================

#ifndef _STLEXPORT_
#define _STLEXPORT_

__BEGIN_NAMESPACE

//fuck it
#if 0

template class QLAPI std::basic_string<TCHAR>;

template class QLAPI std::allocator<BYTE>;
template class QLAPI std::vector<BYTE>;
template class QLAPI std::allocator<class QLGate>;
template class QLAPI std::vector<class QLGate>;
template class QLAPI std::allocator<struct SBasicOperationInGate>;
template class QLAPI std::vector<struct SBasicOperationInGate>;

#endif

__END_NAMESPACE

#endif //#ifndef _STLEXPORT_

//=============================================================================
// END OF FILE
//=============================================================================