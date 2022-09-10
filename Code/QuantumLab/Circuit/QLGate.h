//=============================================================================
// FILENAME : QLGate.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/09/2022 nbale]
//=============================================================================

#ifndef _QLGATE_H_
#define _QLGATE_H_

__BEGIN_NAMESPACE

enum class EBasicOperation : UINT
{
    //hadamard
    EBO_H, 

    //Pauli X
    EBO_X, 

    //Pauli Y
    EBO_Y, 

    //Pauli Z
    EBO_Z, 

    EBO_P,

    EBO_RX,

    EBO_RY,

    EBO_RZ,

    EBO_CN,

    EBO_CP,

    EBO_CRX,

    EBO_CRY,

    EBO_CRZ,

    
};

struct QLAPI SBasicOperation
{
    EBasicOperation m_eOperation;
    std::vector<BYTE> m_lstQubits;
    std::vector<Real> m_lstClassicalParameters;
};

struct QLAPI SClassicalParamters
{
    std::string m_lstParameterNames;
    Real m_lstClassicalParameters;
};

class QLAPI QLGate
{

protected:

public:

    QLGate() {}

    QLGate(EBasicOperation eOp, Real fParam = 0.0);

    virtual ~QLGate();

    QLGate* Controlled(BYTE controlledQubitCount) const;
    QLGate* Dagger() const;
    void AppendGate(const QLGate& toAppend, const std::vector<BYTE>& lstMappingQubits, const std::vector<SClassicalParamters>& addtoinalParameters);
    std::vector<SBasicOperation> GetOperation(const std::vector<BYTE>& lstMappingQubits) const;

    UBOOL m_bBasicOperation;
    std::string m_sName;

    std::vector<BYTE> m_lstQubits;
    std::vector<SClassicalParamters> m_lstClassicalInputs;
    
    std::vector<QLGate*> m_lstSubGates;
    std::vector<SBasicOperation> m_lstOperations;
};


__END_NAMESPACE


#endif //#ifndef _QLGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================