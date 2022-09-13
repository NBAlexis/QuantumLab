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

    /**
    * 1    0
    * 0 exp(it)
    */
    EBO_P,

    /**
    * exp(it)   0
    *    0    exp(it)
    */
    EBO_Phase,

    /**
    *   cos(t/2)    -i sin(t/2)
    * -i sin(t/2)    cos(t/2)
    */
    EBO_RX,

    /** Note that the sign is different from Qiskit
    *   cos(t/2)    sin(t/2)
    *  -sin(t/2)    cos(t/2)
    */
    EBO_RY,

    /**
    *  exp(-it/2)      0
    *  0           exp(it/2)
    */
    EBO_RZ,

    //CNOT, controller first, also is controlled pauli x
    EBO_CX,

    //controlled pauli y
    EBO_CY,

    //controlled pauli z, controlled phase flip
    EBO_CZ,

    //controlled phase shift
    EBO_CP,

    //controlled rotation x
    EBO_CRX,

    //controlled rotation y
    EBO_CRY,

    //controlled rotation z
    EBO_CRZ,



    //controlled collapse
    EBO_CC,

    
};

struct QLAPI SBasicOperation
{
    EBasicOperation m_eOperation;
    TArray<BYTE> m_lstQubits;
    Real m_fClassicalParameter;

    UBOOL operator==(const SBasicOperation& other) const
    {
        return FALSE;
    }
};

struct QLAPI SBasicOperationInGate
{
    EBasicOperation m_eOperation;
    TArray<BYTE> m_lstQubits;

    UBOOL operator==(const SBasicOperationInGate& other) const
    {
        return FALSE;
    }
};

class QLAPI QLGate
{

public:

    QLGate() {}

    QLGate(EBasicOperation eOp, Real fParam = 0.0);

    virtual ~QLGate();

    virtual QLGate Controlled(BYTE controlledQubitCount) const;

    void Dagger();

    /**
    * append a gate as subgate
    * 'lstMappingQubits' are qubits where the gate 'toAppend' acting on
    */
    void AppendGate(QLGate toAppend, const TArray<BYTE>& lstMappingQubits);
    TArray<SBasicOperation> GetOperation(const TArray<BYTE>& lstMappingQubits) const;
    TArray<SBasicOperation> GetOperation() const
    {
        return GetOperation(m_lstQubits);
    }

    static void PerformBasicOperation(const struct Qureg& pReg, const SBasicOperation& op);

    UBOOL m_bBasicOperation;
    UBOOL m_bDagger;
    CCString m_sName;

    TArray<BYTE> m_lstQubits;
    TArray<QLGate> m_lstSubGates;

    Real m_fClassicalParameter;
    TArray<SBasicOperationInGate> m_lstOperations;

    UBOOL operator==(const QLGate& other) const
    {
        return FALSE;
    }

protected:

    void ApplyOnQubits(const TArray<BYTE>& lstMappingQubits);
    TArray<BYTE> ExchangeQubits(const TArray<BYTE>& lstMappingQubits) const;
};


__END_NAMESPACE


#endif //#ifndef _QLGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================