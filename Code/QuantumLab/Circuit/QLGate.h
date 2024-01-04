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

enum { _kPrintName = 12 };

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

    /** Note that the sign is same as Qiskit
    *   cos(t/2)   -sin(t/2)
    *   sin(t/2)    cos(t/2)
    */
    EBO_RY,

    /**
    *  exp(-i sigmaz t / 2)
    * 
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

    //=================================
    //Noises added like gate
    EBO_Noise_Damping,
    EBO_Noise_Dephaseing,
    EBO_Noise_Depolarising,
    EBO_Noise_MixPauliX,
    EBO_Noise_MixPauliY,
    EBO_Noise_MixPauliZ,
    EBO_Noise_MixPauliAll,
    EBO_Noise_TwoQubitDephaseing,
    EBO_Noise_TwoQubitDepolarising,


    //composite gates
    EBO_Composite,
    
};

static BYTE ByteSequnce[64] = 
{
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63
};

struct QLAPI SBasicOperation
{
    EBasicOperation m_eOperation;
    TArray<BYTE> m_lstQubits;
    Real m_fClassicalParameter;

    UBOOL operator==(const SBasicOperation& other) const
    {
        return m_eOperation == other.m_eOperation && m_lstQubits == other.m_lstQubits && abs(m_fClassicalParameter - other.m_fClassicalParameter) < _QL_FLT_EPSILON;
    }

    UBOOL IsNoise() const
    {
        return EBasicOperation::EBO_Noise_Damping == m_eOperation
            || EBasicOperation::EBO_Noise_Dephaseing == m_eOperation
            || EBasicOperation::EBO_Noise_Depolarising == m_eOperation
            || EBasicOperation::EBO_Noise_MixPauliX == m_eOperation
            || EBasicOperation::EBO_Noise_MixPauliY == m_eOperation
            || EBasicOperation::EBO_Noise_MixPauliZ == m_eOperation
            || EBasicOperation::EBO_Noise_MixPauliAll == m_eOperation
            || EBasicOperation::EBO_Noise_TwoQubitDephaseing == m_eOperation
            || EBasicOperation::EBO_Noise_TwoQubitDepolarising == m_eOperation;
    }

    UBOOL IsMeasure() const
    {
        return EBasicOperation::EBO_CC == m_eOperation;
    }

    UBOOL IsSingleQubitGate() const
    {
        return !IsMeasure() && !IsNoise() && 1 == m_lstQubits.Num();
    }

    UBOOL IsTwoQubitGate() const
    {
        return !IsMeasure() && !IsNoise() && 2 == m_lstQubits.Num();
    }

    UBOOL IsBasicGate() const
    {
        return IsSingleQubitGate() || (EBasicOperation::EBO_CX == m_eOperation);
    }
};

struct QLAPI SBasicOperationInGate
{
    EBasicOperation m_eOperation;

    UBOOL operator==(const SBasicOperationInGate& other) const
    {
        return m_eOperation == other.m_eOperation;
    }
};

class QLAPI QLGate
{

public:

    QLGate() 
        : m_eOp(EBasicOperation::EBO_Composite) 
        , m_bDagger(FALSE)
    {
    }

    QLGate(EBasicOperation eOp, Real fParam = 0.0);

    QLGate(const QLGate& other);

    virtual ~QLGate();

    /**
    * return a gate with 0 the controller, others the mapping qubits 
    */
    virtual QLGate CreateControlled() const;

    virtual QLGate Controlled(BYTE controlledQubitCount, const TArray<BYTE>& lstMappingQubits) const;

    void Dagger();

    /**
    * append a gate as subgate
    * 'lstMappingQubits' are qubits where the gate 'toAppend' acting on
    */
    void AppendGate(QLGate toAppend, const TArray<BYTE>& lstMappingQubits);

    void AppendGate(QLGate toAppend, BYTE qubit, ...);

    TArray<SBasicOperation> GetOperation(const TArray<BYTE>& lstMappingQubits) const;
    TArray<SBasicOperation> GetOperation() const
    {
        return GetOperation(m_lstQubits);
    }

    BYTE GetQubitCount() const { return static_cast<BYTE>(m_lstQubits.Num()); }

    static Real PerformBasicOperation(const struct Qureg& pReg, const SBasicOperation& op);

    EBasicOperation m_eOp;
    UBOOL m_bDagger;
    CCString m_sName;

    /**
    * For basic gates, 
    * 
    * operation are applied on m_lstQubits[m_lstOperations[i]]
    * where m_lstOperations are 0,1,2,3,4,...
    * and m_lstQubits could be any, for example, 4,1,3,2,...
    * 
    * For sub gates,
    * 
    * m_lstQubits = sub.m_lstQubits (could be part of parent) = sub.sub.m_lstQubits (could be part of parent) = ...
    * operation are applied on sub.sub.m_lstQubits[m_lstOperations[i]] or, just m_lstQubits[m_lstOperations[i]]
    * 
    */
    TArray<BYTE> m_lstQubits;
    TArray<QLGate> m_lstSubGates;

    /**
    * Only used for draw circuit, almost always empty
    */
    TArray<BYTE> m_lstAdditionalQubitsData;

    Real m_fClassicalParameter;
    TArray<SBasicOperationInGate> m_lstOperations;

    UBOOL operator==(const QLGate& other) const
    {
        return FALSE;
    }

    void AddQubits(BYTE add)
    {
        assert(0 == m_lstQubits.Num());
        if (0 == m_lstQubits.Num())
        {
            m_lstQubits.Append(ByteSequnce, add);
            //for (BYTE i = 0; i < add; ++i)
            //{
            //    m_lstQubits.AddItem(i);
            //}
        }
    }

    const QLGate& operator=(const QLGate& other);

    void DebugPrint(INT iDepth) const;

    void ApplyOnQubits(const TArray<BYTE>& lstMappingQubits);

protected:
    
    TArray<BYTE> ExchangeQubits(const TArray<BYTE>& lstMappingQubits) const;

    static void PrintQubits(BYTE qubitCount);
    void PrintOneGate(INT iDepth, BYTE qubitCount) const;
    void PrintOneCompositeGate(BYTE qubitCount) const;
    static void PrintOneOp(const SBasicOperation& op, BYTE qubitCount);
    static void PrintEmptyLine(BYTE qubitCount);
    static void PrintBasicSingle(BYTE target, BYTE qubitCount, const CCString& sName, const CCString& sInner);
    static void PrintBasicTwo(BYTE ctr, BYTE target, BYTE qubitCount, const CCString& sName, const CCString& sInner);
    static void PrintGateName(const CCString& sName);
};


__END_NAMESPACE


#endif //#ifndef _QLGATE_H_

//=============================================================================
// END OF FILE
//=============================================================================