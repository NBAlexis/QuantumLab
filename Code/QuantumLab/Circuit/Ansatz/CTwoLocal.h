//=============================================================================
// FILENAME : CTwoLocal.h
// 
// DESCRIPTION:
// see:
// https://quantumcomputing.stackexchange.com/questions/26994/references-for-two-local-forms-in-qiskit
// 
// and
// 
// https://docs.quantum.ibm.com/api/qiskit/qiskit.circuit.library.TwoLocal
// 
// which is from:
// https://www.nature.com/articles/nature23879
// 
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#ifndef _CTWOLOCAL_H_
#define _CTWOLOCAL_H_

__BEGIN_NAMESPACE

enum class ESingleLayer : UINT
{
    RY,
    RYRZ,
};

enum class ELinkLayer : UINT
{
    CX,
    CZ,
    CRX,
};

enum class ELinkStyle : UINT
{
    Full,
    Linear,
    Circular,
    PairWise,
    SCA,
};

class QLAPI CTwoLocal : public CAnsatz
{
public:

    CTwoLocal(BYTE qubits, UINT uiLayerCount, 
        ESingleLayer eSingle = ESingleLayer::RYRZ, 
        ELinkLayer eLinkerLayer = ELinkLayer::CZ,
        ELinkStyle eLinkStyle = ELinkStyle::Linear);

    QLGate BuildState(const TArray<Real>& params) const override;

protected:

    static UINT ParamForSingle(ESingleLayer eSingle)
    {
        switch (eSingle)
        {
        case ESingleLayer::RY:
            return 1;
        case ESingleLayer::RYRZ:
            return 2;
        }
        return 0;
    }
    static UINT ParamForLink(ELinkLayer eLinkerLayer)
    {
        switch (eLinkerLayer)
        {
        case ELinkLayer::CX:
            return 0;
        case ELinkLayer::CZ:
            return 0;
        case ELinkLayer::CRX:
            return 1;
        }
        return 0;
    }

    UINT LinkCount(ELinkStyle eLinkStyle, UINT uiLevel)
    {
        switch (eLinkStyle)
        {
        case ELinkStyle::Full:
            return m_byQubits * (m_byQubits - 1) / 2;
        case ELinkStyle::Linear:
            return m_byQubits - 1;
        case ELinkStyle::Circular:
            return m_byQubits;
        case ELinkStyle::PairWise:
            if (1 & m_byQubits)
            {
                //odd, for example, 3
                //0-1 1-2
                return (m_byQubits - 1) / 2;
            }
            else
            {
                //even, for example, 4
                if (0 & uiLevel)
                {
                    //0-1 2-3
                    return m_byQubits / 2;
                }
                else
                {
                    //1-2
                    return (m_byQubits - 2) / 2;
                }
            }
        case ELinkStyle::SCA:
            return m_byQubits;
        }
        return 0;
    }

    void AddSingleLayer(const TArray<Real>& params, UINT& paramIndex, QLGate& gate, BYTE q) const;
    void AddLinkLayer(const TArray<Real>& params, UINT& paramIndex, QLGate& gate, BYTE q1, BYTE q2) const;

    void AddSingleLayerAll(const TArray<Real>& params, UINT& paramIndex, QLGate& gate) const;
    void AddLinkLayerAll(const TArray<Real>& params, UINT& paramIndex, QLGate& gate, UINT level) const;

    UINT m_uiLayerCount;
    ESingleLayer m_eSingle;
    ELinkLayer m_eLinkerLayer;
    ELinkStyle m_eLinkStyle;
};

__END_NAMESPACE


#endif //#ifndef _CTWOLOCAL_H_

//=============================================================================
// END OF FILE
//=============================================================================