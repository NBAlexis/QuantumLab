//=============================================================================
// FILENAME : CFHEA.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

CTwoLocal::CTwoLocal(BYTE qubits, UINT uiLayerCount,
    ESingleLayer eSingle, ELinkLayer eLinkerLayer, ELinkStyle eLinkStyle, EAnsatzInitial eInitial)
    : CAnsatz(qubits)
    , m_uiLayerCount(uiLayerCount)
    , m_eSingle(eSingle)
    , m_eLinkerLayer(eLinkerLayer)
    , m_eLinkStyle(eLinkStyle)
{
    m_uiParameterCount = ParamForSingle(eSingle) * qubits;
    for (UINT level = 0; level < m_uiLayerCount; ++level)
    {
        m_uiParameterCount += ParamForSingle(eSingle) * qubits + ParamForLink(eLinkerLayer) * LinkCount(eLinkStyle, level);
    }

    if (EAnsatzInitial::Small == eInitial && ELinkLayer::CRX == eLinkerLayer)
    {
        appGeneral(_T("EAnsatzInitial::Small does not support ELinkLayer::CRX, go back to random\n"));
        eInitial = EAnsatzInitial::Random;
    }

    if (EAnsatzInitial::MBL == eInitial && (ELinkLayer::CRX == eLinkerLayer || ESingleLayer::RY == eSingle || ESingleLayer::RX == eSingle))
    {
        appGeneral(_T("EAnsatzInitial::MBL does not support ELinkLayer::CRX or ESingleLayer::RY or ESingleLayer::RX, go back to random\n"));
        eInitial = EAnsatzInitial::Random;
    }

    switch (eInitial)
    {
    case EAnsatzInitial::Random:
        appDetailed(_T("initial type: Random\n"));
        for (UINT i = 0; i < m_uiParameterCount; ++i)
        {
            m_lstParameters.AddItem(RandomF() * PI * 2);
        }
        break;
    case EAnsatzInitial::Small:
        {
        appDetailed(_T("initial type: Small\n"));
            const Real upper = PI / (static_cast<Real>(qubits) * static_cast<Real>(m_uiLayerCount));
            for (UINT i = 0; i < m_uiParameterCount; ++i)
            {
                m_lstParameters.AddItem(RandomF() * upper);
            }
        }
        break;
    case EAnsatzInitial::MBL:
        appDetailed(_T("initial type: MBL\n"));
        for (UINT i = 0; i < m_uiParameterCount / 2; ++i)
        {
            m_lstParameters.AddItem(RandomF() * F(0.1));
            m_lstParameters.AddItem(RandomF() * PI * 2 - PI);
        }
        break;
    }

}

void CTwoLocal::AddSingleLayer(const TArray<Real>& params, UINT& paramIndex, QLGate& gate, BYTE qubit) const
{
    switch (m_eSingle)
    {
    case ESingleLayer::RY:
        {
            QLGate ry1(EBasicOperation::EBO_RY, params[paramIndex]);
            gate.AppendGate(ry1, qubit);
            ++paramIndex;
        }
        break;
    case ESingleLayer::RX:
        {
            QLGate rx1(EBasicOperation::EBO_RX, params[paramIndex]);
            gate.AppendGate(rx1, qubit);
            ++paramIndex;
        }
        break;
    case ESingleLayer::RYRZ:
        {
            QLGate ry2(EBasicOperation::EBO_RY, params[paramIndex]);
            QLGate rz2(EBasicOperation::EBO_RZ, params[paramIndex + 1]);
            gate.AppendGate(ry2, qubit);
            gate.AppendGate(rz2, qubit);
            paramIndex += 2;
        }
        break;
    case ESingleLayer::RXRZ:
        {
            QLGate rx3(EBasicOperation::EBO_RX, params[paramIndex]);
            QLGate rz3(EBasicOperation::EBO_RZ, params[paramIndex + 1]);
            gate.AppendGate(rx3, qubit);
            gate.AppendGate(rz3, qubit);
            paramIndex += 2;
        }
        break;
    }
}

void CTwoLocal::AddLinkLayer(const TArray<Real>& params, UINT& paramIndex, QLGate& gate, BYTE q1, BYTE q2) const
{
    switch (m_eLinkerLayer)
    {
    case ELinkLayer::CX:
        {
            QLGate cx(EBasicOperation::EBO_CX);
            gate.AppendGate(cx, q1, q2);
        }
        break;
    case ELinkLayer::CZ:
        {
            QLGate cz(EBasicOperation::EBO_CZ);
            gate.AppendGate(cz, q1, q2);
        }
        break;
    case ELinkLayer::CRX:
        {
            QLGate crx(EBasicOperation::EBO_CRX, params[paramIndex]);
            gate.AppendGate(crx, q1, q2);
            ++paramIndex;
        }
        break;
    }
}

void CTwoLocal::AddSingleLayerAll(const TArray<Real>& params, UINT& paramIndex, QLGate& gate) const
{
    for (BYTE q = 0; q < m_byQubits; ++q)
    {
        AddSingleLayer(params, paramIndex, gate, q);
    }
}

void CTwoLocal::AddLinkLayerAll(const TArray<Real>& params, UINT& paramIndex, QLGate& gate, UINT level) const
{
    switch (m_eLinkStyle)
    {
    case ELinkStyle::Full:
        {
            for (BYTE q1 = 0; q1 < m_byQubits; ++q1)
            {
                for (BYTE q2 = q1 + 1; q2 < m_byQubits; ++q2)
                {
                    AddLinkLayer(params, paramIndex, gate, q1, q2);
                }
            }
        }
        break;
    case ELinkStyle::Linear:
        {
            for (BYTE q1 = 0; q1 < m_byQubits - 1; ++q1)
            {
                AddLinkLayer(params, paramIndex, gate, q1, q1 + 1);
            }
        }
        break;
    case ELinkStyle::Circular:
        {
            for (BYTE q1 = 0; q1 < m_byQubits; ++q1)
            {
                if (q1 == m_byQubits - 1)
                {
                    AddLinkLayer(params, paramIndex, gate, q1, 0);
                }
                else
                {
                    AddLinkLayer(params, paramIndex, gate, q1, q1 + 1);
                }
            }
        }
        break;
    case ELinkStyle::PairWise:
        {
            if (1 & level)
            {
                for (BYTE q1 = 1; q1 < m_byQubits - 1; q1 += 2)
                {
                    AddLinkLayer(params, paramIndex, gate, q1, q1 + 1);
                }
            }
            else
            {
                //0-1, 2-3
                for (BYTE q1 = 0; q1 < m_byQubits - 1; q1 += 2)
                {
                    AddLinkLayer(params, paramIndex, gate, q1, q1 + 1);
                }
            }
        }
        break;
    case ELinkStyle::DoublePairWise:
        {
            for (BYTE q1 = 1; q1 < m_byQubits - 1; q1 += 2)
            {
                AddLinkLayer(params, paramIndex, gate, q1, q1 + 1);
            }

            for (BYTE q1 = 0; q1 < m_byQubits - 1; q1 += 2)
            {
                AddLinkLayer(params, paramIndex, gate, q1, q1 + 1);
            }
        }
        break;
    case ELinkStyle::SCA:
        {
            for (BYTE q1 = 0; q1 < m_byQubits; ++q1)
            {
                //level 0:
                //0 - end
                //1 - 0
                //2 - 1
                //level 1:
                //2 - 1
                //0 - end
                //1 - 0
                INT iq1 = -static_cast<INT>(level) + static_cast<INT>(q1);
                INT iq2 = iq1 - 1;
                while (iq1 < 0)
                {
                    iq1 += m_byQubits;
                }
                while (iq2 < 0)
                {
                    iq2 += m_byQubits;
                }
                while (iq1 >= m_byQubits)
                {
                    iq1 -= m_byQubits;
                }
                while (iq2 >= m_byQubits)
                {
                    iq2 -= m_byQubits;
                }

                AddLinkLayer(params, paramIndex, gate, static_cast<BYTE>(iq1), static_cast<BYTE>(iq2));
            }
        }
        break;
    }
}

QLGate CTwoLocal::BuildState(const TArray<Real>& params) const
{
    QLGate ret;
    ret.AddQubits(m_byQubits);
    ret.m_sName = _T("TwoLocal");

    UINT paramIdx = 0;

    AddSingleLayerAll(params, paramIdx, ret);

    for (UINT i = 0; i < m_uiLayerCount; ++i)
    {
        AddLinkLayerAll(params, paramIdx, ret, i);
        AddSingleLayerAll(params, paramIdx, ret);
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================