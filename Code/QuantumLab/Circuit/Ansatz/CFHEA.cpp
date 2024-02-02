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

CFHEA::CFHEA(BYTE qubits, UINT uiLayerCount, UBOOL bOnlyReal)
    : CAnsatz(qubits)
    , m_uiLayerCount(uiLayerCount)
    , m_bOnlyReal(bOnlyReal)
{
    m_uiParameterCount = qubits + (2 * qubits - 2) * uiLayerCount;
    if (!m_bOnlyReal)
    {
        m_uiParameterCount = m_uiParameterCount * 2;
    }
    for (UINT i = 0; i < m_uiParameterCount; ++i)
    {
        m_lstParameters.AddItem(RandomF() * PI * 2);
    }
}

QLGate CFHEA::BuildStateReal(const TArray<Real>& params) const
{
    QLGate ret;
    ret.AddQubits(m_byQubits);
    ret.m_sName = _T("FHEA");

    for (BYTE q = 0; q < m_byQubits; ++q)
    {
        QLGate ry(EBasicOperation::EBO_RY, params[q]);
        ret.AppendGate(ry, q);
    }

    UINT qubits = static_cast<UINT>(m_byQubits);
    QLGate cz(EBasicOperation::EBO_CZ);

    for (UINT uiL = 0; uiL < m_uiLayerCount; ++uiL)
    {
        for (UINT q = 0; q < qubits; q += 2)
        {
            if (q + 1 < qubits)
            {
                ret.AppendGate(cz, static_cast<BYTE>(q), static_cast<BYTE>(q + 1));
            }
        }

        for (UINT q = 0; q < qubits; ++q)
        {
            QLGate ry(EBasicOperation::EBO_RY, params[qubits + (2 * qubits - 2) * uiL + q]);
            ret.AppendGate(ry, static_cast<BYTE>(q));
        }

        for (UINT q = 1; q < qubits; q += 2)
        {
            if (q + 1 < qubits)
            {
                ret.AppendGate(cz, static_cast<BYTE>(q), static_cast<BYTE>(q + 1));
            }
        }

        for (UINT q = 1; q < qubits - 1; ++q)
        {
            QLGate ry(EBasicOperation::EBO_RY, params[qubits + (2 * qubits - 2) * uiL + qubits + q - 1]);
            ret.AppendGate(ry, static_cast<BYTE>(q));
        }
    }

    return ret;
}

QLGate CFHEA::BuildStateCmp(const TArray<Real>& params) const
{
    QLGate ret;
    ret.AddQubits(m_byQubits);
    ret.m_sName = _T("FHEA");

    for (BYTE q = 0; q < m_byQubits; ++q)
    {
        QLGate ry(EBasicOperation::EBO_RY, params[2 * q]);
        QLGate rz(EBasicOperation::EBO_RZ, params[2 * q + 1]);
        ret.AppendGate(ry, q);
        ret.AppendGate(rz, q);
    }

    UINT qubits = static_cast<UINT>(m_byQubits);
    QLGate cz(EBasicOperation::EBO_CZ);

    for (UINT uiL = 0; uiL < m_uiLayerCount; ++uiL)
    {
        for (UINT q = 0; q < qubits; q += 2)
        {
            if (q + 1 < qubits)
            {
                ret.AppendGate(cz, static_cast<BYTE>(q), static_cast<BYTE>(q + 1));
            }
        }

        for (UINT q = 0; q < qubits; ++q)
        {
            UINT idx = qubits + (2 * qubits - 2) * uiL + q;
            QLGate ry(EBasicOperation::EBO_RY, params[2 * idx]);
            QLGate rz(EBasicOperation::EBO_RZ, params[2 * idx + 1]);
            ret.AppendGate(ry, static_cast<BYTE>(q));
            ret.AppendGate(rz, static_cast<BYTE>(q));
        }

        for (UINT q = 1; q < qubits; q += 2)
        {
            if (q + 1 < qubits)
            {
                ret.AppendGate(cz, static_cast<BYTE>(q), static_cast<BYTE>(q + 1));
            }
        }

        for (UINT q = 1; q < qubits - 1; ++q)
        {
            UINT idx = qubits + (2 * qubits - 2) * uiL + qubits + q - 1;
            QLGate ry(EBasicOperation::EBO_RY, params[2 * idx]);
            QLGate rz(EBasicOperation::EBO_RZ, params[2 * idx + 1]);
            ret.AppendGate(ry, static_cast<BYTE>(q));
            ret.AppendGate(rz, static_cast<BYTE>(q));
        }
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================