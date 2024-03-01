//=============================================================================
// FILENAME : CTwoLocalAdaptive.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [29/02/2024 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

CTwoLocalAdaptive::CTwoLocalAdaptive(BYTE qubits,
    ESingleLayer eSingle, ELinkLayer eLinkerLayer, ELinkStyle eLinkStyle)
    : CTwoLocal(qubits, 0, eSingle, eLinkerLayer, eLinkStyle)
    , m_uiFixedLayerCount(0)
    , m_uiMaxLayer(100)
{
    if (eLinkerLayer == ELinkLayer::CRX)
    {
        appCrucial(_T("ELinkLayer::CRX not supported for CTwoLocalAdaptive!\n"));
    }

    if (eLinkStyle == ELinkStyle::PairWise || eLinkStyle == ELinkStyle::SCA)
    {
        appCrucial(_T("ELinkStyle::PairWise and ELinkStyle::SCA not supported for CTwoLocalAdaptive!\n"));
    }
}

QLGate CTwoLocalAdaptive::BuildState(const TArray<Real>& params) const
{
    QLGate ret;
    ret.AddQubits(m_byQubits);
    ret.m_sName = _T("TwoLocalA");

    UINT paramIdx = 0;

    AddSingleLayerAll(params, paramIdx, ret);

    //paramIdx = 0;
    for (UINT i = 0; i < m_uiFixedLayerCount; ++i)
    {
        //AddLinkLayerAll(m_lstFixedParameters, paramIdx, ret, i);
        //AddSingleLayerAll(m_lstFixedParameters, paramIdx, ret);

        AddLinkLayerAll(params, paramIdx, ret, i);
        AddSingleLayerAll(params, paramIdx, ret);
    }

    return ret;
}

void CTwoLocalAdaptive::SaveParameters(const CCString& sFileName) const
{
    TArray<Real> lstp = m_lstParameters;
    //lstp.Append(m_lstFixedParameters);
    SaveCSVAR(lstp.GetData(), 1, lstp.Num(), sFileName);
}

void CTwoLocalAdaptive::IncreaseAdaptive()
{
    if (m_uiFixedLayerCount >= m_uiMaxLayer)
    {
        appGeneral(_T("Increasing level, but already reach max layer count: %d\n"), m_uiFixedLayerCount);
        return;
    }

    ++m_uiFixedLayerCount;
    //TArray<Real> lstp = m_lstParameters;
    //lstp.Append(m_lstFixedParameters);
    //m_lstFixedParameters = lstp;
    UINT uiSingle = ParamForSingle(m_eSingle) * m_byQubits;
    m_uiParameterCount += uiSingle;

    for (UINT i = 0; i < uiSingle; ++i)
    {
        m_lstParameters.InsertAt(0, F(0.0));
        //m_lstParameters[i] = F(0.0);
    }
    appGeneral(_T("Increase to level: %d\n"), m_uiFixedLayerCount);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================