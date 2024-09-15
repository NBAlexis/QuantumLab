//=============================================================================
// FILENAME : CTwoLocalAdaptive.h
// 
// DESCRIPTION:
// see:
// 
//
// REVISION: [dd/mm/yy]
//  [29/02/2024 nbale]
//=============================================================================

#ifndef _CTWOLOCALADAPTIVE_H_
#define _CTWOLOCALADAPTIVE_H_

__BEGIN_NAMESPACE

class QLAPI CTwoLocalAdaptive : public CTwoLocal
{
public:

    CTwoLocalAdaptive(BYTE qubits,
        ESingleLayer eSingle = ESingleLayer::RYRZ, 
        ELinkLayer eLinkerLayer = ELinkLayer::CZ,
        ELinkStyle eLinkStyle = ELinkStyle::Circular,
        EAnsatzInitial eInitial = EAnsatzInitial::MBL);

    QLGate BuildState(const TArray<Real>& params) const override;

    void SaveParameters(const CCString& sFileName) const override;

    void IncreaseAdaptive() override;

    void SetMaxLayer(UINT uiMaxLayer) 
    {
        m_uiMaxLayer = uiMaxLayer;
    }

protected:

    UINT m_uiFixedLayerCount;
    UINT m_uiMaxLayer;
    TArray<Real> m_lstFixedParameters;
};

__END_NAMESPACE


#endif //#ifndef _CTWOLOCALADAPETIVE_H_

//=============================================================================
// END OF FILE
//=============================================================================