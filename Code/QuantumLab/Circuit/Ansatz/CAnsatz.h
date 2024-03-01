//=============================================================================
// FILENAME : CAnsatz.h
// 
// DESCRIPTION:
// 
//
// REVISION: [dd/mm/yy]
//  [02/02/2024 nbale]
//=============================================================================

#ifndef _CANSATZ_H_
#define _CANSATZ_H_

__BEGIN_NAMESPACE

class QLAPI CAnsatz
{
public:
    
    CAnsatz(BYTE qubits) : m_byQubits(qubits) {}
    virtual ~CAnsatz() {}

    UINT ParameterCount() const
    {
        return m_uiParameterCount;
    }

    BYTE NumberOfQubits() const
    {
        return m_byQubits;
    }

    TArray<Real> GetParameters() const
    {
        return m_lstParameters;
    }

    void SetParameters(const TArray<Real>& param)
    {
        m_lstParameters = param;
    }

    void SetParameterWithOffset(UINT uiParameterCount, Real fOffset)
    {
        if (uiParameterCount >= static_cast<UINT>(m_lstParameters.Num()))
        {
            appCrucial(_T("wrong parameter index!\n"));
            return;
        }
        m_lstParameters[uiParameterCount] = m_lstParameters[uiParameterCount] + fOffset;
    }

    virtual void SaveParameters(const CCString& sFileName) const
    {
        SaveCSVAR(m_lstParameters.GetData(), 1, m_lstParameters.Num(), sFileName);
    }

    virtual QLGate BuildState(const TArray<Real>& params) const = 0;

    QLGate BuildStateWithParam() const
    {
        return BuildState(m_lstParameters);
    }

    QLGate BuildStateAndChangeParam(const TArray<Real>& params)
    {
        m_lstParameters = params;
        return BuildStateWithParam();
    }

    QLGate BuildStateOffset(UINT numberOfParamter, Real fOffset) const
    {
        TArray<Real> params = m_lstParameters;
        params[numberOfParamter] = params[numberOfParamter] + fOffset;
        return BuildState(params);
    }

    virtual void IncreaseAdaptive() {}

protected:

    BYTE m_byQubits;
    UINT m_uiParameterCount;
    TArray<Real> m_lstParameters;
};

__END_NAMESPACE


#endif //#ifndef _CANSATZ_H_

//=============================================================================
// END OF FILE
//=============================================================================