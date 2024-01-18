//=============================================================================
// FILENAME : QLSimulatorMeasure.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [05/10/2022 nbale]
//=============================================================================

#ifndef _QLSIMULATORMEASURE_H_
#define _QLSIMULATORMEASURE_H_

__BEGIN_NAMESPACE

class QLAPI QLSimulatorParametersMeasure : public QLSimulatorParameters
{
public:

    QLSimulatorParametersMeasure()
        : QLSimulatorParameters()
        , m_MasterGate()
        , m_byQubitCount(0)
        , m_bPrint(TRUE)
        , m_iRepeat(1000)
        , m_iMeasureUntil(-1)
    {
        
    }

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLGate m_MasterGate;
    BYTE m_byQubitCount;
    UBOOL m_bPrint;

    UINT m_iRepeat;
    INT m_iMeasureUntil;
    TArray<BYTE> m_lstMeasureBits;
};

class QLAPI QLSimulatorOutputMeasure : public QLSimulatorOutput
{
public:

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLMatrix m_OutputMatrix;
    TArray<UINT> m_lstCounts;
    TArray<UINT> m_lstHist;

    UINT m_uiMeasureUntilRes;
    UINT m_uiMeasureUntilCount;
};

class QLAPI QLSimulatorMeasure : public QLSimulator
{

protected:

public:

    QLSimulatorMeasure()
        : QLSimulator()
    {

    }

    void Simulate(QLSimulatorParameters* params, QLSimulatorOutput* output = NULL) const override;

};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATORMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================