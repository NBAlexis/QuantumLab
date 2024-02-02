//=============================================================================
// FILENAME : QLSimulatorVector.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [01/10/2022 nbale]
//=============================================================================

#ifndef _QLSIMULATORVECTOR_H_
#define _QLSIMULATORVECTOR_H_

__BEGIN_NAMESPACE

class QLAPI QLSimulatorParametersVector : public QLSimulatorParameters
{
public:

    QLSimulatorParametersVector()
        : QLSimulatorParameters()
        , m_MasterGate()
        , m_byQubitCount(0)
        , m_bPrint(TRUE)
        , m_bHasInitial(FALSE)
        , m_pRealWavefunction(NULL)
        , m_pImageWavefunction(NULL)
    {
        
    }

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLGate m_MasterGate;
    BYTE m_byQubitCount;
    UBOOL m_bPrint;

    UBOOL m_bHasInitial;
    Real* m_pRealWavefunction;
    Real* m_pImageWavefunction;
};

class QLAPI QLSimulatorOutputVector : public QLSimulatorOutput
{
public:
    QLSimulatorOutputVector() 
        : m_bOutputToBuffer(FALSE) 
        , m_pRealBuffer(NULL)
        , m_pImageBuffer(NULL)
    {
    }

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLMatrix m_OutputMatrix;
    UBOOL m_bOutputToBuffer;
    Real* m_pRealBuffer;
    Real* m_pImageBuffer;
};

class QLAPI QLSimulatorVector : public QLSimulator
{

protected:

public:

    QLSimulatorVector()
        : QLSimulator()
    {

    }

    void Simulate(QLSimulatorParameters* params, QLSimulatorOutput* output = NULL) const override;

    static QLMatrix ShowState(const QLGate& gate);
};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATORVECTOR_H_

//=============================================================================
// END OF FILE
//=============================================================================