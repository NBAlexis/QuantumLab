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
        : m_MasterGate()
        , m_byQubitCount(0)
        , m_bPrint(TRUE)
    {
        
    }

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLGate m_MasterGate;
    BYTE m_byQubitCount;
    UBOOL m_bPrint;
};

class QLAPI QLSimulatorOutputVector : public QLSimulatorOutput
{
public:

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLMatrix m_OutputMatrix;
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

};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATORVECTOR_H_

//=============================================================================
// END OF FILE
//=============================================================================