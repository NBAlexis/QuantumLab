//=============================================================================
// FILENAME : QLSimulator.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [10/09/2022 nbale]
//=============================================================================

#ifndef _QLSIMULATOR_H_
#define _QLSIMULATOR_H_

__BEGIN_NAMESPACE

enum class ESimulatorType : UINT
{
    //the output is the matrix of the circuit
    EST_MATRIX,

};

class QLAPI QLSimulatorParameters
{
public:
    virtual ~QLSimulatorParameters() {}
};

class QLAPI QLSimulator
{

protected:

public:

    //pOutput will be deleted in the destroy of QLSimulator
    QLSimulator(std::ostream* pOutput, const std::string& sFloatFormat = "%.6f")
        : m_pOutput(pOutput)
        , m_sFloatFormat(sFloatFormat)
    {
        m_pStdOut = new std::ostream(std::cout.rdbuf());
    }
    virtual ~QLSimulator() 
    {
        appSafeDelete(m_pOutput);
        appSafeDelete(m_pStdOut);
    }

    virtual void Simulate(const QLSimulatorParameters* params) const = 0;

    std::ostream* m_pOutput;
    std::ostream* m_pStdOut;
    std::string m_sFloatFormat;
};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================