//=============================================================================
// FILENAME : QLSimulatorMatrix.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [10/09/2022 nbale]
//=============================================================================

#ifndef _QLSIMULATORMATRIX_H_
#define _QLSIMULATORMATRIX_H_

__BEGIN_NAMESPACE

class QLAPI QLSimulatorParametersMatrix : public QLSimulatorParameters
{
public:

    QLSimulatorParametersMatrix()
        : QLSimulatorParameters()
        , m_MasterGate()
        , m_byQubitCount(0)
        , m_bPrint(TRUE)
    {

    }

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLGate m_MasterGate;
    BYTE m_byQubitCount;
    UBOOL m_bPrint;
    TArray<SBasicOperation> m_op;

};

class QLAPI QLSimulatorOutputMatrix : public QLSimulatorOutput
{
public:

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLMatrix m_OutputMatrix;
};

class QLAPI QLSimulatorMatrix : public QLSimulator
{

protected:

public:

    QLSimulatorMatrix()
        : QLSimulator()
    {

    }

    void Simulate(QLSimulatorParameters* params, QLSimulatorOutput* output = NULL) const override;

    static QLMatrix ShowMatrix(const QLGate& gate);
    static QLMatrix ShowMatrix(const TArray<SBasicOperation>& op, BYTE qubit);

};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATORMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================