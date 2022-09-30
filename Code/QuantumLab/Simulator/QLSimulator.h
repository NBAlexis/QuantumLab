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

class QLAPI QLSimulatorOutput
{
public:
    virtual ~QLSimulatorOutput() {}
};

class QLAPI QLSimulator
{

protected:

public:

    //pOutput will be deleted in the destroy of QLSimulator
    QLSimulator() { }

    virtual ~QLSimulator() { }

    virtual void Simulate(const QLSimulatorParameters* params, QLSimulatorOutput* output = NULL) const = 0;

};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================