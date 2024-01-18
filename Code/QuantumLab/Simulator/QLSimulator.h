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

typedef void (*OnGateCall)(UINT numberOfGate, Real fProb, const SBasicOperation& sOp);

class QLAPI QLSimulatorParameters
{
public:

    QLSimulatorParameters() : m_pCallBack(NULL) {}

    virtual ~QLSimulatorParameters() {}

    virtual UINT BuildZeroStart(BYTE byQubit, Real* pReal, Real* pImag) const
    {
        UINT uiVectorLength = 1U << byQubit;

        memset(pReal, 0, sizeof(Real) * uiVectorLength);
        memset(pImag, 0, sizeof(Real) * uiVectorLength);
        pReal[0] = F(1.0);

        return uiVectorLength;
    }

    OnGateCall m_pCallBack;
};

class QLAPI QLSimulatorOutput
{
public:

    QLSimulatorOutput() : m_fProbability(F(1.0)) {}
    virtual ~QLSimulatorOutput() {}
    Real m_fProbability;
};

class QLAPI QLSimulator
{

protected:

public:

    //pOutput will be deleted in the destroy of QLSimulator
    QLSimulator() {}

    virtual ~QLSimulator() {}

    virtual void Simulate(QLSimulatorParameters* params, QLSimulatorOutput* output = NULL) const = 0;

};

__END_NAMESPACE


#endif //#ifndef _QLSIMULATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================