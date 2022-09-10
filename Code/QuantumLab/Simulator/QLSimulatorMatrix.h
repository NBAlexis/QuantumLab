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

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLGate m_MasterGate;

};

class QLAPI QLSimulatorMatrix : public QLSimulator
{

protected:

public:

    QLSimulatorMatrix(std::ostream* pOutput, const std::string& sFloatFormat = "%.6f")
        : QLSimulator(pOutput, sFloatFormat)
    {

    }

    void Simulate(const QLSimulatorParameters* params) const override;

};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATORMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================